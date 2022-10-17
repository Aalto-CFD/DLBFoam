/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.
    https://github.com/Aalto-CFD/DLBFoam
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    PSRTest

Description
    Validation test for constant pressure perfectly stirred reactor (PSR) model
    with reference data obtained from Cantera-2.5.1-dev

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "zeroDimensionalFvMesh.H"
#include "fluidMulticomponentThermo.H"
#include "basicChemistryModel.H"
#include "multicomponentMixture.H"
#include "chemistrySolver.H"
#include "OFstream.H"
#include "basicSpecieMixture.H"
#include "cellModeller.H"
#include "thermoTypeFunctions.H"

#include <iostream> 
#include <vector> 

// pyJac dependencies
extern "C" {
    #include "chem_utils.h"
};


void print_passed(scalar err)
{
    Info << "\e[1;32m" << "PASSED (error " <<  err << ")\e[0;0m" << endl;
}

void print_failed(scalar err)
{
    Info << "\e[1;31m" << "FAILED (error " <<  err << ")\e[0;0m" << endl;

}

double relative_error(double input, double reference){
    return std::abs(input - reference) / std::abs(reference);
}

bool checkTemperature(double input, double reference, double tolerance)
{

    Info << "\e[0;33m" << "\nValidating the combustion end temperature: T=" << input << "\e[0;0m" << endl;

    double err = relative_error(input, reference);
    if( err < tolerance)
    {
        print_passed(err);
        return true;
    }
    else
    {
        print_failed(err);
        return false;
    }
}

// - enthalpy of formation computed by pyJac J/kg
bool validatePyjacEnthalpy(PtrList<volScalarField>& Y, scalar T, scalar referenceEnthalpy, scalar tolerance)
{

    Info << "\e[0;33m" << "\nValidating enthalpy computation at T=" << T << "\e[0;0m" << endl;

    std::vector<scalar> sp_enth_form(Y.size(), 0.0);
    eval_h(T, sp_enth_form.data());

    scalar h = 0.0;
    for (label i = 0; i < Y.size(); i++) 
    { 
        h += sp_enth_form[i]*Y[i][0];
    }

    scalar err = relative_error(h, referenceEnthalpy);

    if( err < tolerance)
    {
        print_passed(err);
        return true;
    }
    else
    {
        print_failed(err);
        return false;
    }
}

void performanceCheck(scalar elapsed, scalar reference)
{
    Info << "\e[0;33m" << "\nPerformance evaluation:" << "\e[0;0m" << endl;

    int xFactor = 20;
    if( elapsed < reference*xFactor)
    {
        Info << "\e[1;32m" << "PASSED (cpu time " << elapsed << "s)\e[0;0m" << endl;
    }
    else
    {
        Info << "\e[0;31m" << "WARNING, cpu time exceeds our reference by " << xFactor << " times."  << "\e[0;0m" << endl;
        Info << "\e[0;33m" << "You may want to check your LAPACK installation and compiler setups for better performance." << "\e[0;0m" << endl;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    
    argList::noParallel();

    #define CREATE_MESH createZeroDimensionalFvMesh.H
    #define NO_CONTROL

    #include "StopWatch.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"

    // Initialisation parameters
    const word constProp = "pressure";
    const word fractionBasis = "mole"; // "mass";

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Test set 1:
    // - Compare pyjac enthalpy evaluation against reference
    // - For consistency, ensure that OpenFoam enthalpy evaluationshould is close to pyJac values. 
    // - See README for reasons why lower tolerance is expected for OpenFOAM routine comparison.
    {
        scalar p0 = 1367890.0;
        scalar T0 = 1000.0;
        dictionary fractions;
        fractions.set("CH4", 0.5);
        fractions.set("O2", 1.0);
        fractions.set("N2", 3.76);
        double endTime = 0.07;

        #include "initPSR.H"       

        bool Hcheck0a = validatePyjacEnthalpy(Y, 298.15, -256579.71591916375, 1e-7);
        bool Hcheck0b = validatePyjacEnthalpy(Y, 298.15, thermo.hc().ref()[0], 1e-6);

        #include "solveChemistry.H"

        bool Hcheck1a = validatePyjacEnthalpy(Y, thermo.T()[0], 587946.4983283795, 1e-5);
        bool Hcheck1b = validatePyjacEnthalpy(Y, thermo.T()[0], thermo.ha().ref()[0], 1e-6);

        bool Tcheck = checkTemperature(thermo.T()[0], 2660.038226, 1e-6);

        performanceCheck(chemistryTime.getTotalTime(), 0.3);

        if(!Hcheck0a || !Hcheck0b || !Hcheck1a || !Hcheck1b || !Tcheck){return 1;}
    }

    return 0;
}

// ************************************************************************* //
