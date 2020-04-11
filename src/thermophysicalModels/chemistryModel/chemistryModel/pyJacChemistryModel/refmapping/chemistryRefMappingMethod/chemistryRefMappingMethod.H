/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::chemistryRefMappingMethod

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryRefMappingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryRefMappingMethod_H
#define chemistryRefMappingMethod_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"



// Mixture fraction headers
#include "atomicWeights.H"
#include "psiReactionThermo.H"
#include "BasicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
class pyJacChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistryRefMappingMethod Declaration
\*---------------------------------------------------------------------------*/

///
///@brief A base class for chemistry reference cell mappings
///
///@tparam CompType NO IDEA
///@tparam ThermoType NO IDEA
///
template<class CompType, class ThermoType>
class chemistryRefMappingMethod
{

public:

    //- Runtime type information
    //TypeName("chemistryRefMappingMethod");

    // Member variables

    int mytest;
    List<scalar> alpha;
    List<scalar> beta_of;

    

    // Constructors

    //- Construct from components
    chemistryRefMappingMethod
    (
        const dictionary& dict
    );


    ///
    ///@brief THIS IS THE BREAD AND BUTTER WHICH SHOULD BE DEFINED IN ALL CLASSES!
    ///       THIS IS NOT PURE AT THIS POINT 
    ///
    // TODO add params
    virtual void applyMapping(pyJacChemistryModel<CompType, ThermoType>& chemistry) const = 0;


    void init_mixture_fraction(basicSpecieMixture& composition){
      
        dictionary mixFracDict(this->coeffsDict_.subDict("mixtureFractionProperties"));

        const wordList& species_ = composition.species();

        List<scalar> alpha_(species_.size(), 0.0);

        forAll(alpha_, i)
        {
            const dictionary& dict = mixFracDict.subDict(species_[i]).subDict("elements");
            scalar a0
            (
                2.0
            *dict.lookupOrDefault<label>("C",0)
            /composition.W(i)
            );
            scalar a1
            (
                0.5
            *dict.lookupOrDefault<label>("H",0)
            /composition.W(i)
            );
            scalar a2
            (
            -1.0
            *dict.lookupOrDefault<label>("O",0)
            /composition.W(i)
            );
            alpha_[i] = a0 + a1 + a2;
        }
        alpha = alpha_;

        //Constant oxidizer and fuel stream variables 
        List<List<scalar> > Yconst( 2, List<scalar>(species_.size(),0.0) );
        forAll(species_, i)
        {
            Yconst[0][i] =
                mixFracDict.subDict
                (
                    "oxidizerMassFractions"
                ).lookupOrDefault<scalar>(species_[i], 0.0);
            Yconst[1][i] =
                mixFracDict.subDict
                (
                    "fuelMassFractions"
                ).lookupOrDefault<scalar>(species_[i], 0.0);
        }

        scalar YoxTot = 0.0;
        scalar YfuTot = 0.0;
        forAll(species_, i)
        {
            YoxTot += Yconst[0][i];
            YfuTot += Yconst[1][i];
        }

        if (mag(1.0 - YoxTot) > SMALL || mag(1.0 - YfuTot) > SMALL)
        {
            FatalErrorIn
            (
                "createMixtureFraction.H :"
            )
            << "oxidizerMassFractions or fuelMassFractions do not sum up to 1.0"
            << abort(FatalError);
        }

        // constructing the beta function at the oxidizer and fuel streams
        List<scalar> beta_of_(2, 0.0);
        forAll(species_, i)
        {
            beta_of_[0] += alpha[i]*Yconst[0][i]; // oxidizer
            beta_of_[1] += alpha[i]*Yconst[1][i]; // fuel
        }
        beta_of = beta_of_;
        // Compute element mass fraction coefficients (only for printing information)
        List<List<scalar> > a(species_.size(), List<scalar>(3,0.0));
        forAll(species_, i)
        {
            const dictionary& dict = mixFracDict.subDict(species_[i]).subDict("elements");
            a[i][0] =
                dict.lookupOrDefault<label>("C",0)
            *atomicWeights["C"]
            /composition.W(i) ;
            a[i][1] =
                dict.lookupOrDefault<label>("H",0)
            *atomicWeights["H"]
            /composition.W(i) ;
            a[i][2] =
                dict.lookupOrDefault<label>("O",0)
            *atomicWeights["O"]
            /composition.W(i) ;
        }

        // Compute element mass fractions (only for printing)
        List<scalar> Zox(species_.size(),0.0);
        List<scalar> Zfu(species_.size(),0.0);
        forAll(Zox, j)
        {
            forAll(species_, i)
            {
                Zox[j] += a[i][j]*Yconst[0][i];
            }
        }
        forAll(Zfu, j)
        {
            forAll(species_, i)
            {
                Zfu[j] += a[i][j]*Yconst[1][i];
            }
        }


        //- stoichiometric mixture fraction
        scalar Z_st =  (0.0 - beta_of_[0])/(beta_of_[1] - beta_of_[0]);




        // Print out some information
        Info<< "    Oxidizer:" << nl
            << "    speciesMassFractions:" << nl;
        forAll(species_, i)
        {
            if (Yconst[0][i] != 0.0)
            {
                Info<< "        Y("<<species_[i]<<") = " << Yconst[0][i] << nl;
            }
        }
        Info<< "    elementMassFractions:" << nl
            << "        Z(C) = " << Zox[0] << nl
            << "        Z(H) = " << Zox[1] << nl
            << "        Z(O) = " << Zox[2] << nl
            << endl;

        Info<< "    Fuel:" << nl
            << "    speciesMassFractions:" << nl;
        forAll(species_, i)
        {
            if (Yconst[1][i] != 0.0)
            {
                Info<< "        Y("<<species_[i]<<") = " << Yconst[1][i] << nl;
            }
        }
        Info<< "    elementMassFractions:" << nl
            << "        Z(C) = " << Zfu[0] << nl
            << "        Z(H) = " << Zfu[1] << nl
            << "        Z(O) = " << Zfu[2] << nl
            << endl;

        // Stoichiometric mixture fraction
        Info<<"    Stoichiomentric mixture fraction Zst = " << Z_st << endl;
    };


    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli)
    {
        scalar beta = 0.0;
        scalar Z;
        forAll(Y, iField)
        {
            const scalarField& Yi = Y[iField];
            beta += alpha[iField]*Yi[celli];
        }
        Z = (beta - beta_of[0])/(beta_of[1] - beta_of[0]);
        if (Z>1e-5)
        {
            return false;
        }
        else
        {
            return true;
        }

    }


    //- Destructor
    virtual ~chemistryRefMappingMethod();


    // Member Functions

    inline bool active()
    {
        return active_;
    }
    


protected:

    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;

private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//TODO: find out when NoRepositiory is turned off
//#ifdef NoRepository
    #include "chemistryRefMappingMethod.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
