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

Class
    Foam::loadBalancedChemistryModel

Description
    Extends chemistryModel model by adding a load balancing algorithm.

SourceFiles
    loadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalancedChemistryModel_H
#define loadBalancedChemistryModel_H

#include "ChemistryProblem.H"
#include "ChemistrySolution.H"
#include "LoadBalancer.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "chemistryModel.H"
#include "clockTime.H"
#include "mixtureFractionRefMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class loadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class loadBalancedChemistryModel
:
    public chemistryModel<ThermoType>
{

private:
    // Private member data

        // Load balancing object
        LoadBalancer balancer_;

        // Reference mapping object
        mixtureFractionRefMapper mapper_;

        // Field containing chemistry CPU time information
        volScalarField cpuTimes_;

        // Field containing reference mapping information
        // 0 -> reference solution
        // 1 -> mapped from reference solution
        // 2 -> solved explicitly
        volScalarField refMap_;

        // A file to output the balancing stats
        autoPtr<OFstream>        cpuSolveFile_;


    // Private Member Functions

        //- Create a reference mapper object
        static mixtureFractionRefMapper createMapper(const fluidReactionThermo& thermo);

        //- Create a load balancer object
        LoadBalancer createBalancer();

        //- Get the list of problems to be solved
        template<class DeltaTType>
        DynamicList<ChemistryProblem> getProblems(const DeltaTType& deltaT);

        //- Solve a list of chemistry problems and return a list of solutions
        DynamicList<ChemistrySolution>
        solveList(UList<ChemistryProblem>& problems) const;

        //- Solve the problem buffer coming from the balancer
        RecvBuffer<ChemistrySolution>
        solveBuffer(RecvBuffer<ChemistryProblem>& problems) const;

        //- Update the reaction rate of cell i
        virtual void
        updateReactionRate(const ChemistrySolution& solution, const label& i);

        //- Update the reaction rates from a list of solutions
        scalar updateReactionRates(const RecvBuffer<ChemistrySolution>& solutions);

        //- Return the main chemistry variable
        //  (concentration for standard, mass fraction for pyJac)
        virtual scalarField getVariable(const scalarField& concentration, const scalarField& massFraction);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Apply the reference cell mapping
        //TODO: This doesnt belong here and should be in reference cell mapping class.
        //      Implement a generic filter class from which the tabulation methods inherit and make
        //      the mapping base class inherit from that.
        void map(DynamicList<ChemistryProblem>& mapped_problems, DynamicList<ChemistryProblem>& solved_problems);



public:

    //- Runtime type information
    TypeName("loadBalanced");


    // Constructors

        //- Construct from thermo
        loadBalancedChemistryModel(const fluidReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        loadBalancedChemistryModel(const loadBalancedChemistryModel&);


    //- Destructor
    virtual ~loadBalancedChemistryModel();


    // Member Functions

        //- Solve a single ChemistryProblem and put the solution to ChemistrySolution
        void solveSingle(ChemistryProblem& problem, ChemistrySolution& solution) const;

        // Chemistry model functions (overriding functions in
        // chemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT) override;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT) override;

        // ODE functions (overriding abstract functions in ODE.H)
        virtual void solve
        (
            scalar&      p,
            scalar&      T,
            scalarField& c,
            const label  li,
            scalar&      deltaT,
            scalar&      subDeltaT
        ) const = 0;

        //- Create and return a log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const loadBalancedChemistryModel&) = delete;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
#include "loadBalancedChemistryModelI.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "loadBalancedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
