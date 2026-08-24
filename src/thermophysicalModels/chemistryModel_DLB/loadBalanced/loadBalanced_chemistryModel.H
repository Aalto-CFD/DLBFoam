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
    Foam::chemistryModels::loadBalanced

Description
    Extends chemistryModels::Standard by adding a load balancing algorithm.

SourceFiles
    loadBalanced_chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalanced_chemistryModel_H
#define loadBalanced_chemistryModel_H

#include "ChemistryProblem.H"
#include "ChemistrySolution.H"
#include "LoadBalancer.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "Standard_chemistryModel.H"
#include "clockTime.H"
#include "mixtureFractionRefMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryModels
{

/*---------------------------------------------------------------------------*\
                        Class loadBalanced Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class loadBalanced
:
    public chemistryModels::Standard<ThermoType>
{

private:
    // Private member data

        // Skip species flag
        Switch  skipSpecies_;

        // Load balancing object
        LoadBalancer balancer_;

        // Reference mapping object
        mixtureFractionRefMapper mapper_;

        // Field containing chemistry CPU time information
        volScalarField cpuTimes_;



        Switch resetSkipSpecies_;

        const scalar skipThreshold_;
        // Field containing reference mapping information
        // 0 -> reference solution
        // 1 -> mapped from reference solution
        // 2 -> solved explicitly
        volScalarField refMap_;

        // A file to output the balancing stats
        autoPtr<OFstream>        cpuSolveFile_;

        //- Tabulation method
        autoPtr<chemistryTabulationMethod> tabulationPtr_;

        //- Tabulation method reference
        chemistryTabulationMethod& tabulation_;


    // Private Member Functions

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
        void updateReactionRate(const ChemistrySolution& solution, const label& i);

        //- Update the reaction rates from a list of solutions
        scalar updateReactionRates(const RecvBuffer<ChemistrySolution>& solutions);

        //- Update the deltaTMin and deltaTChem from a solution
        void updateDeltaT(const ChemistrySolution& solution, scalar& deltaTMin);

        //- Try to retrieve the candidate problem from the ISAT tree
        bool retrieveProblem(ChemistryProblem& problem, scalarField& phiq, scalarField& Rphiq) const;
        void tabulateProblem(ChemistryProblem& problem, scalarField& phiq, scalarField& Rphiq) const;

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
        loadBalanced(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        loadBalanced(const loadBalanced&);


    //- Destructor
    virtual ~loadBalanced();


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

            //- Solve the ODE system.  Forwards to Standard because the solve
            //  overloads above hide it and a using-declaration is not possible:
            //  Standard::solve(const DeltaTType&) is private.
            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const override;

        //- Create and return a log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

    //- Inherit the base class' chemical source term access functions
        using Standard<ThermoType>::RR;
        
	//- Return const access to the chemical source terms
        inline const volScalarField::Internal& RR(const label i) const;

	//- Return non const access to the chemical source terms
        inline volScalarField::Internal& RR(const label i);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const loadBalanced&) = delete;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "loadBalanced_chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "loadBalanced_chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Opts this model into the Standard reaction and reduction registrations of the
// chemistryModel dynamic code template, which tests
// "#if <type>_chemistryModel == Standard_chemistryModel"
#define loadBalanced_chemistryModel 1

#endif

// ************************************************************************* //
