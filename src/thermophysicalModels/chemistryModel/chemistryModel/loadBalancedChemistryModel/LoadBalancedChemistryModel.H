/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::LoadBalancedChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    LoadBalancedChemistryModelI.H
    LoadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef LoadBalancedChemistryModel_H
#define LoadBalancedChemistryModel_H

#include "ChemistryProblem.H"
#include "ChemistrySolution.H"
#include "LoadBalancer.H"
#include "OFstream.H"
#include "StandardChemistryModel.H"
#include "clockTime.H"
#include "mixtureFractionRefMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class LoadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ReactionThermo, class ThermoType>
class LoadBalancedChemistryModel
    : public StandardChemistryModel<ReactionThermo, ThermoType>
{

private:
    LoadBalancer             balancer_;
    mixtureFractionRefMapper mapper_;
    volScalarField           cpuTimes_;
    autoPtr<OFstream>        cpuSolveFile_;

public:
    //- Runtime type information
    TypeName("loadBalanced");

    // Constructors

    //- Construct from thermo
    LoadBalancedChemistryModel(const ReactionThermo& thermo);

    //- Destructor
    virtual ~LoadBalancedChemistryModel();

    //- Disallow default bitwise copy construction
    LoadBalancedChemistryModel(const LoadBalancedChemistryModel&);

    ///
    ///@brief This is injected with the ODE
    ///
    ///@param p pressure
    ///@param T temperature
    ///@param c concentrations
    ///@param li some cell index, no idea why this is required
    ///@param deltaT flow time step
    ///@param subDeltaT ODE time step
    ///
    virtual void solve(
        scalar&      p,
        scalar&      T,
        scalarField& c,
        const label  li,
        scalar&      deltaT,
        scalar&      subDeltaT) const = 0;

    ///
    ///@brief Solve for generic DeltaTType
    ///
    ///@tparam DeltaTType either field or scalar
    ///@param deltaT flow time step
    ///@return scalar min chemistry time step
    ///
    template <class DeltaTType>
    scalar solve(const DeltaTType& deltaT);

    ///
    ///@brief Specialization for scalar deltaT
    ///
    ///@param deltaT flow time step
    ///@return scalar min chemistry time step
    ///
    virtual scalar solve(const scalar deltaT) override;

    ///
    ///@brief Solve for a field of time steps
    ///
    ///@param deltaT cell specific time steps (required for LTS)
    ///@return scalar min chemistry time step
    ///
    virtual scalar solve(const scalarField& deltaT) override;

    ///
    ///@brief Solves a single problem
    ///
    ///@param problem Input problems
    ///@param solution Output solution
    ///
    void
    solveSingle(ChemistryProblem& problem, ChemistrySolution& solution) const;

    //- Disallow default bitwise assignment
    void operator=(const LoadBalancedChemistryModel&) = delete;

    //- Create and return a loadbal log file of the given name
    inline autoPtr<OFstream> logFile(const word& name) const;

private:
    // TODO: remove
    static mixtureFractionRefMapper createMapper(const ReactionThermo& thermo);

    LoadBalancer createBalancer();

    ///
    ///@brief Solves a list
    ///
    ///@param problems
    ///@return DynamicList<ChemistrySolution> A list of solutions
    ///
    DynamicList<ChemistrySolution>
    solveList(UList<ChemistryProblem>& problems) const;

    ///
    ///@brief Solves the whole buffer coming from the balancer
    ///
    ///@param problems All problems this process whill solve
    ///@return RecvBuffer<ChemistrySolution> A buffer of solutions
    ///
    RecvBuffer<ChemistrySolution>
    solveBuffer(RecvBuffer<ChemistryProblem>& problems) const;

    ///
    ///@brief Get the problems of this process
    ///
    ///@tparam DeltaTType unnecessary template, should be a field
    ///@param deltaT Flow time step
    ///@return DynamicList<ChemistryProblem> List of problems of this process
    ///
    template <class DeltaTType>
    DynamicList<ChemistryProblem> getProblems(const DeltaTType& deltaT);

    ///
    ///@brief Updates the reaction rates for the solution buffer
    ///
    ///@param solutions all solutions of this process
    ///@return scalar minimum chemistry time step
    ///
    scalar updateReactionRates(const RecvBuffer<ChemistrySolution>& solutions);

    ///
    ///@brief Updates the reaction rate for cell i
    ///
    ///@param solution solution of cell i
    ///@param i cell index
    ///
    virtual void
    updateReactionRate(const ChemistrySolution& solution, const label& i);

    ///
    ///@brief The concentration definition
    ///
    ///@param rho density
    ///@param i TODO: WTF?
    ///@param celli TODO: WTF?
    ///@return scalar concentation
    ///
    virtual scalar computeConcentration(
        const scalar& rho, const label& i, const label& celli) const;

    ///
    ///@brief Computes the reaction rate
    ///
    ///@param j TODO: WTF?
    ///@param solution the solution for an individual cell
    ///@return scalar reaction rate
    ///
    virtual scalar computeReactionRate(
        const label& j, const ChemistrySolution& solution) const;

    ///
    ///@brief Computes the mass fractions for a given problem
    ///
    ///@param problem problem to compute the mass fraction for
    ///@return scalarField field of mass fractions of individual species
    ///
    virtual ChemistryProblem
    getMassFraction(const ChemistryProblem& problem) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "LoadBalancedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
