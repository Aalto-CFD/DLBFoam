/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM-Aalto library, derived from OpenFOAM.

    https://github.com/blttkgl/OpenFOAM-Aalto

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

\*---------------------------------------------------------------------------*/

#include "LoadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    LoadBalancedChemistryModel(const ReactionThermo& thermo)
    : 
        StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
        balancer_(createBalancer()), 
        mapper_(createMapper(this->thermo())),
        cpuTimes_
        (
            IOobject
            (
                thermo.phasePropertyName("cellCpuTimes"),
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        ),
        refMap_
        (
            IOobject
            (
                thermo.phasePropertyName("referenceMap"),
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    ~LoadBalancedChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::mixtureFractionRefMapper
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::createMapper
(
    const ReactionThermo& thermo
)
{
    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                thermo.phasePropertyName("chemistryProperties"),
                thermo.db().time().constant(),
                thermo.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return mixtureFractionRefMapper(chemistryDict_tmp, thermo.composition());
}


template <class ReactionThermo, class ThermoType>
Foam::LoadBalancer
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::createBalancer()
{

    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                this->thermo().phasePropertyName("chemistryProperties"),
                this->thermo().db().time().constant(),
                this->thermo().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return LoadBalancer(chemistryDict_tmp);
}


template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{

    BasicChemistryModel<ReactionThermo>::correct();

    if(!this->chemistry_)
    {
        return great;
    }

    DynamicList<ChemistryProblem> allProblems = getProblems(deltaT);

    RecvBuffer<ChemistrySolution> incomingSolutions;

    if(balancer_.active())
    {

        balancer_.updateState(allProblems);

        balancer_.printState();

        auto guestProblems = balancer_.balance(allProblems);
        auto ownProblems = balancer_.getRemaining(allProblems);
        auto ownSolutions = solveList(ownProblems);
        auto guestSolutions = solveBuffer(guestProblems);
        incomingSolutions = balancer_.unbalance(guestSolutions);

        incomingSolutions.append(ownSolutions);
    }
    else
    {
        incomingSolutions.append(solveList(allProblems));
    }

    return updateReactionRates(incomingSolutions);
}


template <class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveSingle
(
    ChemistryProblem& problem, ChemistrySolution& solution
) const
{

    scalar timeLeft = problem.deltaT;
    // const scalarList c0       = prob.c;
    scalarField c0 = problem.c;

    // Timer begins
    clockTime time;
    time.timeIncrement();

    // Define a const label to pass as the cell index placeholder
    const label arbitrary = 0;

    // Calculate the chemical source terms
    while(timeLeft > small)
    {
        scalar dt = timeLeft;
        this->solve(
            problem.pi,
            problem.Ti,
            problem.c,
            arbitrary,
            dt,
            problem.deltaTChem);
        timeLeft -= dt;
    }

    solution.c_increment = (problem.c - c0) / problem.deltaT;
    solution.deltaTChem = min(problem.deltaTChem, this->deltaTChemMax_);

    // Timer ends
    solution.cpuTime = time.timeIncrement();

    solution.cellid = problem.cellid;
    solution.rhoi = problem.rhoi;
}


template <class ReactionThermo, class ThermoType>
Foam::scalar
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRates
(
    const RecvBuffer<ChemistrySolution>& solutions
)
{

    scalar deltaTMin = great;

    for(const auto& array : solutions)
    {
        for(const auto& solution : array)
        {

            for(label j = 0; j < this->nSpecie_; j++)
            {
                this->RR_[j][solution.cellid] =
                    solution.c_increment[j] * this->specieThermos_[j].W();
            }

            deltaTMin = min(solution.deltaTChem, deltaTMin);
            
            this->deltaTChem_[solution.cellid] =
                min(solution.deltaTChem, this->deltaTChemMax_);
            
            cpuTimes_[solution.cellid] = solution.cpuTime;
        }
    }

    return deltaTMin;
}


template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min(
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT);
}


template <class ReactionThermo, class ThermoType>
Foam::RecvBuffer<Foam::ChemistrySolution>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveBuffer
(
    RecvBuffer<ChemistryProblem>& problems
) const
{

    // allocate the solutions buffer
    RecvBuffer<ChemistrySolution> solutions;
    
    for(auto& p : problems)
    {
        solutions.append(solveList(p));
    }
    return solutions;
}


template <class ReactionThermo, class ThermoType>
Foam::DynamicList<Foam::ChemistrySolution>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveList
(
    UList<ChemistryProblem>& problems
) const
{

    DynamicList<ChemistrySolution> solutions(
        problems.size(), ChemistrySolution(this->nSpecie_));

    for(label i = 0; i < problems.size(); ++i)
    {
        solveSingle(problems[i], solutions[i]);
    }
    return solutions;
}


template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
Foam::DynamicList<Foam::ChemistryProblem>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getProblems
(
    const DeltaTType& deltaT
)
{

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    bool refCellFound = false;
    DynamicList<ChemistryProblem> chem_problems;
    ChemistrySolution ref_solution(this->nSpecie_);
    mapper_.clear();
    
    forAll(p, celli)
    {
        for(label i = 0; i < this->nSpecie_; i++)
        {
            this->c_[i] = rho[celli] * this->Y_[i][celli] / this->specieThermos_[i].W();
        }

        scalar Ti = T[celli];

        if(Ti > this->Treact())
        {
            // Create and populate the chemisty problem for celli
            ChemistryProblem problem;
            problem.c = this->c_;
            problem.Ti = T[celli];
            problem.pi = p[celli];
            problem.rhoi = rho[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.deltaT = deltaT[celli];
            problem.cpuTime = cpuTimes_[celli];
            problem.cellid = celli;

            // First reference cell is found and solved, following reference
            // cells are mapped from the first reference cell found
            if(mapper_.active() && mapper_.shouldMap(getMassFraction(problem)))
            {
                if(!refCellFound)
                {
                    solveSingle(problem, ref_solution);
                    refCellFound = true;
                    refMap_[celli] = 0;
                }
                else
                {
                    updateReactionRate(ref_solution, celli);
                    refMap_[celli] = 1;
                }
                cpuTimes_[celli] = ref_solution.cpuTime;
                updateReactionRate(ref_solution, celli);

            }
            else
            {
                chem_problems.append(problem);
                refMap_[celli] = 2;
            }
        }
        else
        {
            for(label i = 0; i < this->nSpecie(); i++)
            {
                this->RR_[i][celli] = 0;
            }
        }
    }
    return chem_problems;
}


template <class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRate
(
    const ChemistrySolution& solution, const label& i
)
{
    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][i] = solution.c_increment[j] * this->specieThermos_[j].W();
    }
    this->deltaTChem_[i] = min(solution.deltaTChem, this->deltaTChemMax_);
}


template <class ReactionThermo, class ThermoType>
Foam::ChemistryProblem
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getMassFraction
(
    const ChemistryProblem& problem
) const
{
    ChemistryProblem tmp = problem;
    for(label i = 0; i < this->nSpecie_; i++)
    {
        tmp.c[i] = this->Y_[i][problem.cellid];
    }

    return (tmp);
}
