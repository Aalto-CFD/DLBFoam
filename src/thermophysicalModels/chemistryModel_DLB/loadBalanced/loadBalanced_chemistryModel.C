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

\*---------------------------------------------------------------------------*/

#include "loadBalanced_chemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoType>
Foam::chemistryModels::loadBalanced<ThermoType>::
    loadBalanced(const fluidMulticomponentThermo& thermo)
    :
        chemistryModels::Standard<ThermoType>(thermo),
        skipSpecies_(this->lookupOrDefault("skipSpecies", false)),
        balancer_(this->subOrEmptyDict("loadbalancing")),
        mapper_(this->subOrEmptyDict("refmapping"), this->thermo()),
        cpuTimes_
        (
            IOobject
            (
                thermo.phasePropertyName("cellCpuTimes"),
                this->time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        ),
        resetSkipSpecies_(false),
        skipThreshold_(this->lookupOrDefault("skipThreshold", small)),
        refMap_
        (
            IOobject
            (
                thermo.phasePropertyName("referenceMap"),
                this->time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        ),
        tabulationPtr_(chemistryTabulationMethod::New(*this, *this)),
        tabulation_(*tabulationPtr_)
    {
        if(balancer_.log())
        {
            cpuSolveFile_ = logFile("cpu_solve.out");
            cpuSolveFile_() << "                  time" << tab
                            << "           getProblems" << tab
                            << "           updateState" << tab
                            << "               balance" << tab
                            << "           solveBuffer" << tab
                            << "             unbalance" << tab
                            << "               rank ID" << endl;
        }

    }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ThermoType>
Foam::chemistryModels::loadBalanced<ThermoType>::
    ~loadBalanced()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::chemistryModels::loadBalanced<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    this->zone_.regenerate();

    const label nGlobalZoneCells =
        returnReduce(this->zone_.nCells(), sumOp());
    const bool chemistryActive =
        this->chemistry() && nGlobalZoneCells > 0;
    tabulation_.reset();
    // CPU time analysis
    clockTime timer;
    scalar t_getProblems(0);
    scalar t_updateState(0);
    scalar t_balance(0);
    scalar t_solveBuffer(0);
    scalar t_unbalance(0);

    if (skipSpecies_ && !chemistryActive)
    {
        for(label i = 0; i < this->nSpecie(); i++)
        {
            if (this->thermo().solveSpecie(i))
            {
                const scalar maxY = max(this->Y()[i].oldTime()).value();

                if(maxY < skipThreshold_)
                {
                    this->thermo().setSpecieInactive(i);
                    const_cast<volScalarField&>(this->Y()[i]) == Zero;
                }
            }
        }
    }
    else
    {
        resetSkipSpecies_ = true;
    }

    if (resetSkipSpecies_ && chemistryActive)
    {
        for(label i = 0; i < this->nSpecie(); i++)
        {
            // ensure all species are active
            if (!this->thermo().solveSpecie(i))
            {
                this->thermo().setSpecieActive(i);
            }
        }
        resetSkipSpecies_ = false;
    }
    this->thermo().syncSpeciesActive();

    if(!this->chemistry())
    {
        return great;
    }

    // Set non-zone cells to zero when a zone restriction is active.
    if (!this->zone_.all())
    {
        for(label i = 0; i < this->nSpecie(); i++)
        {
            this->RR(i) = Zero;
        }
    }

    timer.timeIncrement();
    DynamicList<ChemistryProblem> allProblems = getProblems(deltaT);
    t_getProblems = timer.timeIncrement();

    RecvBuffer<ChemistrySolution> incomingSolutions;

    if(balancer_.active())
    {
        timer.timeIncrement();
        balancer_.updateState(allProblems);
        t_updateState = timer.timeIncrement();

        timer.timeIncrement();
        auto guestProblems = balancer_.balance(allProblems);
        auto ownProblems = balancer_.getRemaining(allProblems);
        t_balance = timer.timeIncrement();

        timer.timeIncrement();
        auto ownSolutions = solveList(ownProblems);
        auto guestSolutions = solveBuffer(guestProblems);
        t_solveBuffer = timer.timeIncrement();

        timer.timeIncrement();
        incomingSolutions = balancer_.unbalance(guestSolutions);
        incomingSolutions.append(ownSolutions);
        t_unbalance = timer.timeIncrement();
    }
    else
    {
        timer.timeIncrement();
        incomingSolutions.append(solveList(allProblems));
        t_solveBuffer = timer.timeIncrement();
    }

    if(balancer_.log())
    {
        if(balancer_.active())
        {
            balancer_.printState();
        }
        cpuSolveFile_() << setw(22)
                        << this->time().userTimeValue()<<tab
                        << setw(22) << t_getProblems<<tab
                        << setw(22) << t_updateState<<tab
                        << setw(22) << t_balance<<tab
                        << setw(22) << t_solveBuffer<<tab
                        << setw(22) << t_unbalance<<tab
                        << setw(22) << Pstream::myProcNo()
                        << endl;
    }
    tabulation_.update();

    return updateReactionRates(incomingSolutions);
}


template <class ThermoType>
void Foam::chemistryModels::loadBalanced<ThermoType>::solveSingle
(
    ChemistryProblem& problem, ChemistrySolution& solution
) const
{
    scalar timeLeft = problem.deltaT;
    scalarField Y0 = problem.Y;
    solution.cellid = problem.cellid;

    // Timer begins
    clockTime time;
    time.timeIncrement();
    if(problem.procNo == Pstream::myProcNo())
    {
        // Composition vector (Yi, T, p, deltaT)
        scalarField phiq(this->nEqns() + 1);
        scalarField Rphiq(this->nEqns() + 1);
        solution.retrieved = retrieveProblem(problem, phiq, Rphiq);
        if(!solution.retrieved)
        {
            // Calculate the chemical source terms
            while(timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(
                    problem.pi,
                    problem.Ti,
                    problem.Y,
                    problem.cellid,
                    dt,
                    problem.deltaTChem);
                timeLeft -= dt;
            }
            tabulateProblem(problem, phiq, Rphiq);
            solution.deltaTChem = problem.deltaTChem;
        }
        else
        {
            solution.deltaTChem = this->deltaTChem_[solution.cellid];
        }
    }
    else
    {
        // Calculate the chemical source terms
        while(timeLeft > small)
        {
            scalar dt = timeLeft;
            this->solve(
                problem.pi,
                problem.Ti,
                problem.Y,
                problem.cellid,
                dt,
                problem.deltaTChem);
            timeLeft -= dt;
        }
        solution.deltaTChem = problem.deltaTChem;
    }
    solution.rr = (problem.Y - Y0) * problem.rhoi / problem.deltaT;
    // Timer ends
    solution.cpuTime = time.timeIncrement();
}


template <class ThermoType>
Foam::scalar
Foam::chemistryModels::loadBalanced<ThermoType>::updateReactionRates
(
    const RecvBuffer<ChemistrySolution>& solutions
)
{
    scalar deltaTMin = great;

    for(const auto& array : solutions)
    {
        for(const auto& solution : array)
        {

            updateReactionRate(solution, solution.cellid);

            updateDeltaT(solution, deltaTMin);

            cpuTimes_[solution.cellid] = solution.cpuTime;
        }
    }
    return deltaTMin;
}


template <class ThermoType>
void Foam::chemistryModels::loadBalanced<ThermoType>::updateDeltaT
(
    const ChemistrySolution& solution, scalar& deltaTMin
)
{
    if(!solution.retrieved)
    {
        deltaTMin = min(solution.deltaTChem, deltaTMin);
        this->deltaTChem_[solution.cellid] = min(solution.deltaTChem, this->deltaTChemMax_);
    }
}


template <class ThermoType>
Foam::scalar Foam::chemistryModels::loadBalanced<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template <class ThermoType>
Foam::scalar Foam::chemistryModels::loadBalanced<ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min(
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT);
}


template <class ThermoType>
void Foam::chemistryModels::loadBalanced<ThermoType>::solve
(
    scalar& p,
    scalar& T,
    scalarField& Y,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    chemistryModels::Standard<ThermoType>::solve
    (
        p,
        T,
        Y,
        li,
        deltaT,
        subDeltaT
    );
}


template <class ThermoType>
Foam::RecvBuffer<Foam::ChemistrySolution>
Foam::chemistryModels::loadBalanced<ThermoType>::solveBuffer
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


template <class ThermoType>
Foam::DynamicList<Foam::ChemistrySolution>
Foam::chemistryModels::loadBalanced<ThermoType>::solveList
(
    UList<ChemistryProblem>& problems
) const
{
    DynamicList<ChemistrySolution> solutions(
        problems.size(), ChemistrySolution(this->nSpecie()));

    for(label i = 0; i < problems.size(); ++i)
    {
        solveSingle(problems[i], solutions[i]);
    }
    return solutions;
}



template <class ThermoType>
template<class DeltaTType>
Foam::DynamicList<Foam::ChemistryProblem>
Foam::chemistryModels::loadBalanced<ThermoType>::getProblems
(
    const DeltaTType& deltaT
)
{
    const scalarField& T = this->thermo().T().oldTime();
    const scalarField& p = this->thermo().p().oldTime();
    const volScalarField& rho0vf =
        this->mesh().template lookupObject<volScalarField>
        (
            this->thermo().phasePropertyName("rho")
        ).oldTime();



    DynamicList<ChemistryProblem> solved_problems;
    DynamicList<ChemistryProblem> mapped_problems;

    const label nZoneCells = this->zone_.nCells();

    solved_problems.resize(nZoneCells, ChemistryProblem(this->nSpecie()));
    refMap_ = -1;

    scalarField massFraction(this->nSpecie());

    label counter = 0;
    for(label zci = 0; zci < nZoneCells; zci++)
    {
            const label celli = this->zone_.celli(zci);

            for(label i = 0; i < this->nSpecie(); i++)
            {
                massFraction[i] = this->Y()[i].oldTime()[celli];
            }

            ChemistryProblem problem;
            problem.Y = massFraction;
            problem.Ti = T[celli];
            problem.pi = p[celli];
            problem.rhoi = rho0vf[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.deltaT = deltaT[celli];
            problem.cpuTime = cpuTimes_[celli];
            problem.cellid = celli;
            problem.procNo = Pstream::myProcNo();

            // This check can only be done based on the concentration as the
            // reference temperature is not known
            if (mapper_.shouldMap(massFraction))
            {
                mapped_problems.append(problem);
                refMap_[celli] = 1;
            }

            else
            {
                solved_problems[counter] = problem;
                counter++;
                refMap_[celli] = 2;
            }


    }

    //the real size is set here
    solved_problems.setSize(counter);

    runtime_assert
    (
        solved_problems.size() + mapped_problems.size() == nZoneCells,
        "getProblems fails"
    );

    this->map(mapped_problems, solved_problems);


    return solved_problems;
}


template <class ThermoType>
void Foam::chemistryModels::loadBalanced<ThermoType>::map
(
    DynamicList<ChemistryProblem>& mapped_problems,
    DynamicList<ChemistryProblem>& solved_problems
)
{
    if (mapped_problems.size() > 0)
    {

        ChemistryProblem refProblem = mapped_problems[0];
        scalar refTemperature = refProblem.Ti;

        ChemistrySolution refSolution(this->nSpecie());
        solveSingle(refProblem, refSolution);
        refMap_[refProblem.cellid] = 0;

        for (auto& problem : mapped_problems){

            // Check that the refmap temperature condition is also fullfilled
            if (mapper_.temperatureWithinRange(problem.Ti, refTemperature))
            {
                updateReactionRate(refSolution, problem.cellid);
                cpuTimes_[problem.cellid] = refSolution.cpuTime;
            }
            // Otherwise solve
            else
            {
                solved_problems.append(problem);
                refMap_[problem.cellid] = 2;
            }
        }


    }
}

template <class ThermoType>
void Foam::chemistryModels::loadBalanced<ThermoType>::updateReactionRate
(
    const ChemistrySolution& solution, const label& i
)
{
    for(label j = 0; j < this->nSpecie(); j++)
    {
        this->RR(j)[i] = solution.rr[j];
    }
}


template <class ThermoType>
bool Foam::chemistryModels::loadBalanced<ThermoType>::retrieveProblem
(
    ChemistryProblem& problem, scalarField& phiq, scalarField& Rphiq
) const
{
    for (label i=0; i<this->nEqns()-2; i++)
    {
        phiq[i] = problem.Y[i];
    }
    phiq[phiq.size()-3] = problem.Ti;
    phiq[phiq.size()-2] = problem.pi;
    phiq[phiq.size()-1] = problem.deltaT;

    if(tabulation_.retrieve(phiq,Rphiq))
    {
        // Retrieved solution stored in Rphiq
        scalar Ysum = 0.0;
        for (label i=0; i<this->nEqns()-2; i++)
        {
            problem.Y[i] = Rphiq[i];
            Ysum+=Rphiq[i];
        }
        if(this->nEqns() == this->nSpecie()+1)      // check if you are operating with pyJac solver (nEqns=nSpecie+1)
        {
            problem.Y[this->nSpecie()-1] = 1.0-Ysum;
        }
        return true;
    }
    else
    {
        return false;
    }
}
template <class ThermoType>
void Foam::chemistryModels::loadBalanced<ThermoType>::tabulateProblem
(
    ChemistryProblem& problem, scalarField& phiq, scalarField& Rphiq
) const
{
    if (tabulation_.tabulates())
    {

        for (label i=0; i<this->nEqns()-2; i++)
        {
            Rphiq[i] = problem.Y[i];
        }
        Rphiq[Rphiq.size()-3] = problem.Ti;
        Rphiq[Rphiq.size()-2] = problem.pi;
        Rphiq[Rphiq.size()-1] = problem.deltaT;

        tabulation_.add
        (
            phiq,
            Rphiq,
            this->nEqns()-2,   // TODO: What should be the value here, when we have the pyJac array Nsp-1+T+p???
            problem.cellid,
            problem.deltaT
        );
    }
}
