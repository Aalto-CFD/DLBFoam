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

#include "loadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoType>
Foam::loadBalancedChemistryModel<ThermoType>::
    loadBalancedChemistryModel(const fluidMulticomponentThermo& thermo)
    :
        chemistryModel<ThermoType>(thermo),
        balancer_(createBalancer()),
        mapper_(createMapper(this->thermo())),
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
Foam::loadBalancedChemistryModel<ThermoType>::
    ~loadBalancedChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoType>
Foam::mixtureFractionRefMapper
Foam::loadBalancedChemistryModel<ThermoType>::createMapper
(
    const fluidMulticomponentThermo& thermo
)
{
    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                thermo.phasePropertyName("chemistryProperties"),
                thermo.mesh().time().constant(),
                thermo.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return mixtureFractionRefMapper(chemistryDict_tmp, thermo);
}


template <class ThermoType>
Foam::LoadBalancer
Foam::loadBalancedChemistryModel<ThermoType>::createBalancer()
{
    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                this->thermo().phasePropertyName("chemistryProperties"),
                this->thermo().mesh().time().constant(),
                this->thermo().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return LoadBalancer(chemistryDict_tmp);
}


template <class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::loadBalancedChemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    tabulation_.reset();
    // CPU time analysis
    clockTime timer;
    scalar t_getProblems(0);
    scalar t_updateState(0);
    scalar t_balance(0);
    scalar t_solveBuffer(0);
    scalar t_unbalance(0);

    if(!this->chemistry_)
    {
        return great;
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
void Foam::loadBalancedChemistryModel<ThermoType>::solveSingle
(
    ChemistryProblem& problem, ChemistrySolution& solution
) const
{
    scalar timeLeft = problem.deltaT;
    scalarField c0 = problem.c;
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
                    problem.c,
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
                problem.c,
                problem.cellid,
                dt,
                problem.deltaTChem);
            timeLeft -= dt;
        }
        solution.deltaTChem = problem.deltaTChem;
    }
    solution.rr = (problem.c - c0) * problem.rhoi / problem.deltaT;
    // Timer ends
    solution.cpuTime = time.timeIncrement();
}


template <class ThermoType>
Foam::scalar
Foam::loadBalancedChemistryModel<ThermoType>::updateReactionRates
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
void Foam::loadBalancedChemistryModel<ThermoType>::updateDeltaT
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
Foam::scalar Foam::loadBalancedChemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template <class ThermoType>
Foam::scalar Foam::loadBalancedChemistryModel<ThermoType>::solve
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
Foam::RecvBuffer<Foam::ChemistrySolution>
Foam::loadBalancedChemistryModel<ThermoType>::solveBuffer
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
Foam::loadBalancedChemistryModel<ThermoType>::solveList
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
Foam::loadBalancedChemistryModel<ThermoType>::getProblems
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

    solved_problems.resize(p.size(), ChemistryProblem(this->nSpecie()));

    scalarField massFraction(this->nSpecie());

    label counter = 0;
    forAll(T, celli)
    {
            for(label i = 0; i < this->nSpecie(); i++)
            {
                massFraction[i] = this->Y()[i].oldTime()[celli];
            }

            ChemistryProblem problem;
            problem.c = massFraction;
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

    runtime_assert(solved_problems.size() + mapped_problems.size() == p.size(), "getProblems fails");

    this->map(mapped_problems, solved_problems);


    return solved_problems;
}


template <class ThermoType>
void Foam::loadBalancedChemistryModel<ThermoType>::map
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
void Foam::loadBalancedChemistryModel<ThermoType>::updateReactionRate
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
bool Foam::loadBalancedChemistryModel<ThermoType>::retrieveProblem
(
    ChemistryProblem& problem, scalarField& phiq, scalarField& Rphiq
) const
{
    for (label i=0; i<this->nEqns()-2; i++)
    {
        phiq[i] = problem.c[i];
    }
    phiq[phiq.size()-3] = problem.Ti;
    phiq[phiq.size()-2] = problem.pi;
    phiq[phiq.size()-1] = problem.deltaT;

    if(tabulation_.retrieve(phiq,Rphiq))
    {
        // Retrieved solution stored in Rphiq
        scalar csum = 0.0;
        for (label i=0; i<this->nEqns()-2; i++)
        {
            problem.c[i] = Rphiq[i];
            csum+=Rphiq[i];
        }
        if(this->nEqns() == this->nSpecie()+1)      // check if you are operating with pyJac solver (nEqns=nSpecie+1)
        {
            problem.c[this->nSpecie()-1] = 1.0-csum;
        }
        return true;
    }
    else
    {
        return false;
    }
}
template <class ThermoType>
void Foam::loadBalancedChemistryModel<ThermoType>::tabulateProblem
(
    ChemistryProblem& problem, scalarField& phiq, scalarField& Rphiq
) const
{
    if (tabulation_.tabulates())
    {

        for (label i=0; i<this->nEqns()-2; i++)
        {
            Rphiq[i] = problem.c[i];
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
