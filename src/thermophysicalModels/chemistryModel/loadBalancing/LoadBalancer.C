/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
\*---------------------------------------------------------------------------*/

#include "LoadBalancer.H"

namespace Foam
{

void
LoadBalancer::updateState(
    const DynamicList<ChemistryProblem>& problems)
{

    auto myLoad   = computeLoad(problems);
    auto allLoads = allGather(myLoad);

    // Info << globalMean << endl;

    auto operations = getOperations(allLoads, myLoad);
    auto info       = operationsToInfo(operations, problems, myLoad);

    setState(info);

    convertNew();
}

LoadBalancerBase::BalancerState
LoadBalancer::operationsToInfo(
    const std::vector<Operation>&        operations,
    const DynamicList<ChemistryProblem>& problems,
    const ChemistryLoad&                 myLoad)
{

    BalancerState info;

    if(isSender(operations, myLoad.rank))
    {
        double              sum = 0.0;
        std::vector<double> times;
        for(const auto& op : operations)
        {
            info.destinations.push_back(op.to);
            sum += op.value;
            times.push_back(op.value);
        }
        info.nProblems = timesToProblemCounts(times, problems);
    }

    // receiver
    else
    {

        for(const auto& op : operations)
        {
            info.sources.push_back(op.from);
        }
        info.nProblems = {problems.size()};
    }

    info.destinations.push_back(Pstream::myProcNo());
    info.sources.push_back(Pstream::myProcNo());

    return info;
}

std::vector<label>
LoadBalancer::timesToProblemCounts(
    const std::vector<scalar>&           times,
    const DynamicList<ChemistryProblem>& problems)
{

    std::vector<int> counts;
    counts.reserve(times.size() + 1);
    auto begin = problems.begin();

    for(const auto& time : times)
    {

        scalar sum(0);
        auto   operation = [&](const ChemistryProblem& problem) {
            sum += problem.cpuTime;
            return sum <= time;
        };
        auto count = count_while(begin, problems.end(), operation);
        begin += count;
        counts.push_back(count);
    }

    int remaining =
        problems.size() - std::accumulate(counts.begin(), counts.end(), 0);
    counts.push_back(remaining);

    return counts;
}

std::vector<LoadBalancer::Operation>
LoadBalancer::getOperations(
    DynamicList<ChemistryLoad>& loads, const ChemistryLoad& myLoad)
{

    double globalMean = getMean(loads);

    std::vector<Operation> operations;

    std::sort(loads.begin(), loads.end());

    auto sender   = loads.end() - 1;
    auto receiver = loads.begin();

    while(sender != receiver)
    {

        double send_value = std::min(
            sender->value - globalMean, globalMean - receiver->value);
        Operation operation{sender->rank, receiver->rank, send_value};
        if(sender->rank == myLoad.rank || receiver->rank == myLoad.rank)
        {
            operations.push_back(operation);
        }
        sender->value -= send_value;
        receiver->value += send_value;

        if(std::abs(sender->value - globalMean) < SMALL)
        {
            sender--;
        }

        else
        {
            receiver++;
        }
    }

    // explicitly filter very small operations
    std::vector<Operation> large;
    for(const auto& op : operations)
    {
        if(op.value > 0.01 * globalMean)
        {
            large.push_back(op);
        }
    }

    runtime_assert(
        !((isSender(operations, myLoad.rank) &&
           isReceiver(operations, myLoad.rank))),
        "Only sender or receiver should be possible.");

    runtime_assert(
        std::abs(getMean(loads) - globalMean) < 1E-7, "Vanishing load");

    return large;

    // return operations;
}

bool
LoadBalancer::isSender(
    const std::vector<Operation>& operations, int rank)
{

    if(operations.size() == 0)
    {
        return false;
    }

    for(const auto& op : operations)
    {
        if(op.from != rank)
        {
            return false;
        }
    }
    return true;
}

bool
LoadBalancer::isReceiver(
    const std::vector<Operation>& operations, int rank)
{

    for(const auto& op : operations)
    {
        if(op.to != rank)
        {
            return false;
        }
    }
    return true;
}

} // namespace Foam
