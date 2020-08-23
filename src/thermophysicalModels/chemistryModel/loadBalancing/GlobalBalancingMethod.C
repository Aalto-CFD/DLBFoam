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

#include "GlobalBalancingMethod.H"

namespace Foam
{

void
GlobalBalancingMethod::updateState(
    const DynamicList<chemistryProblem>& problems)
{

    auto myLoad   = compute_myLoad(problems);
    auto allLoads = all_gather(myLoad);

    // Info << globalMean << endl;

    auto operations = getOperations(allLoads, myLoad);
    auto info       = operationsToInfo(operations, problems, myLoad);

    set_state(info);
}

chemistryLoadBalancingMethod::sendRecvInfo
GlobalBalancingMethod::operationsToInfo(
    const std::vector<Operation>&        operations,
    const DynamicList<chemistryProblem>& problems,
    const chemistryLoad&                 myLoad)
{

    sendRecvInfo info;

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
        info.number_of_problems = timesToProblemCounts(times, problems);
    }

    // receiver
    else
    {

        for(const auto& op : operations)
        {
            info.sources.push_back(op.from);
        }
        info.number_of_problems = {problems.size()};
    }

    info.destinations.push_back(Pstream::myProcNo());
    info.sources.push_back(Pstream::myProcNo());

    return info;
}

std::vector<label>
GlobalBalancingMethod::timesToProblemCounts(
    const std::vector<scalar>&           times,
    const DynamicList<chemistryProblem>& problems)
{

    std::vector<int> counts;
    counts.reserve(times.size() + 1);
    auto begin = problems.begin();

    for(const auto& time : times)
    {

        scalar sum(0);
        auto   operation = [&](const chemistryProblem& problem) {
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
    runtime_assert(
        std::accumulate(counts.begin(), counts.end(), 0) == problems.size(),
        "Mismatch in the sliced problem count and original problem count.");

    return counts;
}

std::vector<GlobalBalancingMethod::Operation>
GlobalBalancingMethod::getOperations(
    DynamicList<chemistryLoad>& loads, const chemistryLoad& myLoad)
{

    double globalMean = get_mean(loads);

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
        std::abs(get_mean(loads) - globalMean) < 1E-7, "Vanishing load");

    return large;

    // return operations;
}

bool
GlobalBalancingMethod::isSender(
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
GlobalBalancingMethod::isReceiver(
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
