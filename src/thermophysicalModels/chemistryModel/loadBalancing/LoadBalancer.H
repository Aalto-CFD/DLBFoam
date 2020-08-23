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
Class
    Foam::LoadBalancer
Description
    Extends the base class LoadBalancerBase by implementing a
    balancing algorithm which tries to set the global mean load to each rank.
SourceFiles
    LoadBalancer.C
\*---------------------------------------------------------------------------*/

#ifndef LoadBalancer_H
#define LoadBalancer_H

#include "IOdictionary.H"
#include "Switch.H"
#include "algorithms_aalto.H"
#include "LoadBalancerBase.H"
#include "ChemistryProblem.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

#include <algorithm>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class LoadBalancer : public LoadBalancerBase
{

public:
    struct Operation
    {
        int    from, to; // ranks
        double value;
    };

public:
    LoadBalancer() = default;

    virtual ~LoadBalancer() = default;

    ///
    ///@brief Given a list of all problems, updates the BalancerState member
    ///
    ///@param problems a list of all problems on a rank
    ///
    virtual void updateState(const DynamicList<ChemistryProblem>& problems);

protected:
    ///
    ///@brief Get the operations of this process which minimize the deviation
    ///from global mean
    ///
    ///@param loads Loads on all processes
    ///@param myLoad Load on this process
    ///@return std::vector<LoadBalancer::Operation> Operations of this
    ///rank
    ///
    static std::vector<LoadBalancer::Operation> getOperations(
        DynamicList<ChemistryLoad>& loads, const ChemistryLoad& myLoad);

    ///
    ///@brief Converts the operations to send and receive info for the base
    ///class to handle balancing
    ///
    ///@param operations Operations of this process
    ///@param problems Problems of this process
    ///@param myLoad   Load on this process
    ///@return BalancerState An object describing how many problems to send/recv
    ///from different ranks
    ///
    static BalancerState operationsToInfo(
        const std::vector<Operation>&        operations,
        const DynamicList<ChemistryProblem>& problems,
        const ChemistryLoad&                 myLoad);

    ///
    ///@brief Converts the vector of cpu times to number of problems
    ///
    ///@param times Vector of cpu times
    ///@param problems Problems of this process
    ///@return std::vector<label>
    ///
    static std::vector<label> timesToProblemCounts(
        const std::vector<scalar>&           times,
        const DynamicList<ChemistryProblem>& problems);

private:
    ///
    ///@brief Based on the given operations checks if a process is a sender or
    ///receiver.
    ///       Zero operations are marked as a receiver rank.
    ///
    ///@param operations Operations to check.
    ///@param rank       Rank to check.
    ///@return true      If all operations have a from = rank
    ///@return false     If operations are null or any of the operations from !=
    ///rank
    ///
    static bool isSender(const std::vector<Operation>& operations, int rank);

    ///
    ///@brief Based on the given operations checks if a process is a sender or
    ///receiver.
    ///       Zero operations are marked as a receiver rank.
    ///
    ///@param operations Operations to check
    ///@param rank       Rank to check
    ///@return true      If all operations have a to = rank
    ///@return false     Otherwise
    ///
    static bool isReceiver(const std::vector<Operation>& operations, int rank);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam

#endif