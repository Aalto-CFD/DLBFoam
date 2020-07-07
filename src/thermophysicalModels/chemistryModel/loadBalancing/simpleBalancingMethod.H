#pragma once

#include <memory> //std::shared_ptr
#include <algorithm>

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

#include "loadTree.H"
#include "chemistryProblem.H"
#include "chemistryLoadBalancingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class simpleBalancingMethod : public chemistryLoadBalancingMethod
{

public:


    simpleBalancingMethod() = default;


    virtual void update_state(const DynamicList<chemistryProblem>& problems);



    chemistryLoad get_my_load(const DynamicList<chemistryProblem>& problems) const;

    std::pair<std::vector<chemistryLoad>, std::vector<chemistryLoad>>
    divide(const DynamicList<chemistryLoad>& in, double treshold=0.3) const;

    node_ptr build_tree(const DynamicList<chemistryLoad>& loads) const;


private:


    double compute_mean_load(const DynamicList<chemistryLoad>& loads) const; 



    std::vector<int> times_to_problem_counts(const std::vector<double>& times, const DynamicList<chemistryProblem>& problems) const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam

