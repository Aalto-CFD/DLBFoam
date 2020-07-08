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

    struct ComputeSendValue{

        double operator()(const node_ptr& node, const chemistryLoad& load) const {
            return 0.3 * node->value.value;
        }
    };

    struct FindCandidate{

        node_ptr operator()(const std::vector<node_ptr>& candidates) const{

                        
            node_ptr correct_parent;
            double largest = 0.0;
            for (const auto candidate : candidates) {

                if (candidate->value > largest) {
                    largest        = candidate->value.value; // TODO FIX!
                    correct_parent = candidate;
                }
            }
            return correct_parent;
        }

    };


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam

