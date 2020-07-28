#pragma once

#include <algorithm>
#include <memory> //std::shared_ptr

#include "IOdictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

#include "algorithms.H"
#include "chemistryLoadBalancingMethod.H"
#include "chemistryProblem.H"
#include "loadTree.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class simpleBalancingMethod : public chemistryLoadBalancingMethod {

private:

    // this value controls the separation of possible sender and receiver ranks
    // a value of 0.5 indicates that load values of higher than 1/0.5 = 2.5 times the
    // mean load on all ranks are to be balanced
    static constexpr scalar HIGH_LOW_LOAD_SEPARATION_TRESHOLD = 0.8;

public:
    simpleBalancingMethod() = default;

    ///
    ///@brief Given a list of all problems, updates the sendRecvInfo member
    ///
    ///@param problems a list of all problems on a rank
    ///
    virtual void update_state(const DynamicList<chemistryProblem>& problems);

protected:

    ///
    ///@brief Computes the sendRecvInfo from a pointer to my node in the load tree
    ///       and the total problem count
    ///
    ///@param my_node pointer to my node in the load tree
    ///@param problems all problems on this rank
    ///@return chemistryLoadBalancingMethod::sendRecvInfo 
    ///
    static chemistryLoadBalancingMethod::sendRecvInfo
    compute_info(const node_ptr& my_node, const DynamicList<chemistryProblem>& problems);


    ///
    ///@brief Builds the load tree from a list of loads on all processes
    ///
    ///@param loads list of loads on all processes 
    ///@return node_ptr pointer to the root of the tree
    ///
    static node_ptr build_tree(const DynamicList<chemistryLoad>& loads);


    ///
    ///@brief Applies a balancing operation to the readily built tree such that 
    ///       each subtree below the root has as close to even mean load 
    ///
    ///@param root of the load tree
    ///
    static void balance_tree(node_ptr& root);


    ///
    ///@brief Applies balancing to the parent assuming that the parent is the sender
    ///       and its children are receivers.
    ///
    ///@param parent root node of the _subtree_
    ///
    static void balance_subtree(node_ptr& parent);



    ///
    ///@brief Computes the load of this process from a given list of problems on the process 
    ///
    ///@param problems list of problems on this process
    ///@return chemistryLoad the load of this process
    ///
    static chemistryLoad compute_my_load(const DynamicList<chemistryProblem>& problems);


    ///
    ///@brief Computes a vector of send counts for a given sender node in the load tree and a list of problems
    ///       on that process
    ///
    ///@param sender the sender node pointer
    ///@param problems list of problems on that process
    ///@return std::vector<int> vector (of min size 1) of integers corresponding to the send counts.
    ///        the first index always corresponds to the amount of problems remaining on the process 
    ///
    static std::vector<int> compute_send_counts(const node_ptr&                      sender,
                                                  const DynamicList<chemistryProblem>& problems);

    ///
    ///@brief Converts of a vector of cpu times to send to a vector of problem counts
    ///
    ///@param times a vector of cpu times to convert
    ///@param problems the list of problems corresponding to the cpu times
    ///@return std::vector<int> a vector of problem counts
    ///
    static std::vector<int> times_to_problem_counts(const std::vector<scalar>&           times,
                                                    const DynamicList<chemistryProblem>& problems);


    ///
    ///@brief Splits the list of loads to two based on the input treshold
    ///
    ///@param in list of loads to split 
    ///@param treshold a treshold separating high values and small values
    ///@return std::pair<std::vector<chemistryLoad>, std::vector<chemistryLoad>>
    ///         a pair of vectors with high values (first) and small values
    ///
    static std::pair<std::vector<chemistryLoad>, std::vector<chemistryLoad>>
    divide(const DynamicList<chemistryLoad>& in, scalar treshold);

private:


    ///
    ///@brief This class is used for estimation of the load balance in the tree build-up
    ///       phase. After the tree has been built, the values are rebalanced.
    ///
    ///
    struct ComputeSendValue {

        //TODO: this value is currently somewhat arbitrary
        scalar operator()(const node_ptr& node, const chemistryLoad& load) const {
            return 0.3 * node->balanced_load.value;
        }
    };


    ///
    ///@brief Finds a largest load value from the vector of candidates.
    ///
    ///
    struct FindCandidate {

        node_ptr operator()(const std::vector<node_ptr>& candidates) const {

            node_ptr correct_parent;
            scalar   largest = 0.0;
            for (const auto candidate : candidates) {

                if (candidate->balanced_load > largest) {
                    largest        = candidate->balanced_load.value;
                    correct_parent = candidate;
                }
            }
            return correct_parent;
        }
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam
