#pragma once

#include <algorithm>
#include <vector>
#include <memory> //std::shared_ptr

#include "IOdictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

#include "algorithms_aalto.H"
#include "chemistryLoadBalancingMethod.H"
#include "chemistryProblem.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class globalBalancingMethod : public chemistryLoadBalancingMethod {

public:
    struct Operation {
        int from, to; //ranks
        double        value;
    };

public:
    globalBalancingMethod() = default;

    ~globalBalancingMethod() {}

    ///
    ///@brief Given a list of all problems, updates the sendRecvInfo member
    ///
    ///@param problems a list of all problems on a rank
    ///
    virtual void update_state(const DynamicList<chemistryProblem>& problems);

protected:




    std::vector<globalBalancingMethod::Operation> get_operations(DynamicList<chemistryLoad>& loads,
                                                                 double global_mean, const chemistryLoad& my_load) const;

    static sendRecvInfo operations_to_info(const std::vector<Operation>&        operations,
                                           const DynamicList<chemistryProblem>& problems,
                                           const chemistryLoad&                 my_load);


    static std::vector<int> times_to_problem_counts(const std::vector<scalar>&           times,
                                                    const DynamicList<chemistryProblem>& problems);

private:


    static bool is_sender(const std::vector<Operation>& operations, int rank);
    static bool is_receiver(const std::vector<Operation>& operations, int rank);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam
