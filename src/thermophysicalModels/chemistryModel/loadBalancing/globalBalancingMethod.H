#pragma once

#include <algorithm>
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
        chemistryLoad from, to;
        double        value;
    };

public:
    globalBalancingMethod() = default;

    ///
    ///@brief Given a list of all problems, updates the sendRecvInfo member
    ///
    ///@param problems a list of all problems on a rank
    ///
    virtual void update_state(const DynamicList<chemistryProblem>& problems);

protected:


    double get_balancer_tolerance() const { return m_tolerance; }

    static chemistryLoad get_min(const DynamicList<chemistryLoad>& vec);
    static chemistryLoad get_max(const DynamicList<chemistryLoad>& vec);

    static Operation
    get_operation(const chemistryLoad& sender, const chemistryLoad& receiver, double mean);

    static void apply_operation(DynamicList<chemistryLoad>& loads, const Operation& operation);

    std::vector<globalBalancingMethod::Operation> get_operations(DynamicList<chemistryLoad>& loads,
                                                                 double global_mean, const chemistryLoad& my_load) const;

    static sendRecvInfo operations_to_info(const std::vector<Operation>&        operations,
                                           const DynamicList<chemistryProblem>& problems,
                                           const chemistryLoad&                 my_load);


    static std::vector<int> times_to_problem_counts(const std::vector<scalar>&           times,
                                                    const DynamicList<chemistryProblem>& problems);

private:
    double m_tolerance = 0.05;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam
