#include "simpleBalancingMethod.H"

namespace Foam {

void simpleBalancingMethod::update_state(const DynamicList<chemistryProblem>& problems) {

    auto my_load   = get_my_load(problems);
    auto all_loads = all_gather(my_load);

    auto root = build_tree(all_loads);

    Info << problems.size() << endl;

    sendRecvInfo info;
    info.destinations       = {Pstream::myProcNo()};
    info.sources            = {Pstream::myProcNo()};
    info.number_of_problems = {problems.size()};

    set_state(info);

}


node_ptr simpleBalancingMethod::build_tree(const DynamicList<chemistryLoad>& loads) const{

    double treshold = 0.5 * compute_mean_load(loads);


    auto [big, small] = loadTree::divide(loads, treshold);

    
    std::reverse(big.begin(), big.end());

    auto root = loadTree::new_node(chemistryLoad(-1, -1));

        
    for (const auto& v : big){
        auto child = loadTree::new_node(v);
        loadTree::add_child(root, child);
    }

    for (const auto& v : small) {
        loadTree::add_child(root->children, v, FindCandidate(), ComputeSendValue());
    }

    return root;

}



std::vector<int> simpleBalancingMethod::times_to_problem_counts(
    const std::vector<double>& times, const DynamicList<chemistryProblem>& problems) const {

    std::vector<int> counts;
    /*    double sum = 0.0;

        for (const auto& time : times){

            for (auto it = problems.rbegin(); it != problems.rend(); ++it)
            {
                if (sum > time) {
                    sum = 0.0;
                    counts.push_back(count);
                    count = 0;
                }
                count++;
                sum+=problems[it].cpuTime;
            }
        }

    */
    return counts;
}

chemistryLoad
simpleBalancingMethod::get_my_load(const DynamicList<chemistryProblem>& problems) const {

    auto   lambda = [](double sum, const chemistryProblem& rhs) { return sum + rhs.cpuTime; };
    double sum    = std::accumulate(problems.begin(), problems.end(), 0.0, lambda);
    return chemistryLoad(Pstream::myProcNo(), sum);
}

double simpleBalancingMethod::compute_mean_load(const DynamicList<chemistryLoad>& loads) const {

    auto   lambda = [](double sum, const chemistryLoad& rhs) { return sum + rhs.value; };
    double sum    = std::accumulate(loads.begin(), loads.end(), 0.0, lambda);
    return sum / loads.size();
}

} // namespace Foam