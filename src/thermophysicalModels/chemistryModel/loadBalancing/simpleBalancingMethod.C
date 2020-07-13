#include "simpleBalancingMethod.H"

namespace Foam {




std::vector<int>
simpleBalancingMethod::compute_problem_count(const node_ptr&                      sender,
                                             const DynamicList<chemistryProblem>& problems) const {

    std::vector<double> send_times;
    double              total = 0.0;
    for (const auto& child : sender->children) {

        double time_diff = child->value.value - child->org_value.value;
        runtime_assert(time_diff > 0, "A receiver has more load before balancing.");
        send_times.push_back(time_diff);
        total += time_diff;
    }

    double remaining = sender->org_value.value - total;
    send_times.insert(send_times.begin(), remaining);


    return times_to_problem_counts(send_times, problems);


}

chemistryLoadBalancingMethod::sendRecvInfo
simpleBalancingMethod::convert(const node_ptr&                      my_node,
                               const DynamicList<chemistryProblem>& problems) const {

    sendRecvInfo info;

    info.destinations       = {Pstream::myProcNo()};
    info.sources            = {Pstream::myProcNo()};
    info.number_of_problems = {problems.size()};
    // return info;

    // receiver
    if (my_node->parent->value.rank != -1) {

        info.sources.push_back(my_node->parent->value.rank);
        // info.number_of_problems.push_back(1);
        return info;
    } else {

        for (const auto& child : my_node->children) {
            info.destinations.push_back(child->value.rank);
            // info.number_of_problems.push_back(1);
            // info.number_of_problems[0] -= 1;
        }
        info.number_of_problems = compute_problem_count(my_node, problems);
        // auto temp = compute_problem_count(my_node, problems);
    }

    return info;
}

void simpleBalancingMethod::update_state(const DynamicList<chemistryProblem>& problems) {

    auto my_load   = get_my_load(problems);
    auto all_loads = all_gather(my_load);

    auto root = build_tree(all_loads);
    loadTree::print(root);

    auto my_node = loadTree::find(root, my_load.rank);

    runtime_assert(my_node != nullptr, "My node not found from the node tree");

    auto info = convert(my_node, problems);

    set_state(info);
}

node_ptr simpleBalancingMethod::build_tree(const DynamicList<chemistryLoad>& loads) const {

    double treshold = 1.0 * compute_mean_load(loads);

    auto [big, small] = loadTree::divide(loads, treshold);

    std::reverse(big.begin(), big.end());

    auto root = loadTree::new_node(chemistryLoad(-1, -1));

    for (const auto& v : big) {
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