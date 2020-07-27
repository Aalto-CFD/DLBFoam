#include "simpleBalancingMethod.H"

namespace Foam {




void simpleBalancingMethod::update_state(const DynamicList<chemistryProblem>& problems) {

    auto my_load   = compute_my_load(problems);
    auto all_loads = all_gather(my_load);

    auto root = build_tree(all_loads);

    //loadTree::print(root);

    auto my_node = loadTree::find(root, Pstream::myProcNo());
    runtime_assert(my_node != nullptr, "My node not found from the node tree");
    auto info = compute_info(my_node, problems);

    set_state(info);
}





std::vector<int>
simpleBalancingMethod::compute_send_counts(const node_ptr&                      sender,
                                             const DynamicList<chemistryProblem>& problems) {

    std::vector<double> send_times; 
    send_times.reserve(sender->children.size() + 1);
    
    double              total = 0.0;
    for (const auto& child : sender->children) {

        double time_diff = child->balanced_load.value - child->original_load.value;
        runtime_assert(time_diff > 0, "A receiver has more load before balancing.");
        send_times.push_back(time_diff);
        total += time_diff;
    }

    double remaining = sender->original_load.value - total;
    send_times.insert(send_times.begin(), remaining);

    auto   sum_op = [](double sum, const chemistryProblem& rhs) { return sum + rhs.cpuTime; };
    runtime_assert(
        std::accumulate(send_times.begin(), send_times.end(), 0.0) == 
        std::accumulate(problems.begin(), problems.end(), 0.0, sum_op),
        "Original problem time not matching the split send times"
    );


    return times_to_problem_counts(send_times, problems);
}

chemistryLoadBalancingMethod::sendRecvInfo
simpleBalancingMethod::compute_info(const node_ptr&                      my_node,
                               const DynamicList<chemistryProblem>& problems) {


    sendRecvInfo info;

    info.destinations       = {Pstream::myProcNo()};
    info.sources            = {Pstream::myProcNo()};
    info.number_of_problems = {problems.size()};

    // receiver
    if (my_node->parent->balanced_load.rank != -1) {
        info.sources.push_back(my_node->parent->balanced_load.rank);
        return info;
    } 
    // sender
    else {

        for (const auto& child : my_node->children) {
            info.destinations.push_back(child->balanced_load.rank);
        }
        info.number_of_problems = compute_send_counts(my_node, problems);
    }

    return info;
}

node_ptr simpleBalancingMethod::build_tree(const DynamicList<chemistryLoad>& loads) {


    auto sum_op = [](double sum, const chemistryLoad& rhs) { return sum + rhs.value; };
    double mean_load = std::accumulate(loads.begin(), loads.end(), 0.0, 
                           sum_op) / loads.size();

    double treshold = 1.0 * mean_load;

    auto [big, small] = divide(loads, treshold);
    std::sort(big.begin(), big.end());
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
    const std::vector<double>& times, const DynamicList<chemistryProblem>& problems) {

    std::vector<int> counts; counts.reserve(times.size());
    auto             begin = problems.begin();

    for (const auto& time : times) {

        double sum       = 0.0;
        auto   operation = [&](const chemistryProblem& problem) {
            sum += problem.cpuTime;
            return sum <= time;
            //return sum <= (time + 0.01 * time); // TODO: fix the 0.01
        };
        auto count = count_while(begin, problems.end(), operation);
        begin += count;
        counts.push_back(count);
    }

    return counts;
}

chemistryLoad
simpleBalancingMethod::compute_my_load(const DynamicList<chemistryProblem>& problems) {

    auto   lambda = [](double sum, const chemistryProblem& rhs) { return sum + rhs.cpuTime; };
    double sum    = std::accumulate(problems.begin(), problems.end(), 0.0, lambda);
    return chemistryLoad(Pstream::myProcNo(), sum);
}

std::pair<std::vector<chemistryLoad>, std::vector<chemistryLoad>>
simpleBalancingMethod::divide(const DynamicList<chemistryLoad>& in, double treshold) {

    std::vector<chemistryLoad> big, small;

    big.reserve(in.size());
    small.reserve(in.size());

    for (const auto& v : in) {
        if (v >= treshold) {
            big.push_back(v);
        } else {
            small.push_back(v);
        }
    }

    return {big, small};
}

} // namespace Foam