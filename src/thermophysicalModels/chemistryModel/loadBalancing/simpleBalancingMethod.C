#include "simpleBalancingMethod.H"

namespace Foam {

void simpleBalancingMethod::update_state(const DynamicList<chemistryProblem>& problems) {

    auto my_load   = compute_my_load(problems);
    auto all_loads = all_gather(my_load);

    auto root = build_tree(all_loads);

    loadTree::set_balanced_to_original(root);

    //break early if no balancing is required
    if (root->children.size() == 0) {
    
        sendRecvInfo info;
        info.destinations       = {Pstream::myProcNo()};
        info.sources            = {Pstream::myProcNo()};
        info.number_of_problems = {problems.size()};
        set_state(info);
        return;
    }


    balance_tree(root);

    loadTree::print(root);

    auto my_node = loadTree::find(root, Pstream::myProcNo());
    runtime_assert(my_node != nullptr, "My node not found from the node tree");
    auto info = compute_info(my_node, problems);

    set_state(info);
}

void simpleBalancingMethod::balance_tree(node_ptr& root) {

    auto senders = root->children;

    for (const auto& sender : senders) {

        sender->balanced_load = sender->original_load;

        scalar load_avg(sender->original_load.value);
        for (const auto& child : sender->children) { load_avg += child->original_load.value; }

        //load_avg /= sender->children.size();
        load_avg /= (sender->children.size() + 1);

        for (const auto& child : sender->children) {

            scalar send_value = load_avg - child->original_load.value;

            if ((sender->balanced_load.value - send_value) > load_avg) {

                sender->balanced_load.value -= send_value;
                child->balanced_load.value += send_value;
            }

            else {
                child->balanced_load = child->original_load;
                break;
            }
        }
    }
}

std::vector<int>
simpleBalancingMethod::compute_send_counts(const node_ptr&                      sender,
                                           const DynamicList<chemistryProblem>& problems) {

    std::vector<scalar> send_times;
    size_t              n_children = sender->children.size();
    send_times.reserve(n_children + 1);

    // count the send times towards children
    for (const auto& child : sender->children) {

        scalar time_diff = child->balanced_load.value - child->original_load.value;
        send_times.push_back(time_diff);
    }
    
    // count the remaining send time on the sender
    scalar send_total = std::accumulate(send_times.begin(), send_times.end(), scalar(0));
    scalar remaining  = sender->original_load.value - send_total;
    send_times.insert(send_times.begin(), remaining);

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

    auto   sum_op    = [](scalar sum, const chemistryLoad& rhs) { return sum + rhs.value; };
    scalar mean_load = std::accumulate(loads.begin(), loads.end(), scalar(0), sum_op) / loads.size();

    scalar treshold = (1.0 / HIGH_LOW_LOAD_SEPARATION_TRESHOLD) * mean_load;

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

std::vector<int>
simpleBalancingMethod::times_to_problem_counts(const std::vector<scalar>&           times,
                                               const DynamicList<chemistryProblem>& problems) {

    std::vector<int> counts;
    counts.reserve(times.size());
    auto begin = problems.begin();

    for (const auto& time : times) {

        scalar sum(0);
        auto   operation = [&](const chemistryProblem& problem) {
            sum += problem.cpuTime;
            return sum <= time;
        };
        auto count = count_while(begin, problems.end(), operation);
        begin += count;
        counts.push_back(count);
    }

    //Add any remaining problems to this rank. This should not be required...
    // TODO: fix
    counts[0] += (problems.size() - std::accumulate(counts.begin(), counts.end(), 0));


    runtime_assert(std::accumulate(counts.begin(), counts.end(), 0) == problems.size(),
                   "Mismatch in the sliced problem count and original problem count.");

    return counts;
}

chemistryLoad
simpleBalancingMethod::compute_my_load(const DynamicList<chemistryProblem>& problems) {

    auto   lambda = [](scalar sum, const chemistryProblem& rhs) { return sum + rhs.cpuTime; };
    scalar sum    = std::accumulate(problems.begin(), problems.end(), scalar(0), lambda);
    return chemistryLoad(Pstream::myProcNo(), sum);
}

std::pair<std::vector<chemistryLoad>, std::vector<chemistryLoad>>
simpleBalancingMethod::divide(const DynamicList<chemistryLoad>& in, scalar treshold) {

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