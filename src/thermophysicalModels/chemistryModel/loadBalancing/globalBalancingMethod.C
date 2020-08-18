#include "globalBalancingMethod.H"

namespace Foam {

void globalBalancingMethod::update_state(const DynamicList<chemistryProblem>& problems) {

    auto my_load   = compute_my_load(problems);
    auto all_loads = all_gather(my_load);

    double global_mean = get_mean(all_loads);

    // Info << global_mean << endl;

    auto operations = get_operations(all_loads, global_mean, my_load);
    auto info       = operations_to_info(operations, problems, my_load);

    set_state(info);
}

chemistryLoadBalancingMethod::sendRecvInfo globalBalancingMethod::operations_to_info(const std::vector<Operation>&        operations,
                                           const DynamicList<chemistryProblem>& problems,
                                           const chemistryLoad&                 my_load){

    sendRecvInfo info;

    bool sender = false;
    bool receiver = false;

    for (const auto& op: operations) {
        if (op.from.rank == my_load.rank){
            sender = true;
            break;
        }
    } 

    for (const auto& op: operations) {
        if (op.to.rank == my_load.rank){
            receiver = true;
            break;
        }
    } 

    runtime_assert(!(sender && receiver), "Only sender or receiver should be possible.");


    if (sender) {
        double sum = 0.0;
        std::vector<double> times;
        for (const auto& op : operations){
            info.destinations.push_back(op.to.rank);
            sum += op.value;
            times.push_back(op.value);
        }
        info.number_of_problems = times_to_problem_counts(times, problems);

    }

    else if (receiver) {

        for (const auto& op : operations){
            info.sources.push_back(op.from.rank);
        }

    }

    info.destinations.push_back(Pstream::myProcNo());
    info.sources.push_back(Pstream::myProcNo());


    return info;
    


}

std::vector<int>
globalBalancingMethod::times_to_problem_counts(const std::vector<scalar>&           times,
                                               const DynamicList<chemistryProblem>& problems) {

    std::vector<int> counts;
    counts.reserve(times.size()+1);
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

    // Add any remaining problems to this rank. This should not be required...
    // TODO: fix
    //counts[0] += (problems.size() - std::accumulate(counts.begin(), counts.end(), 0));
    int left = problems.size() - std::accumulate(counts.begin(),counts.end(), 0);
    counts.push_back(left);
    runtime_assert(std::accumulate(counts.begin(), counts.end(), 0) == problems.size(),
                   "Mismatch in the sliced problem count and original problem count.");
    
    return counts;
}



std::vector<globalBalancingMethod::Operation> globalBalancingMethod::get_operations(
    DynamicList<chemistryLoad>& loads, double global_mean, const chemistryLoad& my_load) const {

    std::vector<Operation> operations;

    std::sort(loads.begin(), loads.end());

    auto sender   = loads.end() - 1;
    auto receiver = loads.begin();

    while (sender != receiver) {

        double    send_value = std::min(sender->value - global_mean, global_mean - receiver->value);
        Operation operation{*sender, *receiver, send_value};
        if (sender->rank == my_load.rank || receiver->rank == my_load.rank) {
            operations.push_back(operation);
        }
        sender->value -= send_value;
        receiver->value += send_value;

        if (std::abs(sender->value - global_mean) < SMALL) {
            sender--;
        }

        else if (std::abs(receiver->value - global_mean) < SMALL) {
            receiver++;
        }
    }
    
    return operations;
}

void globalBalancingMethod::apply_operation(DynamicList<chemistryLoad>&             loads,
                                            const globalBalancingMethod::Operation& operation) {

    loads[operation.from.rank].value -= operation.value;
    loads[operation.to.rank].value += operation.value;
}

globalBalancingMethod::Operation globalBalancingMethod::get_operation(const chemistryLoad& sender,
                                                                      const chemistryLoad& receiver,
                                                                      double               mean) {

    double diff1 = mean - receiver.value;
    double diff2 = sender.value - mean;

    runtime_assert(!(diff1 < 0), "Receiver value larger than mean.");
    runtime_assert(!(diff2 < 0), "Sender value smaller than mean.");

    double send_value = std::min(diff1, diff2);

    return Operation{sender, receiver, send_value};
}

chemistryLoad globalBalancingMethod::get_min(const DynamicList<chemistryLoad>& vec) {

    auto comp = [](const chemistryLoad& lhs, const chemistryLoad& rhs) {
        return lhs.value < rhs.value;
    };
    return *std::min_element(vec.begin(), vec.end(), comp);
}

chemistryLoad globalBalancingMethod::get_max(const DynamicList<chemistryLoad>& vec) {
    auto comp = [](const chemistryLoad& lhs, const chemistryLoad& rhs) {
        return lhs.value < rhs.value;
    };
    return *std::max_element(vec.begin(), vec.end(), comp);
}

} // namespace Foam
