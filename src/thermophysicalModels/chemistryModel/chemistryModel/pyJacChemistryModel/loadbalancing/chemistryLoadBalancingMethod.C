#include "chemistryLoadBalancingMethod.H"

namespace Foam {

void chemistryLoadBalancingMethod::apply_balancing(const chemistryRefMappingMethod* mapper,
                                                   PtrList<volScalarField>&         Y) const {

    auto loads = this->get_loads();

    sendRecvInfo state = this->determine_state(loads);

    // WHATTODO state = this->determine_state(loads);

    /*
    //move to a function
    switch(state) {
        case WHATTODO::e_SENDER:
            auto asd = get_send_info();
            //send_stuff(asd);
            break;
        case WHATTODO::e_RECEIVER:
            auto asd = get_recv_info();
            //receive_stuff(asd);
            break;

        default: //This is the DONOTHING
            break;
    }
    */
}

DynamicList<chemistryLoad> chemistryLoadBalancingMethod::get_loads() const {

    label nprocs = Pstream::nProcs();

    DynamicList<chemistryLoad> loads(nprocs, chemistryLoad());

    loads[Pstream::myProcNo()] = this->get_my_load();
    int tag                    = 1;
    // TODO: Call MPI_Allgather instead, these functions appear to be doing something crazy
    Pstream::gatherList(loads, tag);
    Pstream::scatterList(loads, tag);

    std::sort(loads.begin(), loads.end());

    return loads;
}

DynamicList<DynamicList<chemistryProblem>> chemistryLoadBalancingMethod::get_send_buffer(
    const DynamicList<chemistryProblem>& all_problems) const {

    int n_problems = all_problems.size();

    int total_send_count = std::accumulate(
        m_current_state.number_of_problems.begin(), m_current_state.number_of_problems.end(), 0);
    // TODO: only in debug mode
    if (n_problems < total_send_count) {
        throw "number of send problems more than the current problem count";
    }

    // number of own problems this process will solve
    int n_my_problems = n_problems - total_send_count;

    DynamicList<DynamicList<chemistryProblem>> send_buffer;

    int start = n_my_problems;
    for (const auto& n : m_current_state.number_of_problems) {

        int end = start + n;
        DynamicList<chemistryProblem> slice(n, chemistryProblem());
        std::copy(all_problems.begin() + start, all_problems.begin() + end, slice.begin());
        send_buffer.append(slice);
        start = end;
    }
    return send_buffer;
}


chemistryLoadBalancingMethod::buffer_t<chemistryProblem>
chemistryLoadBalancingMethod::get_problems(const DynamicList<chemistryProblem>& problems) const {

    auto sources      = m_current_state.sources;
    auto destinations = m_current_state.destinations;
    auto counts       = m_current_state.number_of_problems;

    // TODO: Check only in debug mode
    if (sources.size() != 0 && destinations.size() != 0) {
        throw "A process can be either a source, destination or a do nothinger";
    }


    auto send_buffer = get_send_buffer(problems);

    //The unittests currently test such that everyone sends to master which requires a blocking commstype
    //in reality we wont block at this phase
    //TODO: change to nonBlocking
    return send_recv<chemistryProblem, Pstream::commsTypes::blocking>(send_buffer, sources, destinations);


}

DynamicList<chemistrySolution>
chemistryLoadBalancingMethod::get_solutions(const chemistryLoadBalancingMethod::buffer_t<chemistrySolution>& solutions) const{

    auto sources      = m_current_state.sources;
    auto destinations = m_current_state.destinations;
    auto counts       = m_current_state.number_of_problems;

    // TODO: Check only in debug mode
    if (sources.size() != 0 && destinations.size() != 0) {
        throw "A process can be either a source, destination or a do nothinger";
    }


    //this probably has to be always blocking
    auto recv_buffer = send_recv<chemistrySolution, Pstream::commsTypes::blocking>(solutions, destinations, sources);


    DynamicList<chemistrySolution> ret;
    return ret;


}


} // namespace Foam
