#include "chemistryLoadBalancingMethod.H"

namespace Foam {


chemistryLoad chemistryLoadBalancingMethod::get_my_load(const DynamicList<chemistryProblem>& problems) const{

    //this sets value = n_active cells
    return chemistryLoad(Pstream::myProcNo(), double(problems.size()), problems.size());

}


size_t chemistryLoadBalancingMethod::get_my_load_index(const DynamicList<chemistryLoad>& loads) const {

    for (size_t i = 0; i < loads.size(); ++i){
        if (loads[i].rank == Pstream::myProcNo()) return i;
    }

    throw "Could not find my rank from loads.";
}



DynamicList<chemistryLoad> chemistryLoadBalancingMethod::get_loads(const DynamicList<chemistryProblem>& problems) const {

    label nprocs = Pstream::nProcs();

    DynamicList<chemistryLoad> loads(nprocs, chemistryLoad());

    loads[Pstream::myProcNo()] = this->get_my_load(problems);



    int tag                    = 1;
    // TODO: Call MPI_Allgather instead, these functions appear to be doing something crazy
    Pstream::gatherList(loads, tag);
    Pstream::scatterList(loads, tag);

    std::sort(loads.begin(), loads.end());

    return loads;
}

DynamicList<DynamicList<chemistryProblem>> chemistryLoadBalancingMethod::get_send_buffer(
    const DynamicList<chemistryProblem>& all_problems) const {

    DynamicList<DynamicList<chemistryProblem>> send_buffer;

    // Only calculate send buffer for sender ranks
    if(m_current_state.sources.size()==0 && m_current_state.destinations.size()!=0)
    {
        int n_problems = all_problems.size();

        int total_send_count = std::accumulate(
            m_current_state.number_of_problems.begin(), m_current_state.number_of_problems.end(), 0);


        // TODO: only in debug mode, and if sender
        if (n_problems < total_send_count) {
            throw "number of send problems more than the current problem count";
        }

        // number of own problems this process will solve
        int n_my_problems = n_problems - total_send_count;

        int start = n_my_problems;
        for (const auto& n : m_current_state.number_of_problems) {

            int                           end = start + n;
            DynamicList<chemistryProblem> slice(n, chemistryProblem());
            std::copy(all_problems.begin() + start, all_problems.begin() + end, slice.begin());
            send_buffer.append(slice);
            start = end;
        }
    }
    return send_buffer;
}

} // namespace Foam
