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
    Pstream::gatherList(loads, tag);
    Pstream::scatterList(loads, tag);

    std::sort(loads.begin(), loads.end());

    return loads;
}




} // namespace Foam
