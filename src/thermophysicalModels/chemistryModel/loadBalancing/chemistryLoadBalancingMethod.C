#include "chemistryLoadBalancingMethod.H"

namespace Foam {

/*
chemistryLoad chemistryLoadBalancingMethod::get_my_load(const DynamicList<chemistryProblem>& problems) const{

    //this sets value = n_active cells
    return chemistryLoad(Pstream::myProcNo(), double(problems.size()), problems.size());

}
*/

size_t chemistryLoadBalancingMethod::rank_to_load_idx(const DynamicList<chemistryLoad>& loads, int rank) const {

    /*
    auto iter = std::find_if(loads.begin(), loads.end(), 
                            [](const chemistryLoad& l){return l.rank == rank;} );

    return iter - loads.begin();
    */
    for (size_t i = 0; i < size_t(loads.size()); ++i){
        if (loads[i].rank == rank) return i;
    }

    throw "Could not find my rank from loads.";
}


/*
DynamicList<chemistryLoad> chemistryLoadBalancingMethod::get_loads(const DynamicList<chemistryProblem>& problems) const {
    
    auto loads = all_gather(this->get_my_load(problems));
    std::sort(loads.begin(), loads.end());

    return loads;
}
*/



} // namespace Foam
