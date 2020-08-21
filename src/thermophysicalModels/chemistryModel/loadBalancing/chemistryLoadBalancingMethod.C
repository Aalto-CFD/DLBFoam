#include "chemistryLoadBalancingMethod.H"

namespace Foam {

/*
chemistryLoad chemistryLoadBalancingMethod::get_my_load(const DynamicList<chemistryProblem>&
problems) const{

    //this sets value = n_active cells
    return chemistryLoad(Pstream::myProcNo(), double(problems.size()), problems.size());

}
*/

scalar chemistryLoadBalancingMethod::get_mean(const DynamicList<chemistryLoad>& loads) {

    auto op = [](scalar sum, const chemistryLoad& load) { return sum + load.value; };
    return std::accumulate(loads.begin(), loads.end(), 0.0, op) / loads.size();
}
chemistryLoad chemistryLoadBalancingMethod::get_min(const DynamicList<chemistryLoad>& vec) {

    auto comp = [](const chemistryLoad& lhs, const chemistryLoad& rhs) {
        return lhs.value < rhs.value;
    };
    return *std::min_element(vec.begin(), vec.end(), comp);
}

chemistryLoad chemistryLoadBalancingMethod::get_max(const DynamicList<chemistryLoad>& vec) {
    auto comp = [](const chemistryLoad& lhs, const chemistryLoad& rhs) {
        return lhs.value < rhs.value;
    };
    return *std::max_element(vec.begin(), vec.end(), comp);
}

label chemistryLoadBalancingMethod::rank_to_load_idx(const DynamicList<chemistryLoad>& loads,
                                                     label                             rank) const {

    /*
    auto iter = std::find_if(loads.begin(), loads.end(),
                            [](const chemistryLoad& l){return l.rank == rank;} );

    return iter - loads.begin();
    */
    for (label i = 0; i < loads.size(); ++i) {
        if (loads[i].rank == rank) return i;
    }

    throw "Could not find my rank from loads.";
}

/*
DynamicList<chemistryLoad> chemistryLoadBalancingMethod::get_loads(const
DynamicList<chemistryProblem>& problems) const {

    auto loads = all_gather(this->get_my_load(problems));
    std::sort(loads.begin(), loads.end());

    return loads;
}
*/

} // namespace Foam
