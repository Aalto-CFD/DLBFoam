#include "simpleLoadBalancing.H"

namespace Foam{


chemistryLoadBalancingMethod::sendRecvInfo simpleLoadBalancing::determine_state(const DynamicList<chemistryLoad>& loads) const{

    //std::vector<int> sources;            // ranks which send to this process
    //    std::vector<int> destinations;       // ranks to which this process sends to
    //    std::vector<int> number_of_problems; // number of problems which this rank sends/receivs

    size_t my_idx = get_my_load_index(loads);

    if (my_idx == 0){
        sendRecvInfo ret;
        auto load_pair = loads[loads.size() - 1];

        ret.sources = { load_pair.rank }; //receive from busiest
        ret.destinations = {};
        ret.number_of_problems = {3};
        return ret;
    }

    if (my_idx == loads.size() - 1){
        sendRecvInfo ret;
        auto load_pair = loads[0];
        ret.sources = { };
        ret.destinations = { load_pair.rank }; //send to least busy
        ret.number_of_problems = {3};
        return ret;
    }

    //others do nothing
    sendRecvInfo ret;
    ret.sources = {};
    ret.destinations = {};
    ret.number_of_problems = {};
    return ret;
    
 


}

} //namespace Foam
