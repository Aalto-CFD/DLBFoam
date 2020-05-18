#include "simpleLoadBalancing.H"

namespace Foam{


chemistryLoadBalancingMethod::sendRecvInfo simpleLoadBalancing::determine_state(const DynamicList<chemistryLoad>& loads) const{

    //std::vector<int> sources;            // ranks which send to this process
    //    std::vector<int> destinations;       // ranks to which this process sends to
    //    std::vector<int> number_of_problems; // number of problems which this rank sends/receivs
    sendRecvInfo ret;
    ret.sources = {};
    ret.destinations = {};
    ret.number_of_problems = {};

    return sendRecvInfo{};
 


}

} //namespace Foam
