#include "simpleLoadBalancing.H"

namespace Foam{



chemistryLoad simpleLoadBalancing::get_my_load() const{

    throw "NOT IMPLEMENTED";

    chemistryLoad load;
    return load;

}


chemistryLoadBalancingMethod::sendRecvInfo simpleLoadBalancing::determine_state(const DynamicList<chemistryLoad>& loads) const{

    throw "NOT IMPLEMENTED";

    /*
    size_t n = loads.size();

    
    int my_order = 0;
    for (int i = 0; i < n; ++i){
        if (loads[i].rank == Pstream::Pstream::myProcNo()){
            my_order = i;
            break;
        }

    }

    //TODO:: make do nothing possible
    if (my_order < 0.5 * n){
        return WHATTODO::e_SENDER;
    }
    */
    return sendRecvInfo{};
 


}

} //namespace Foam
