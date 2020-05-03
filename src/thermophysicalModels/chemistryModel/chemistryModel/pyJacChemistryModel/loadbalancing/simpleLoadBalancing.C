#include "simpleLoadBalancing.H"

namespace Foam{



chemistryLoad simpleLoadBalancing::get_my_load() const{

    //This is currently arbitrary but defined such that testing is possible

    chemistryLoad load;

    
    load.rank = Pstream::myProcNo();
    load.number_of_active_cells = load.rank * 2;
    load.value = double(3.0 * load.rank);

    
    

    return load;

}


WHATTODO simpleLoadBalancing::determine_state(const DynamicList<chemistryLoad>& loads) const{

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

    return WHATTODO::e_RECEIVER;
 


}

} //namespace Foam
