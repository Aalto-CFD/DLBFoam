#include "chemistryLoadBalancingMethod.H"

namespace Foam{


void chemistryLoadBalancingMethod::apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const{

    auto loads = this->get_loads();

    
    WHATTODO state = this->determine_state(loads);

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




//TODO: consider breaking to smaller functions
DynamicList<chemistryLoad> chemistryLoadBalancingMethod::get_loads() const{


    label nprocs = Pstream::nProcs();

    DynamicList<chemistryLoad> loads(nprocs, chemistryLoad());
    
    loads[Pstream::myProcNo()] = this->get_my_load();
    int tag = 1;
    Pstream::gatherList(loads, tag);

    std::sort(loads.begin(), loads.end());

    return loads;

}





} //namespace Foam



