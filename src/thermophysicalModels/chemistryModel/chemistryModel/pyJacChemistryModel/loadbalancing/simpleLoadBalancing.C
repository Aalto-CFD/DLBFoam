#include "simpleLoadBalancing.H"

namespace Foam{

void simpleLoadBalancing::apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const{

    auto loads = get_loads(mapper, Y);

    WHATTODO state = determine_state(loads);

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
    


}

chemistryLoad simpleLoadBalancing::get_load(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const{


    int active_cells = mapper->compute_active_cells(Y);
    /*
        do other stuff
    */
    double load_val = Y[0].size() / active_cells;
    int my_rank = 3;

    return chemistryLoad(my_rank, load_val);

}

//TODO make work
std::vector<chemistryLoad> simpleLoadBalancing::get_loads(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const{


    //auto my_load = get_load(mapper, Y);



    //TODO replace with some
    //int nprocs = get_world_size();
    int nprocs = Pstream::nProcs();

    std::vector<chemistryLoad> ret;
    ret.reserve(nprocs);

    // CALL MPI_ALLGATHER-like function. The Pstream::gather seems to be templated to do all kinds of things and may be useful.


    for (int i = 0; i < nprocs; ++i){
        chemistryLoad load;
        load.rank = i;
        load.number_of_active_cells = i * 42;
        load.value = 321.0 * i;
        ret.push_back(load);
    }

    return ret;
}



WHATTODO simpleLoadBalancing::determine_state(const std::vector<chemistryLoad>& loads) const{

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
