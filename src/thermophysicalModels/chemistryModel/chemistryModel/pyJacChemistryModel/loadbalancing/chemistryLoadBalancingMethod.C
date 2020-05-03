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
std::vector<chemistryLoad> chemistryLoadBalancingMethod::get_loads() const{

    
    int nprocs = Pstream::nProcs();

    List<int>      ranks(nprocs);
    List<int>      active_cells(nprocs);
    List<double>   vals(nprocs);

    auto my_load = this->get_my_load();

    ranks[Pstream::myProcNo()] = my_load.rank;
    active_cells[Pstream::myProcNo()] = my_load.number_of_active_cells;
    vals[Pstream::myProcNo()] = my_load.value;

    //This is now three calls to MPI_Allgather, if a List<chemistryLoads> was communicated it could be done with a single call
    //That would require serialization of the chemistryLoad class
    //TODO: serialize chemistryLoad and pass the objects directly
    int tag = 1;
    Pstream::gatherList(ranks, tag);
    Pstream::gatherList(active_cells, tag);
    Pstream::gatherList(vals, tag);

    std::vector<chemistryLoad> loads;
    loads.reserve(nprocs);

    for (int i = 0; i < nprocs; ++i){

        chemistryLoad load;
        load.rank = ranks[i];
        load.number_of_active_cells = active_cells[i];
        load.value = vals[i];

        loads.emplace_back(load);
    }    

 
    std::sort(loads.begin(), loads.end());

    return loads;
    
}



void chemistryLoadBalancingMethod::send_recv(const DynamicList<chemistryProblem>& problems, int source, int destination){



    /*
    PstreamBuffers pBufs(Pstream::commsTypes::blocking);
    UOPstream toNbr(destination, pBufs);

    toNbr  << problems;

    pBufs.finishedSends();

    UIPstream fromNeighb(source, pBufs);
    */

    //fromNeighb >> nbrdAlpha >> nbrpFaceCells >> nbrp_rghPross >> nbrCellC ;


    /*
    Pstream::exchange(
        problems,

    );
    */



}


} //namespace Foam



