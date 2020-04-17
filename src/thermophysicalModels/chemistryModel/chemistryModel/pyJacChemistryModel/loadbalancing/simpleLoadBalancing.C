#include "simpleLoadBalancing.H"

namespace Foam{

void simpleLoadBalancing::apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const{

    auto loads = get_loads(mapper, Y);


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

} //namespace Foam