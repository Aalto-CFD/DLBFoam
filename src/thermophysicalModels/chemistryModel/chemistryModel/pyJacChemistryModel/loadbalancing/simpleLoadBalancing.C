#include "simpleLoadBalancing.H"

namespace Foam{

void simpleLoadBalancing::apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const{

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


std::vector<chemistryLoad> simpleLoadBalancing::get_loads() const{

}

} //namespace Foam