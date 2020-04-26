/*
    Write header
*/

#ifndef simpleLoadBalancing_H
#define simpleLoadBalancing_H



#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"


#include "chemistryLoadBalancingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class simpleLoadBalancing : public chemistryLoadBalancingMethod
{

public:


    simpleLoadBalancing() = default;

    void apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const;

    void print_parameters() const {
        std::cout << "I am balancing load" << std::endl;
    }

    bool active() const{return true;}
   
    
    std::vector<chemistryLoad> get_loads() const;

protected:

    // None

private:

    ///
    ///@brief Gets the chemistry load of this process
    ///
    ///@param a reference cell mapper that tells how many cells are active 
    ///@param Y mass fractions
    ///@return chemistryLoad load information
    ///
    chemistryLoad get_load(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const;
    

    ///
    ///@brief Determines if a process is a sender receiver or do-nothinger
    ///
    ///@param loads a sorted array of loads 
    ///@return WHATTODO state information on what the process should do
    ///
    WHATTODO determine_state(const std::vector<chemistryLoad>& loads) const;



    //TODO: implement
    int get_send_info() const {return 2;}
    int get_recv_info() const {return 1;}



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
