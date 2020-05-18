#pragma once

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"


#include "chemistryLoadBalancingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class bulutLoadBalancing : public chemistryLoadBalancingMethod
{

public:


    bulutLoadBalancing() = default;



private:


    //TODO: make user input
    static constexpr scalar m_tolerance_to_balance = 0.02;


    sendRecvInfo determine_state(const DynamicList<chemistryLoad>& loads) const override;

    
    void update_load(chemistryLoad& sender_load, chemistryLoad& receiver_load) const;
    
    int ncells_to_send(chemistryLoad& sender_load, chemistryLoad& receiver_load, int nCellsSentTot) const;

    bool large_sender(const chemistryLoad& load) const;
    bool small_receiver(const chemistryLoad& load) const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


