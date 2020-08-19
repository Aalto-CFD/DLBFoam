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

    virtual ~bulutLoadBalancing() {}

    bool active() const {return true;}

private:

    ///
    ///@brief Given a list of all problems, updates the sendRecvInfo member
    ///
    ///@param problems a list of all problems on a rank
    ///
    virtual void update_state(const DynamicList<chemistryProblem>& problems);

    //sendRecvInfo determine_state(DynamicList<chemistryLoad>& loads,const DynamicList<chemistryProblem>& all_problems) const override;



private:


    //TODO: make user input
    static constexpr scalar m_tolerance_to_balance = 0.1;

        
    scalar cpu_time_to_send(chemistryLoad& sender_load, chemistryLoad& receiver_load, scalar& meanCPUT) const;
    void update_load(chemistryLoad& sender_load, chemistryLoad& receiver_load, scalar& meanCPUT) const;

    int cpu_time_to_ncells(const scalar& load_to_send, const DynamicList<chemistryProblem>& all_problems, int& index) const;


    bool large_sender(const chemistryLoad& load, const scalar& meanCPUT) const;
    bool small_receiver(const chemistryLoad& load, const scalar& meanCPUT) const;
    scalar get_mean_load(DynamicList<chemistryLoad>& loads) const;
    sendRecvInfo initialize_to_nothinger(const DynamicList<chemistryProblem>& problems) const;
    void update_self_send(sendRecvInfo& ret) const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam