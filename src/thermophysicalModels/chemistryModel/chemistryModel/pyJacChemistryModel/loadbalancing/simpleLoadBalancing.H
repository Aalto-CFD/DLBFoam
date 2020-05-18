#pragma once


#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

#include "chemistryProblem.H"
#include "chemistryLoadBalancingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class simpleLoadBalancing : public chemistryLoadBalancingMethod
{

public:


    simpleLoadBalancing() = default;



private:

    ///
    ///@brief Determines if a process is a sender receiver or do-nothinger
    ///
    ///@param loads a sorted array of loads 
    ///@return WHATTODO state information on what the process should do
    ///
    sendRecvInfo determine_state(const DynamicList<chemistryLoad>& loads) const override;






};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam

