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



private:

    ///
    ///@brief Determines the load of this process
    ///
    ///@return chemistryLoad load of this process
    ///
    chemistryLoad get_my_load() const;

    ///
    ///@brief Determines if a process is a sender receiver or do-nothinger
    ///
    ///@param loads a sorted array of loads 
    ///@return WHATTODO state information on what the process should do
    ///
    WHATTODO determine_state(const std::vector<chemistryLoad>& loads) const;






};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
