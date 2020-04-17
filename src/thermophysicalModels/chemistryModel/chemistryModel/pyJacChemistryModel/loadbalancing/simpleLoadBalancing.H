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
#include "load.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class simpleLoadBalancing : public chemistryLoadBalancingMethod
{

public:


    simpleLoadBalancing() = default;


    void print_parameters() const {
        std::cout << "I am balancing load" << std::endl;
    }

    bool active() const{return true;}
   

protected:

    // None

private:

    
    Load get_load() const;



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
