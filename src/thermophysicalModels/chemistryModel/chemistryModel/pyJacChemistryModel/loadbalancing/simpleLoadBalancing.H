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
    ///@brief Returns a sorted vector of loads such that the largest load is at index 0
    ///
    ///@return std::vector<chemistryLoad> a sorted vector of loads in descending order. 
    ///        The size of the vector is N=number of processes.  
    ///
    std::vector<chemistryLoad> get_loads() const;

    //std::vector<chemistryLoad> get_loads() const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
