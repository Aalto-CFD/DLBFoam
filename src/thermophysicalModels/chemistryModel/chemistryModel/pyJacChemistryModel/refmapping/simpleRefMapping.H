/*
    Write header
*/

#ifndef simpleRefMapping_H
#define simpleRefMapping_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"


#include "chemistryRefMappingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
class pyJacChemistryModel;


template<class CompType, class ThermoType>
class simpleRefMapping : public chemistryRefMappingMethod<CompType, ThermoType>{

public:


    simpleRefMapping(const dictionary& dict) : chemistryRefMappingMethod<CompType, ThermoType>(dict) {
        
    }
    
    void derived_hello() {
        std::cout << "Hello from a derived simpleRefMapping::function" << std::endl;    
    }


    int some_function() {return 32;}


    //Note this is defined here empty.
    virtual void applyMapping(pyJacChemistryModel<CompType, ThermoType>& chemistry) const {

    };

protected:

    // None

private:

    // None


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
