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

class simpleRefMapping : public chemistryRefMappingMethod{

public:


    simpleRefMapping(const dictionary& dict, const wordList& species) : chemistryRefMappingMethod(dict, species) {
        
    }
    

    //Note this is defined here empty.
    bool applyMapping();

        ///
    ///@brief Dont know.
    ///
    ///@param composition 
    ///
    //this is a bit useless and should be moved to constructor
    void init_mixture_fraction(basicSpecieMixture& composition);
    

    ///
    ///@brief Checks if a reference cell or not
    ///
    ///@param Y 
    ///@param celli 
    ///@return true 
    ///@return false 
    ///
    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli);

    

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
