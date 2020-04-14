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


    simpleRefMapping(const dictionary& dict, const wordList& species) : chemistryRefMappingMethod(dict, species) 
    {
        
    }
    

    bool applyMapping(PtrList<volScalarField>& Y, const label celli);
    void init_mixture_fraction(basicSpecieMixture& composition);
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
