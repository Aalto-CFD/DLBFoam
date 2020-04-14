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

class simpleRefMapping : public chemistryRefMappingMethod
{

public:


    simpleRefMapping(const dictionary& dict, basicSpecieMixture& composition) :
    chemistryRefMappingMethod(),
    dict_(dict),
    coeffsDict_(dict.subDict("refmapping")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
    mixture_fraction_(coeffsDict_.subDict("mixtureFractionProperties"), composition.species())
    {
        //TODO: move to the mixture fraction constructor and dont call anything here
        mixture_fraction_.update(composition);
    }
    

    bool applyMapping(PtrList<volScalarField>& Y, const label celli);
    
    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli);

    bool active() const {return active_;}
    
protected:

    // None

private:

    // None
    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;


    mixtureFraction mixture_fraction_;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
