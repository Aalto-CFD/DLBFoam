/*
    Write header
*/

#ifndef simpleRefMappingMethod_H
#define simpleRefMappingMethod_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"


#include "chemistryRefMappingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class simpleRefMappingMethod : public chemistryRefMappingMethod
{

public:


    simpleRefMappingMethod() = default;

    simpleRefMappingMethod(const dictionary& dict, const basicSpecieMixture& composition) :
    chemistryRefMappingMethod(),
    dict_(dict),
    coeffsDict_(dict.subDict("refmapping")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
    tolerance_(coeffsDict_.lookupOrDefault<scalar>("tolerance",1e-4)),
    mixture_fraction_(coeffsDict_.subDict("mixtureFractionProperties"), composition.species())
    {
        //TODO: move to the mixture fraction constructor and dont call anything here
        mixture_fraction_.update(composition);
    }
    


    bool shouldMap(const scalarField& mass_fraction) const;
    

    bool active() const {return active_;}
    
protected:

    // None

private:

    //TODO: for runtime switching, these need to be references
    const dictionary dict_;
    const dictionary coeffsDict_;

    //- Is reference mapping active?
    Switch active_;

    //- Tolerance for reference mapping
    scalar tolerance_;


    mixtureFraction mixture_fraction_;

    bool check_if_refcell(const scalarField& mass_fraction) const;



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
