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


    simpleRefMapping() = default;

    simpleRefMapping(const dictionary& dict, basicSpecieMixture& composition) :
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
    

    int compute_active_cells(PtrList<volScalarField>& Y) const;

    bool shouldMap(chemistryProblem& problem);
    
    bool check_if_refcell(chemistryProblem& problem);

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


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
