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


    simpleLoadBalancing(const dictionary& dict, basicSpecieMixture& composition) :
    chemistryLoadBalancingMethod(),
    dict_(dict),
    coeffsDict_(dict.subDict("loadbalancing")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
    mpiInfoTableSize_(coeffsDict_.lookupOrDefault<scalar>("mpiInfoTableSize",5000)),
    mpiBufferLimit_(coeffsDict_.lookupOrDefault<scalar>("mpiBufferLimit",50000)),
    chemCPUTimeLimit_(coeffsDict_.lookupOrDefault<scalar>("chemCPUTimeLimit",5))
    
    {
        //TODO: move to the mixture fraction constructor and dont call anything here
        //mixture_fraction_.update(composition);
    }
    

    bool applyLoadBalancing();
    
    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli);

    labelListList compute_stats();
    //labelListList compute_stats();

    bool active() const {return active_;}

    void print_parameters() 
    {
        Info<<"Load balancing is initialized with:\n "<<
        "\t mpiInfoTableSize: "<<mpiInfoTableSize_<<
        "\n\t mpiBufferLimit: "<<mpiBufferLimit_<<
        "\n\t chemCPUTimeLimit: "<<chemCPUTimeLimit_<<" seconds."<<endl;
    }
    
protected:

    // None

private:

    // None
    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is reference mapping active?
    Switch active_;

    //- Tolerance for reference mapping
    scalar tolerance_;

    label mpiInfoTableSize_;
    label mpiBufferLimit_;
    scalar chemCPUTimeLimit_;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
