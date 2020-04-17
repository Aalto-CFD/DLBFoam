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
    chemCPUTimeLimit_(coeffsDict_.lookupOrDefault<scalar>("chemCPUTimeLimit",5)),
    nPs_(Pstream::nProcs()),
    nVar_(5),
    maxCPUt(0),
    minCPUt(0),
    meanCPUt(0),  
    overheadT(nPs_, 0.0),
    nc_rm(nPs_, 0.0),
    t_cpu_list(nPs_,-SMALL),
    ncells_list(nPs_,0.0),
    nActiveCells_list(nPs_,0.0),
    receiverInfo_all(nPs_, List<label>(mpiInfoTableSize_,0.0) ),
    senderInfo_all(nPs_, List<label>(mpiInfoTableSize_,0.0) ),
    receiverInfo(mpiInfoTableSize_,0.0),
    senderInfo(mpiInfoTableSize_,0.0)
    {
        //TODO: move to the mixture fraction constructor and dont call anything here
        //mixture_fraction_.update(composition);
    }
    

    bool applyLoadBalancing();

    labelList getLoadBalStats(scalarField& t_cpu_list, labelField& ncells_list, labelField& nActiveCells_list);

    
    void compute_global_stats();

    void get_overhead();
    void do_balancing_calc();
    void get_balancing_info();

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

    label nPs_;
    label nVar_;

    scalar maxCPUt;
    scalar minCPUt;
    scalar meanCPUt; 
    scalarList overheadT;
    scalarList nc_rm;

    scalarField t_cpu_list;
    labelField ncells_list;
    labelField nActiveCells_list;
    List<List<label>> receiverInfo_all;
    List<List<label>> senderInfo_all;
    labelList senderInfo;
    labelList receiverInfo;

    

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
