#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "Random.H"



//#include "aaltoChemistryModelBase.H"
#include "loadBalancedChemistryModel.H"
#include "StandardChemistryModel.H"
#include "noChemistrySolver.H"
#include "ode.H"
#include "seulex.H"
#include "laminar.H"


using namespace Foam;


enum class parallel{on, off};

enum class chemistryModelName{standard, loadBalanced};



template<parallel flag, chemistryModelName modelName>
struct ChemState{

    
    ChemState() : 
    //args(create_args()), 
    runTime(Foam::Time::controlDictName, args), 
    mesh(   Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
            ),
    pimple(mesh)
    {}


    ~ChemState() {delete model; }


    

    //////////////////////
    //#include "addCheckCaseOptions.H"
    argList args = create_args();
    Time runTime;
    fvMesh mesh;
    pimpleControl pimple;

    //readTimeControls.H
    bool adjustTimeStep = false;
    scalar maxCo = 1.0;
    scalar maxDeltaT = great;
    scalar cumulativeContErr = 0;


    ////////////////// Specific requirements //////////////////////


    volScalarField T        
    =    volScalarField(
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ, //
                IOobject::NO_WRITE //
            ),
            //thermo.rho()
            mesh,
            dimensionedScalar("some", dimTemperature , 1.0)
        );


    volScalarField rho        
    =    volScalarField(
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ, //
                IOobject::NO_WRITE //
            ),
            //thermo.rho()
            mesh,
            dimensionedScalar("some", dimDensity , 1.0)
        );




    volVectorField U = 
    volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("some", dimVelocity, vector(0,0,0))
        );

    volScalarField p        
    =    volScalarField(
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ, //
                IOobject::NO_WRITE //
            ),
            //thermo.rho()
            mesh,
            dimensionedScalar("some", dimPressure , 1.0)
        );


    surfaceScalarField phi = 
    surfaceScalarField(
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );


    autoPtr<psiReactionThermo> pThermo = psiReactionThermo::New(mesh);
    psiReactionThermo& thermo = pThermo();

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();


    autoPtr<compressible::turbulenceModel> turbulence = 
            compressible::turbulenceModel::New
            (
                rho,
                U,
                phi,
                thermo
            );


    //for pyJac, ensure that the correct ODE is used

    
    
    BasicChemistryModel<psiReactionThermo>* model = create_chem_model(thermo); 




    ///////////////////////////////////////////////////////////////


    static argList create_args() {
        /*int argc;
        char** argv;

        #include "setRootCaseLists.H" 
        */
        if (flag == parallel::on){
            int argc = 2;
            char* argv[] = {"benchmarkChemistry", "-parallel"}; // exactly as it is defined in main
            char** ptr = argv;

            return argList(argc, ptr);
        }
        int argc = 1;
        char* argv[] = {"benchmarkChemistry"}; // exactly as it is defined in main
        char** ptr = argv;

        return argList(argc, ptr);
    }

    static BasicChemistryModel<psiReactionThermo>* create_chem_model(psiReactionThermo& thermo) {

        if (modelName == chemistryModelName::standard){
            return new ode<StandardChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
        }

        return new ode<loadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);


    }


};
