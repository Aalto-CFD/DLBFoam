#define PICOBENCH_IMPLEMENT
#include "../third_party/picobench.hpp"


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
    
    //////////////////////

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

    autoPtr<psiReactionThermo> pThermo = psiReactionThermo::New(mesh);
    psiReactionThermo& thermo = pThermo();

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    

    volScalarField rho        
    =    volScalarField(
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh
            ),
            thermo.rho()
        );


    volVectorField U = 
    volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

    volScalarField& p = thermo.p();


    surfaceScalarField phi = 
    surfaceScalarField(
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );


    autoPtr<compressible::turbulenceModel> turbulence = 
            compressible::turbulenceModel::New
            (
                rho,
                U,
                phi,
                thermo
            );


    //for pyJac, ensure that the correct ODE is used

    //BasicChemistryModel<psiReactionThermo>* model1 
    //= new ode<loadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
    

    BasicChemistryModel<psiReactionThermo>* model1 = nullptr;

    //BasicChemistryModel<psiReactionThermo>* model1 
    //= new ode<aaltoChemistryModelBase<psiReactionThermo, gasHThermoPhysics>>(thermo);
    
    BasicChemistryModel<psiReactionThermo>* model2 
    = new ode<StandardChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);



    volScalarField dpdt = 
    volScalarField (
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );

    volScalarField K = 0.5*magSqr(U);


    //multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    //forAll(Y, i)
    //{
    //    fields.add(Y[i]);
    //}
    //fields.add(thermo.he());


    volScalarField Qdot
    = 
    volScalarField
    (
        IOobject
        (
            "Qdot",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
    );


    ///////////////////////////////////////////////////////////////


    static argList create_args() {

        /*
        if (!Pstream::parRun()){
                
            int argc = 1;
            char* argv[] = {"benchmarkChemistry"}; // exactly as it is defined in main
            char** ptr = argv;

            return argList(argc, ptr);
        }
        */

        

            int argc = 2;
            char* argv[] = {"benchmarkChemistry", "-parallel"}; // exactly as it is defined in main
            char** ptr = argv;

            return argList(argc, ptr);
        

        /*

            int argc = 1;
            char* argv[] = {"benchmarkChemistry"}; // exactly as it is defined in main
            char** ptr = argv;

            return argList(argc, ptr);

        */

    }


};


static ChemState global_chemistry_state;

struct InitialCondition1{

    static void set() {

        Random  randomizer1(541); 
        Random  randomizer2(44); 
        Random  randomizer3(10);

        forAll(global_chemistry_state.rho, i) {
            global_chemistry_state.rho[i] = 1.2 + randomizer1.sample01<scalar>() * 3;
            global_chemistry_state.thermo.T()[i] = 500 + randomizer2.sample01<scalar>() * 100;
        }

        forAll(global_chemistry_state.rho, celli) {

            for (label i = 0; i < global_chemistry_state.Y.size(); i++) { 
                global_chemistry_state.Y[i][celli] = randomizer3.sample01<scalar>();
                //c_[i] = Y_[i][celli]; 
            }

        }


    }


};




static void standard_chemistry_solve(picobench::state& s)
{

    InitialCondition1::set();


    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){

        global_chemistry_state.model2->solve(1E-5);
    }

}

PICOBENCH(standard_chemistry_solve);

static void pyjac_chemistry_solve(picobench::state& s)
{

    InitialCondition1::set();

    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){
        global_chemistry_state.model1->solve(1E-5);
    }
}

PICOBENCH(pyjac_chemistry_solve);





int main(int argc, char* argv[]){


    //ode<aaltoChemistryModelBase<psiReactionThermo, gasHThermoPhysics>> asd(global_chemistry_state.thermo);

    picobench::runner runner;
    // Disregard command-line for simplicity

    // Two sets of iterations
    runner.set_default_state_iterations({10, 50});

    // One sample per benchmark because the huge numbers are expected to compensate
    // for external factors
    runner.set_default_samples(1);

    // Run the benchmarks with some seed which guarantees the same order every time
//    auto report = runner.run_benchmarks(123);

    return runner.run();

    //runner.run_benchmarks(123);

    // Output to some file
//    report.to_csv(std::ofstream("my.csv"));

    return 0;


//    reaction->correct();


    





}





