#include <vector>

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


///
#include "thermoPhysicsTypes.H"
#include "Random.H"
#include "loadBalancedChemistryModel.H"
#include "StandardChemistryModel.H"
#include "noChemistrySolver.H"
#include "ode.H"

enum class ModelType{standard, balanced};

struct Result{

    Result(const std::vector<double>& times)
    : average(calc_avg(times)), standard_dev(calc_std(times))
    {}

    double calc_avg(const std::vector<double>& times) const{
        return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    }

    double calc_std(const std::vector<double>& times) const {
        return 0.0;
    }

    double average, standard_dev;

};

struct Benchmark{

    Benchmark(ModelType model_type, psiReactionThermo& thermo) {
        if (model_type == ModelType::standard){
            m_model = new ode<StandardChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
        }
        else {
            m_model = new ode<loadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
        }
    }

    ~Benchmark() {
        delete m_model;
    }


    BasicChemistryModel<psiReactionThermo>* get_model() {return m_model;}

private:
    BasicChemistryModel<psiReactionThermo>* m_model;    

};

struct Runner{


    template<class BenchmarkType>
    static Result run(BenchmarkType benchmark, size_t n_times) {
        
        std::vector<double> times;
        clockTime time;
        for (size_t i = 0; i < n_times; ++i){
            time.timeIncrement();
            benchmark.run();
            times.push_back(time.timeIncrement());
        }
        return Result(times);
    }

};




struct Benchmark1 : public Benchmark{

    Benchmark1(ModelType model_type, psiReactionThermo& thermo) : Benchmark(model_type, thermo) {}



    void run() {
        this->get_model()->solve(1E-3);
    }


};





int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"


    auto result1 = Runner::run(Benchmark1(ModelType::standard, thermo), 10);

    Info << "result 1 mean: " << result1.average << endl;


    //Runner r1(B1<ModelType::standard>(m1), 10);


    /*
    auto times = Benchmark(B1(m1), 10).times;

    for (auto t : times) {
        Info << t << endl;
    }
    */

    return 0;


}

