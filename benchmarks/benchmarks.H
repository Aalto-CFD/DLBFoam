#pragma once

#include <string>

//#include "thermoPhysicsTypes.H"
#include "LoadBalancedChemistryModel.H"
#include "StandardChemistryModel.H"
#include "noChemistrySolver.H"
#include "ode.H"
#include "initial_conditions.H"
#include "benchmark_info.H"


enum class ModelType{standard, balanced};


struct Benchmark{

    Benchmark(BenchmarkInfo info, ModelType model_type, psiReactionThermo& thermo)
    :
    m_info(info)
    {
        
        if (model_type == ModelType::standard){
            m_model = new ode<StandardChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
        }
        else {
            m_model = new ode<LoadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
        }
    }

    ~Benchmark() {
        delete m_model;
    }


    const BenchmarkInfo& get_info() const {return m_info;}


    BasicChemistryModel<psiReactionThermo>* get_model() {return m_model;}

private:
    BenchmarkInfo m_info;
    BasicChemistryModel<psiReactionThermo>* m_model;    

};



struct BenchmarkSolve : public Benchmark{

    BenchmarkSolve(BenchmarkInfo info, ModelType model_type, psiReactionThermo& thermo) : Benchmark(info, model_type, thermo) {}


    void run() {
        this->get_model()->solve(detail::FLOWTIMESTEP);
    }


};





struct BenchmarkSolveSingle{

    static constexpr size_t n_problems = 10;

    BenchmarkSolveSingle(BenchmarkInfo info, psiReactionThermo& thermo, bool heavy=true) 
    : m_info(info), m_problems(create_problems(heavy)), m_solutions(create_solutions())
    {
        m_model = new ode<LoadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
    }


    ~BenchmarkSolveSingle() {delete m_model; }


    const BenchmarkInfo& get_info() const {return m_info;}


    std::vector<ChemistryProblem> create_problems(bool heavy) const {

        if (heavy) {
            return detail::create_n_random_heavy(n_problems);
        }

        return detail::create_n_random_light(n_problems);

    }

    std::vector<ChemistrySolution> create_solutions() const {
        std::vector<ChemistrySolution> solutions;
        for (size_t i=0; i < n_problems; ++i){
            solutions.push_back(detail::create_empty_solution());
        }
        return solutions;
    }

    void run() {

        for (size_t i = 0; i < m_problems.size(); ++i){
            m_model->solveSingle(m_problems[i], m_solutions[i]);
        }

    }

private:

    BenchmarkInfo m_info;
    std::vector<ChemistryProblem> m_problems;
    std::vector<ChemistrySolution> m_solutions;

    ode<LoadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>* m_model;

};
