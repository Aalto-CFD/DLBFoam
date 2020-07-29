#pragma once

#include <string>

#include "thermoPhysicsTypes.H"
#include "loadBalancedChemistryModel.H"
#include "StandardChemistryModel.H"
#include "noChemistrySolver.H"
#include "ode.H"
#include "initial_conditions.H"

enum class ModelType{standard, balanced};


struct Benchmark{

    template<class InitialCondition>
    Benchmark(std::string name, ModelType model_type, psiReactionThermo& thermo, InitialCondition ic)
    :
    m_name(name)
    {
        
        InitialCondition::set(thermo);

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

    std::string name() const {
        return m_name; 
    } 


    BasicChemistryModel<psiReactionThermo>* get_model() {return m_model;}

private:
    std::string m_name;
    BasicChemistryModel<psiReactionThermo>* m_model;    

};


struct BenchmarkSolveSingle{

    static constexpr size_t n_problems = 200;

    BenchmarkSolveSingle(std::string name, psiReactionThermo& thermo, bool heavy=true) 
    : m_name(name), m_problems(create_problems(heavy)), m_solutions(create_solutions())
    {
        m_model = new ode<loadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
    }


    ~BenchmarkSolveSingle() {delete m_model; }

    std::string name() const {
        return m_name; 
    } 


    std::vector<chemistryProblem> create_problems(bool heavy) const {

        if (heavy) {
            return detail::create_n_random_heavy(n_problems);
        }

        return detail::create_n_random_light(n_problems);

    }

    std::vector<chemistrySolution> create_solutions() const {
        std::vector<chemistrySolution> solutions;
        for (size_t i=0; i < n_problems; ++i){
            solutions.push_back(detail::create_empty_solution());
        }
        return solutions;
    }

    void run() {

        for (size_t i = 0; i < m_problems.size(); ++i){
            m_model->solve_single(m_problems[i], m_solutions[i]);
        }

    }

private:

    std::string m_name;
    std::vector<chemistryProblem> m_problems;
    std::vector<chemistrySolution> m_solutions;

    ode<loadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>* m_model;

};


struct BenchmarkSolve : public Benchmark{

    template<class InitialCondition>
    BenchmarkSolve(std::string name, ModelType model_type, psiReactionThermo& thermo, InitialCondition ic) : Benchmark(name, model_type, thermo, ic) {}


    void run() {
        this->get_model()->solve(1E-6);
    }


};