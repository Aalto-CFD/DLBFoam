/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "loadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {

template <class ReactionThermo, class ThermoType>
loadBalancedChemistryModel<ReactionThermo, ThermoType>::loadBalancedChemistryModel(
    ReactionThermo& thermo)
    : StandardChemistryModel<ReactionThermo, ThermoType>(thermo)
    , cpu_times_(this->mesh().cells().size(), 0.0)
    , load_balancer_(create_balancer()) {

    Info << "Running with a load balanced" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
loadBalancedChemistryModel<ReactionThermo, ThermoType>::~loadBalancedChemistryModel(){
    delete load_balancer_; // TODO: use a smart pointer
}



template <class ReactionThermo, class ThermoType>
chemistryLoadBalancingMethod* loadBalancedChemistryModel<ReactionThermo, ThermoType>::create_balancer(){

    return new globalBalancingMethod();
    //return new bulutLoadBalancing();
    //return new simpleBalancingMethod();

}


template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
scalar loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve(const DeltaTType& deltaT) {


    using problem_buffer_t  = chemistryLoadBalancingMethod::buffer_t<chemistryProblem>;
    using solution_buffer_t = chemistryLoadBalancingMethod::buffer_t<chemistrySolution>;


    BasicChemistryModel<ReactionThermo>::correct();

    if (!this->chemistry_) { return great; }


    
    DynamicList<chemistryProblem> all_problems = get_problems(this->Y_, deltaT);

    //I dont think this should be done as this operation is very expensive for large problem counts
    // and ruins up the cellid order
    /* 
    std::sort(all_problems.begin(), all_problems.end(), 
    [](const chemistryProblem& lhs, const chemistryProblem& rhs){
        return lhs.cpuTime > rhs.cpuTime;}
    );
    }
    */ 


    load_balancer_->update_state(all_problems);
    load_balancer_->print_state();
    problem_buffer_t balanced_problems = load_balancer_->balance(all_problems);

    //Pstream::waitRequests(); 

    solution_buffer_t balanced_solutions = solve_buffer(balanced_problems);
    solution_buffer_t my_solutions = load_balancer_->unbalance(balanced_solutions);

    
    //Pstream::waitRequests(); 

    return update_reaction_rates(my_solutions);




}


template <class ReactionThermo, class ThermoType>
void loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve_single(
    chemistryProblem& prob, chemistrySolution& soln) const {

    // TODO: Make this work so that chemistryProblem is a const &




    scalar           timeLeft = prob.deltaT;
    const scalarList c0       = prob.c;



    //Timer begin
    clockTime time;
    time.timeIncrement();

    // Calculate the chemical source terms
    while (timeLeft > small) {
        scalar dt = timeLeft;
        this->solve(prob.c, prob.Ti, prob.pi, dt, prob.deltaTChem);
        timeLeft -= dt;
    }




    soln.c_increment = (prob.c - c0) / prob.deltaT;
    soln.deltaTChem = prob.deltaTChem;
    //Timer end
    soln.cpuTime = time.timeIncrement();
    soln.cellid = prob.cellid;


}



template <class ReactionThermo, class ThermoType>
scalar loadBalancedChemistryModel<ReactionThermo, ThermoType>::update_reaction_rates(const  chemistryLoadBalancingMethod::buffer_t<chemistrySolution>& solutions){


    scalar deltaTMin = great;


    for (const auto& array : solutions){
    for (const auto& solution : array){

        
        for (label j = 0; j < this->nSpecie_; j++) { 
            this->RR_[j][solution.cellid] = solution.c_increment[j] * this->specieThermo_[j].W();                 
        }

        deltaTMin = min(solution.deltaTChem, deltaTMin);
        deltaTMin = min(deltaTMin, this->deltaTChemMax_);
        this->deltaTChem_[solution.cellid] = deltaTMin;
        this->cpu_times_[solution.cellid]  = solution.cpuTime;
    }

    }

    return deltaTMin;


}

template <class ReactionThermo, class ThermoType>
scalar loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve(const scalarField& deltaT) {
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class ThermoType>
scalar loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}








template <class ReactionThermo, class ThermoType>
chemistryLoadBalancingMethod::buffer_t<chemistrySolution>
loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve_buffer(
    chemistryLoadBalancingMethod::buffer_t<chemistryProblem>& problems) const {

    //allocate the solutions buffer
    chemistryLoadBalancingMethod::buffer_t<chemistrySolution> solutions;
    for (int i = 0; i < problems.size(); ++i){
        DynamicList<chemistrySolution> sublist(problems[i].size(), chemistrySolution(this->nSpecie_));
        solutions.append(sublist);
    }


    for (int i = 0; i < solutions.size(); ++i){
        for (int j = 0; j < solutions[i].size(); ++j ){
            
            solve_single(problems[i][j], solutions[i][j]);
        }
    }

    return solutions;
}




template <class ReactionThermo, class ThermoType>
template<class DeltaTType>
DynamicList<chemistryProblem>
loadBalancedChemistryModel<ReactionThermo, ThermoType>::get_problems(PtrList<volScalarField>& Y_,
                                                              const DeltaTType&            deltaT) {

    // TODO: Add refcell and Treact as conditions to get problems.

    const scalarField&            T = this->thermo().T();
    const scalarField&            p = this->thermo().p();
    tmp<volScalarField>           trho(this->thermo().rho());
    const scalarField&            rho          = trho();

    DynamicList<chemistryProblem> chem_problems;
    
    //TODO: reserve
    //chem_problems.reserve(T.size());


    forAll(p, celli) {


            for (label i = 0; i < this->nSpecie_; i++) { 
                this->c_[i] = rho[celli] * this->Y_[i][celli]/this->specieThermo_[i].W(); 
            }

            scalar Ti = T[celli];

            if (Ti > this->Treact()) {

                // Create the problem
                chemistryProblem problem;
                problem.c          = this->c_;
                problem.Ti         = T[celli];
                problem.pi         = p[celli];
                problem.rhoi       = rho[celli];
                problem.deltaTChem = this->deltaTChem_[celli];
                problem.deltaT     = deltaT[celli];
                problem.cpuTime    = this->cpu_times_[celli];
                problem.cellid     = celli;

                chem_problems.append(problem);
            }

            else {
                for (label i = 0; i < this->nSpecie(); i++) { this->RR_[i][celli] = 0; }
            }


    }
    return chem_problems;
}





} // namespace Foam
