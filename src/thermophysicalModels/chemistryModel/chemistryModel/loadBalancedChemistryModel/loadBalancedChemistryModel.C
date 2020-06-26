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
    : StandardChemistryModel<ReactionThermo, ThermoType>(thermo) {

    Info << "loadBalancedChemistryModel: Number of species = " << this->nSpecie()
         << " and reactions = " << this->nReaction() << endl;

    load_balancer_ = new simpleBalancingMethod();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
loadBalancedChemistryModel<ReactionThermo, ThermoType>::~loadBalancedChemistryModel() {
    delete load_balancer_; // TODO: use a smart pointer
}

template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
scalar loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve(const DeltaTType& deltaT) {


    using problem_buffer_t  = chemistryLoadBalancingMethod::buffer_t<chemistryProblem>;
    using solution_buffer_t = chemistryLoadBalancingMethod::buffer_t<chemistrySolution>;


    Info << "HELLO FROM HERE" << endl;

    if (!this->chemistry_) { return great; }


    //TODO: What happens here???
    const scalar SOMEOTHER_DELTAT = deltaT[0]; 
    // This assumes no refMapping, all the cells are in problems list
    DynamicList<chemistryProblem> all_problems = get_problems(this->Y_, SOMEOTHER_DELTAT);


    int original_problem_count = all_problems.size();

    load_balancer_->update_state(all_problems);


    load_balancer_->print_state();
    
    problem_buffer_t balanced_problems = load_balancer_->balance(all_problems);


    solution_buffer_t balanced_solutions = solve_buffer(balanced_problems);

    solution_buffer_t my_solutions = load_balancer_->unbalance(balanced_solutions);




    int solution_count = 0;

    scalar deltaTMin = great;
    for (size_t i = 0; i < my_solutions.size(); ++i){

        scalar temp = update_reaction_rates(my_solutions[i]);
        if (temp < deltaTMin) {
            deltaTMin = temp;
        }
        solution_count += my_solutions[i].size(); 
    }



    if (solution_count != original_problem_count){
        throw error("Solution count differs from problem count.");
    }


    return deltaTMin;


    /*
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_) { return deltaTMin; }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField&  rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    auto& RR = this->RR();

    scalarField c0(this->nSpecie());

    forAll(rho, celli) {
        scalar Ti = T[celli];

        if (Ti > this->Treact()) {
            const scalar rhoi = rho[celli];
            scalar       pi   = p[celli];

            for (label i = 0; i < this->nSpecie(); i++) {
                this->c_[i] = rhoi * this->Y_[i][celli] / this->specieThermo_[i].W();
                c0[i]       = this->c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small) {
                scalar dt = timeLeft;
                this->solve(this->c_, Ti, pi, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] = min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i = 0; i < this->nSpecie(); i++) {
                RR[i][celli] = (this->c_[i] - c0[i]) * this->specieThermo_[i].W() / deltaT[celli];
            }
        } else {
            for (label i = 0; i < this->nSpecie(); i++) { RR[i][celli] = 0; }
        }
    }

    return deltaTMin;
    */
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
void loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve_single(
    chemistryProblem& prob, chemistrySolution& soln) const {

    // TODO: Make this work so that chemistryProblem is a const &

    scalar           timeLeft = prob.deltaT;
    const scalarList c0       = prob.c;

    // Calculate the chemical source terms
    while (timeLeft > small) {
        scalar dt = timeLeft;
        this->solve(prob.c, prob.Ti, prob.pi, dt, prob.deltaTChem);
        timeLeft -= dt;
    }

    soln.c = prob.c;
    soln.RR         = prob.rhoi * (soln.c - c0) / prob.deltaT;
    soln.cellid     = prob.cellid;
    soln.deltaTChem = min(prob.deltaTChem, this->deltaTChemMax_);
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
DynamicList<chemistryProblem>
loadBalancedChemistryModel<ReactionThermo, ThermoType>::get_problems(PtrList<volScalarField>& Y_,
                                                              const scalar&            deltaT) {

    // TODO: Add refcell and Treact as conditions to get problems.

    const scalarField&            T = this->thermo().T();
    const scalarField&            p = this->thermo().p();
    tmp<volScalarField>           trho(this->thermo().rho());
    const scalarField&            rho          = trho();
    bool                          refCellFound = false;

    DynamicList<chemistryProblem> chem_problems;
    chem_problems.reserve(T.size());

    chemistrySolution             ref_soln(this->nSpecie_);

    forAll(p, celli) {

        for (label i = 0; i < this->nSpecie_; i++) { this->c_[i] = this->Y_[i][celli]; }
        // Create the problem
        chemistryProblem problem;
        problem.pi         = p[celli];
        problem.Ti         = T[celli];
        problem.c          = this->c_;
        problem.rhoi       = rho[celli];
        problem.deltaTChem = this->deltaTChem_[celli];
        problem.cellid     = celli;
        problem.deltaT     = deltaT;


        chem_problems.append(problem);
    }
    return chem_problems;
}






template <class ReactionThermo, class ThermoType>
scalar loadBalancedChemistryModel<ReactionThermo, ThermoType>::update_reaction_rates(
    const DynamicList<chemistrySolution>& solutions) {

    scalar deltaTMin = great;
    for (int i = 0; i < solutions.size(); i++) {
        update_reaction_rate(solutions[i],solutions[i].cellid);
        deltaTMin                = min(solutions[i].deltaTChem, deltaTMin);
    }

    return deltaTMin;
}




template <class ReactionThermo, class ThermoType>
void loadBalancedChemistryModel<ReactionThermo, ThermoType>::update_reaction_rate(
    const chemistrySolution& solution, const label& cellid){

        for (label j = 0; j < this->nSpecie_; j++) { this->RR_[j][cellid] = solution.RR[j]; }
        this->deltaTChem_[cellid] = min(solution.deltaTChem, this->deltaTChemMax_);
}






} // namespace Foam
