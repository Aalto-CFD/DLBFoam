/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Base OF-dev git commit : 09e6c93d4a7f2abb20243da013952c8a9369a9f2
Base OF-dev file path :
src/thermophysicalModels/chemistryModel/chemistryModel/StandardChemistryModel/StandardChemistryModel.C

\*---------------------------------------------------------------------------*/

#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "pyJacChemistryModel.H"
#include "reactingMixture.H"

#include "clockTime.H"

namespace Foam {

template <class ReactionThermo, class ThermoType>
pyJacChemistryModel<ReactionThermo, ThermoType>::pyJacChemistryModel(ReactionThermo& thermo)
    : BasicChemistryModel<ReactionThermo>(thermo)
    , ODESystem()
    , Y_(this->thermo().composition().Y())
    , reactions_(dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo()))
    , specieThermo_(dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo()).speciesData())
    ,

    nSpecie_(Y_.size())
    , nReaction_(FWD_RATES)
    , // reads from pyJac library
    Treact_(BasicChemistryModel<ReactionThermo>::template lookupOrDefault<scalar>("Treact", 0))
    , RR_(nSpecie_)
    , c_(nSpecie_)
    , dcdt_(nSpecie_)
    , sp_enth_form(nSpecie_)
    , nActiveCells(0)
    , chemCPUT(0.0) {

    // Create the fields for the chemistry sources
    forAll(RR_, fieldi) {
        RR_.set(fieldi,
                new volScalarField::Internal(
                    IOobject("RR." + Y_[fieldi].name(),
                             this->mesh().time().timeName(),
                             this->mesh(),
                             IOobject::NO_READ,
                             IOobject::NO_WRITE),
                    thermo.p().mesh(),
                    dimensionedScalar("zero", dimMass / dimVolume / dimTime, 0)));
    }

    // TODO: this is probobably not correct
    const IOdictionary SOMEDICT_DIDNT_CHECK_IF_CORRECT(
        IOobject(thermo.phasePropertyName("chemistryProperties"),
                 thermo.db().time().constant(),
                 thermo.db(),
                 IOobject::MUST_READ,
                 IOobject::NO_WRITE,
                 false));

    refcell_mapper_ = new simpleRefMapping(SOMEDICT_DIDNT_CHECK_IF_CORRECT, thermo.composition());

    // load_balancer_ = new simpleLoadBalancing(SOMEDICT_DIDNT_CHECK_IF_CORRECT,
    // thermo.composition());

    load_balancer_ = new simpleLoadBalancing();
    //load_balancer_ = new bulutLoadBalancing();

    if (this->chemistry_) {

        //- Enthalpy of formation for all species
        // sp_enth_form_ "note underscore" is a tmp variable to avoid problems
        // in pyJac function call.
        std::vector<double> sp_enth_form_(nSpecie_, 0.0);
        //- Enthalpy of formation is taken from pyJac at T-standard
        eval_h(298.15, sp_enth_form_.data());
        for (int i = 0; i < nSpecie_; i++) { sp_enth_form[i] = sp_enth_form_[i]; }
    }

    if (refcell_mapper_->active()) { Info << "Reference cell mapping is active." << endl; }

    if (load_balancer_->active()) {
        // load_balancer_->print_parameters();
        if (Pstream::parRun()) {
            Info << "Load balancing is active and running on " << Pstream::nProcs() << " cores."
                 << endl;
        } else {
            /*
            FatalErrorIn
            (
                "chemistryModel::New"
                "(const fvMesh& mesh)"
            )   << "You are trying to do load balancing on a single core."
                << exit(FatalError);
            */
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
pyJacChemistryModel<ReactionThermo, ThermoType>::~pyJacChemistryModel() {
    // TODO: use a smart pointer so this can be removed!
    delete refcell_mapper_;
    delete load_balancer_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
void pyJacChemistryModel<ReactionThermo, ThermoType>::omega(const scalarField& c,
                                                            const scalar       T,
                                                            const scalar       p,
                                                            scalarField&       dcdt) const {
    notImplemented("pyJacChemistryModel::omega"
                   "("
                   "scalarField&, "
                   "scalar, "
                   "scalar, "
                   "scalarField& "
                   ") const");
}

template <class ReactionThermo, class ThermoType>
scalar pyJacChemistryModel<ReactionThermo, ThermoType>::omegaI(const label        index,
                                                               const scalarField& c,
                                                               const scalar       T,
                                                               const scalar       p,
                                                               scalar&            pf,
                                                               scalar&            cf,
                                                               label&             lRef,
                                                               scalar&            pr,
                                                               scalar&            cr,
                                                               label&             rRef) const {
    notImplemented("pyJacChemistryModel::omegaI"
                   "("
                   "label, "
                   "scalarField&, "
                   "scalar, "
                   "scalar, "
                   "scalar&, "
                   "scalar&, "
                   "label&, "
                   "scalar&, "
                   "scalar&, "
                   "label& "
                   ") const");

    return (0);
}

template <class ReactionThermo, class ThermoType>
scalar pyJacChemistryModel<ReactionThermo, ThermoType>::omega(const Reaction<ThermoType>& R,
                                                              const scalarField&          c,
                                                              const scalar                T,
                                                              const scalar                p,
                                                              scalar&                     pf,
                                                              scalar&                     cf,
                                                              label&                      lRef,
                                                              scalar&                     pr,
                                                              scalar&                     cr,
                                                              label& rRef) const {
    notImplemented("pyJacChemistryModel::omega"
                   "("
                   "Reaction<ThermoType>&, "
                   "scalarField&, "
                   "scalar, "
                   "scalar, "
                   "scalar&, "
                   "scalar&, "
                   "label&, "
                   "scalar&, "
                   "scalar&, "
                   "label& "
                   ") const");

    return (0);
}

template <class ReactionThermo, class ThermoType>
void pyJacChemistryModel<ReactionThermo, ThermoType>::derivatives(const scalar       time,
                                                                  const scalarField& c,
                                                                  scalarField&       dcdt) const {

    // TODO: pre-allocate me!
    // the '()' operator at the end sets the chunk to zero
    std::vector<double> yToPyJac(nSpecie_ + 1, 0.0);
    std::vector<double> dy(nSpecie_, 0.0);

    const scalar T = c[0];
    const scalar p = c[nSpecie_];

    scalar csum = 0.0;
    for (label i = 0; i < nSpecie_; i++) {
        c_[i] = max(c[i + 1], 0.0);
        csum += c_[i];
    }
    c_[nSpecie_ - 1] = 1.0 - csum; // The last specie

    yToPyJac[0] = T;
    // i=1->nSpecie are mass fractions
    for (label i = 1; i < nSpecie_; i++) { yToPyJac[i] = c_[i - 1]; }
    // The last specie
    yToPyJac[nSpecie_] = c_[nSpecie_ - 1];

    // call pyJac RHS function
    dydt(0, p, yToPyJac.data(), dy.data());
    for (label i = 0; i < nSpecie_; i++) { dcdt[i] = dy[i]; }
    // dp/dt = 0
    dcdt[nSpecie_] = 0.0;
}

template <class ReactionThermo, class ThermoType>
void pyJacChemistryModel<ReactionThermo, ThermoType>::jacobian(const scalar        t,
                                                               const scalarField&  c,
                                                               scalarField&        dcdt,
                                                               scalarSquareMatrix& J) const {
    std::vector<double> yToPyJac(nSpecie_ + 1, 0.0);
    std::vector<double> jac(nSpecie_ * nSpecie_, 0.0);

    J              = Zero;
    dcdt           = Zero;
    const scalar T = c[0];
    const scalar p = c[nSpecie_];

    scalar csum = 0.0;

    for (label i = 0; i < nSpecie_; i++) {
        c_[i] = max(c[i + 1], 0);
        csum += c_[i];
    }
    c_[nSpecie_ - 1] = 1.0 - csum; // The last specie

    yToPyJac[0] = T;
    // i=1->nSpecie are mass fractions
    for (label i = 1; i < nSpecie_; i++) { yToPyJac[i] = c_[i - 1]; }
    // The last specie
    yToPyJac[nSpecie_] = c_[nSpecie_ - 1];

    // call pyJac Jacobian evaluation
    eval_jacob(0, p, yToPyJac.data(), jac.data());

    int k = 0;
    for (label j = 0; j < nSpecie_; j++) {
        for (label i = 0; i < nSpecie_; i++) { J[i][j] = jac[k + i]; }
        k += nSpecie_;
    }

    // Last row and column to zero
    for (label j = 0; j < nSpecie_ + 1; j++) {
        J[nSpecie_][j] = 0.0;
        J[j][nSpecie_] = 0.0;
    }
}

template <class ReactionThermo, class ThermoType>
tmp<volScalarField> pyJacChemistryModel<ReactionThermo, ThermoType>::tc() const {
    notImplemented("pyJacChemistryModel::tc() const");

    tmp<volScalarField> ttc(new volScalarField(IOobject("tc",
                                                        this->time().timeName(),
                                                        this->mesh(),
                                                        IOobject::NO_READ,
                                                        IOobject::NO_WRITE,
                                                        false),
                                               this->mesh(),
                                               dimensionedScalar("zero", dimTime, small),
                                               extrapolatedCalculatedFvPatchScalarField::typeName));

    return ttc;
}

template <class ReactionThermo, class ThermoType>
tmp<volScalarField> pyJacChemistryModel<ReactionThermo, ThermoType>::Qdot() const {
    tmp<volScalarField> tQdot(
        new volScalarField(IOobject("Qdot",
                                    this->mesh_.time().timeName(),
                                    this->mesh_,
                                    IOobject::NO_READ,
                                    IOobject::NO_WRITE,
                                    false),
                           this->mesh_,
                           dimensionedScalar("zero", dimEnergy / dimVolume / dimTime, 0)));

    if (this->chemistry_) {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i) {
            forAll(Qdot, celli) { Qdot[celli] -= sp_enth_form[i] * RR_[i][celli]; }
        }
    }

    return tQdot;
}

template <class ReactionThermo, class ThermoType>
tmp<DimensionedField<scalar, volMesh>>
pyJacChemistryModel<ReactionThermo, ThermoType>::calculateRR(const label ri, const label si) const {

    notImplemented("pyJacChemistryModel::calculateRR"
                   "("
                   "label, "
                   "label "
                   ") const");

    tmp<volScalarField::Internal> tRR(
        new volScalarField::Internal(IOobject("RR",
                                              this->mesh().time().timeName(),
                                              this->mesh(),
                                              IOobject::NO_READ,
                                              IOobject::NO_WRITE),
                                     this->mesh(),
                                     dimensionedScalar("zero", dimMass / dimVolume / dimTime, 0)));

    return tRR;
}

template <class ReactionThermo, class ThermoType>
void pyJacChemistryModel<ReactionThermo, ThermoType>::calculate() {
    if (!this->chemistry_) { return; }

    notImplemented("pyJacChemistryModel::calculateRR() const");
}

template <class ReactionThermo, class ThermoType>
DynamicList<chemistryProblem>
pyJacChemistryModel<ReactionThermo, ThermoType>::get_problems(PtrList<volScalarField>& Y_,
                                                              const scalar&            deltaT) {

    // TODO: Add refcell and Treact as conditions to get problems.

    DynamicList<chemistryProblem> chem_problems;
    const scalarField&            T = this->thermo().T();
    const scalarField&            p = this->thermo().p();
    tmp<volScalarField>           trho(this->thermo().rho());
    const scalarField&            rho = trho();
    bool refCellFound = false;
    chemistrySolution ref_soln;

    forAll(p, celli) {   
        if (refcell_mapper_->active())
        {
            for (label i = 0; i < nSpecie_; i++) { c_[i] = Y_[i][celli]; }

            // Create the problem
            chemistryProblem problem;
            problem.pi         = p[celli];
            problem.Ti         = T[celli];
            problem.c          = c_;
            problem.rhoi       = rho[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.cellid     = celli;
            problem.deltaT     = deltaT;
         
            if(refcell_mapper_-> shouldMap(problem))
            {  
                 // TODO: Also update deltaTmin from reference solution
                if(!refCellFound)
                {
                    ref_soln = solve_single(problem);
                    update_reaction_rate(ref_soln,ref_soln.cellid);
                    refCellFound = true;
                }
                else
                {
                    update_reaction_rate(ref_soln,celli);
                }
            }
            else
            {
                chem_problems.append(problem);
            }
        }
        else
        {
                for (label i = 0; i < nSpecie_; i++) { c_[i] = Y_[i][celli]; }

                // Create the problem
                chemistryProblem problem;
                problem.pi         = p[celli];
                problem.Ti         = T[celli];
                problem.c          = c_;
                problem.rhoi       = rho[celli];
                problem.deltaTChem = this->deltaTChem_[celli];
                problem.cellid     = celli;
                problem.deltaT     = deltaT;

                chem_problems.append(problem); 
        }
    }
    return chem_problems;
}

template <class ReactionThermo, class ThermoType>
chemistrySolution
pyJacChemistryModel<ReactionThermo, ThermoType>::solve_single(chemistryProblem& prob) const {

    chemistrySolution soln;
    soln.cellid = prob.cellid;
    //- Initialize the time progress
    scalar timeLeft = prob.deltaT;

    scalarList c0 = prob.c;
    // Calculate the chemical source terms
    while (timeLeft > small) {
        scalar dt = timeLeft;

        this->solve(prob.c, prob.Ti, prob.pi, dt, prob.deltaTChem);
        timeLeft -= dt;
    }
    scalarList RRtmp(this->nSpecie_, 0.0);
    soln.c = prob.c;
    // RR_ should be weighted by density, due to NS formulation. rr = Y0-Y1 / dt, but RR_ =
    // rho(Y0-Y1)/dt , the effect at Sh() function too!
    for (label i = 0; i < this->nSpecie_; i++) {
        RRtmp[i] = prob.rhoi * (soln.c[i] - c0[i]) / prob.deltaT;
    }
    soln.RR = RRtmp;

    soln.deltaTChem = min(prob.deltaTChem, this->deltaTChemMax_);

    return soln;
}

template <class ReactionThermo, class ThermoType>
DynamicList<chemistrySolution> pyJacChemistryModel<ReactionThermo, ThermoType>::solve_list(
    DynamicList<chemistryProblem>& problems) const {

    DynamicList<chemistrySolution> chem_solns(problems.size(), chemistrySolution());

    for (size_t i = 0; i < chem_solns.size(); i++) { chem_solns[i] = solve_single(problems[i]); }
    return chem_solns;
}

template <class ReactionThermo, class ThermoType>
chemistryLoadBalancingMethod::buffer_t<chemistrySolution>
pyJacChemistryModel<ReactionThermo, ThermoType>::solve_buffer(
    chemistryLoadBalancingMethod::buffer_t<chemistryProblem>& problems) const {

    chemistryLoadBalancingMethod::buffer_t<chemistrySolution> ret;

    for (size_t i = 0; i < problems.size(); ++i) { ret.append(solve_list(problems[i])); }
    return ret;
}

template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
scalar pyJacChemistryModel<ReactionThermo, ThermoType>::solve(const DeltaTType& deltaT) {
    
    using problem_buffer_t  = chemistryLoadBalancingMethod::buffer_t<chemistryProblem>;
    using solution_buffer_t = chemistryLoadBalancingMethod::buffer_t<chemistrySolution>;

    clockTime timer;

    double interval = timer.timeIncrement();    
    BasicChemistryModel<ReactionThermo>::correct();
    Info << "correct() took: " << timer.timeIncrement() << endl;


    if (!this->chemistry_) { return great; }
    const scalar deltaT_ = deltaT[0]; // TODO: Check if this is true (all cells have same deltaT)

    

    
    timer.timeIncrement();
    // This assumes no refMapping, all the cells are in problems list
    DynamicList<chemistryProblem> all_problems = get_problems(Y_, deltaT_);
    Info << "get_problems() took: " << timer.timeIncrement() << endl;


    int original_problem_count = all_problems.size();

    timer.timeIncrement();
    load_balancer_->update_state(all_problems);
    Info << "update_state() took: " << timer.timeIncrement() << endl;


    load_balancer_->print_state();
    
    timer.timeIncrement();
    problem_buffer_t balanced_problems = load_balancer_->balance(all_problems);
    Info << "balance() took: " << timer.timeIncrement() << endl;


    timer.timeIncrement();
    solution_buffer_t balanced_solutions = solve_buffer(balanced_problems);
    Info << "solve_buffer() took: " << timer.timeIncrement() << endl;

    timer.timeIncrement();
    solution_buffer_t my_solutions = load_balancer_->unbalance(balanced_solutions);
    Info << "unbalance() took: " << timer.timeIncrement() << endl;



    timer.timeIncrement();

    int solution_count = 0;

    scalar deltaTMin = great;
    for (size_t i = 0; i < my_solutions.size(); ++i){

        scalar temp = update_reaction_rates(my_solutions[i]);
        if (temp < deltaTMin) {
            deltaTMin = temp;
        }
        solution_count += my_solutions[i].size(); 
    }

    Info << "update_reaction_rates() took: " << timer.timeIncrement() << endl;


    if (solution_count != original_problem_count){
        throw error("Solution count differs from problem count.");
    }


    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
scalar pyJacChemistryModel<ReactionThermo, ThermoType>::update_reaction_rates(
    const DynamicList<chemistrySolution>& solutions) {

    scalar deltaTMin = great;
    for (int i = 0; i < solutions.size(); i++) {
        update_reaction_rate(solutions[i],solutions[i].cellid);
        deltaTMin                = min(solutions[i].deltaTChem, deltaTMin);
    }

    return deltaTMin;
}

template <class ReactionThermo, class ThermoType>
void pyJacChemistryModel<ReactionThermo, ThermoType>::update_reaction_rate(
    const chemistrySolution& solution, const label& cellid){

        for (label j = 0; j < nSpecie_; j++) { RR_[j][cellid] = solution.RR[j]; }
        this->deltaTChem_[cellid] = min(solution.deltaTChem, this->deltaTChemMax_);
}




template <class ReactionThermo, class ThermoType>
scalar pyJacChemistryModel<ReactionThermo, ThermoType>::solve(const scalar deltaT) {
    // Don't allow the time-step to change more than a factor of 2
    return min(this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)), 2 * deltaT);
}

template <class ReactionThermo, class ThermoType>
scalar pyJacChemistryModel<ReactionThermo, ThermoType>::solve(const scalarField& deltaT) {
    return this->solve<scalarField>(deltaT);
}

} // namespace Foam

// ************************************************************************* //
