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
Base OF-dev file path : src/thermophysicalModels/chemistryModel/chemistryModel/StandardChemistryModel/StandardChemistryModel.C

\*---------------------------------------------------------------------------*/

#include "pyJacChemistryModel.H"
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

#include "clockTime.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::pyJacChemistryModel
(
    ReactionThermo& thermo
)
:
    BasicChemistryModel<ReactionThermo>(thermo),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(FWD_RATES), // reads from pyJac library
    Treact_
    (
        BasicChemistryModel<ReactionThermo>::template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_),
    sp_enth_form(nSpecie_),
    nActiveCells(0),
    chemCPUT(0.0)
{


    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.p().mesh(),
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
            )
        );
    }
    
    //TODO: this is probobably not correct
    const IOdictionary SOMEDICT_DIDNT_CHECK_IF_CORRECT
        (
            IOobject
            (
                thermo.phasePropertyName("chemistryProperties"),
                thermo.db().time().constant(),
                thermo.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
    

    refcell_mapper_ = new simpleRefMapping(SOMEDICT_DIDNT_CHECK_IF_CORRECT, thermo.composition());

    //load_balancer_ = new simpleLoadBalancing(SOMEDICT_DIDNT_CHECK_IF_CORRECT, thermo.composition());

    load_balancer_ = new simpleLoadBalancing();

    if (this->chemistry_)
    {


        //- Enthalpy of formation for all species
        // sp_enth_form_ "note underscore" is a tmp variable to avoid problems 
        // in pyJac function call.
        std::vector<double> sp_enth_form_(nSpecie_,0.0);
        //- Enthalpy of formation is taken from pyJac at T-standard     
        eval_h(298.15,sp_enth_form_.data());         
        for(int i=0; i<nSpecie_; i++)
        {
            sp_enth_form[i] = sp_enth_form_[i];
        }
    }

    if (refcell_mapper_->active())
    {
        Info<<"Reference cell mapping is active."<<endl;
    }

    if (load_balancer_->active())
    {
        //load_balancer_->print_parameters();
        if(Pstream::parRun())
        {
            Info<<"Load balancing is active and running on "<<Pstream::nProcs()<<" cores."<<endl;
        }
        else
        {
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

template<class ReactionThermo, class ThermoType>
Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::
~pyJacChemistryModel()
{
    //TODO: use a smart pointer so this can be removed!
    delete refcell_mapper_;
    delete load_balancer_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::omega
( 
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    notImplemented
    (
        "pyJacChemistryModel::omega"
        "("
            "scalarField&, "
            "scalar, "
            "scalar, "
            "scalarField& "
        ") const"
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    notImplemented
    (
        "pyJacChemistryModel::omegaI"
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
        ") const"
    );

    return(0);
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    notImplemented
    (
        "pyJacChemistryModel::omega"
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
        ") const"
    );

    return(0);
}


template<class ReactionThermo, class ThermoType>
void Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{

    //TODO: pre-allocate me!
    //the '()' operator at the end sets the chunk to zero
    std::vector<double> yToPyJac(nSpecie_+1,0.0);
    std::vector<double> dy(nSpecie_,0.0);

   

    const scalar T = c[0];
    const scalar p = c[nSpecie_];

    scalar csum = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        c_[i] = max(c[i+1], 0.0);
	    csum += c_[i];
    }
    c_[nSpecie_-1] = 1.0 - csum; // The last specie

    yToPyJac[0] = T;
    // i=1->nSpecie are mass fractions
    for (label i=1; i<nSpecie_; i++)
    {
   	    yToPyJac[i] = c_[i-1];
    }
    // The last specie
    yToPyJac[nSpecie_] = c_[nSpecie_-1];    

    // call pyJac RHS function
    dydt(0,p,yToPyJac.data(),dy.data());
    for (label i=0; i<nSpecie_; i++)
    {
    	dcdt[i] = dy[i];
    }
    // dp/dt = 0
    dcdt[nSpecie_] = 0.0;


}


template<class ReactionThermo, class ThermoType>
void Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    std::vector<double> yToPyJac(nSpecie_+1,0.0);
    std::vector<double> jac(nSpecie_*nSpecie_,0.0);

    J = Zero;
    dcdt = Zero;
    const scalar T = c[0];
    const scalar p = c[nSpecie_];

    scalar csum = 0.0;

    for (label i=0; i<nSpecie_; i++)
    {
        c_[i] = max(c[i+1], 0);
	    csum += c_[i];
    }
    c_[nSpecie_-1] = 1.0 - csum; // The last specie

    yToPyJac[0] = T;
    // i=1->nSpecie are mass fractions
    for (label i=1; i<nSpecie_; i++)
    {
   	    yToPyJac[i] = c_[i-1];
    }
    // The last specie
    yToPyJac[nSpecie_] = c_[nSpecie_-1];
    
    // call pyJac Jacobian evaluation
    eval_jacob (0, p, yToPyJac.data(), jac.data());

    int k=0;
    for (label j=0; j<nSpecie_; j++)
    {
        for (label i=0; i<nSpecie_; i++)
        {
            J[i][j] = jac[k+i];
        }
        k += nSpecie_;
    }

    // Last row and column to zero
    for (label j=0; j<nSpecie_+1; j++)
    {
        J[nSpecie_][j] = 0.0;
        J[j][nSpecie_] = 0.0;
    }

}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::tc() const
{
    notImplemented
    (
        "pyJacChemistryModel::tc() const"
    );
    
    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    return ttc;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                Qdot[celli] -= sp_enth_form[i]*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{

    notImplemented
    (
        "pyJacChemistryModel::calculateRR"
        "("
            "label, "
            "label "
        ") const"
    );
    
    tmp<volScalarField::Internal> tRR
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
        )
    );

    return tRR;
}


template<class ReactionThermo, class ThermoType>
void Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    notImplemented
    (
        "pyJacChemistryModel::calculateRR() const"
    );
}

template<class ReactionThermo, class ThermoType>
Foam::List<Foam::chemistrySolution> Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::solve_problems(Foam::DynamicList<Foam::chemistryProblem>& problems)
{
    List<chemistrySolution> chem_solns(problems.size());

    for (int i = 0; i < chem_solns.size(); i++)
    {
        chem_solns[i] = callODE(problems[i]);
    }
    return chem_solns;

}


template<class ReactionThermo, class ThermoType>
Foam::DynamicList<Foam::chemistryProblem> Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::get_problems(PtrList<volScalarField>& Y_,const scalar& deltaT)
{

    // TODO: Add refcell and Treact as conditions to get problems.
    DynamicList<chemistryProblem> chem_problems;
    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
   
    forAll(p,celli)
    {

       for (label i=0; i<nSpecie_; i++)
        {
            c_[i] = Y_[i][celli];
        }

        // Create the problem
        chemistryProblem problem;
        problem.pi = p[celli];
        problem.Ti = T[celli];
        problem.c = c_;
        problem.rhoi = rho[celli];
        problem.deltaTChem = this->deltaTChem_[celli];
        problem.cellid = celli;
        problem.deltaT = deltaT;
      
        if (refcell_mapper_->active())
        {
            if(!refcell_mapper_-> applyMapping(Y_,celli))
            {
                chem_problems.append(problem);
            }
        }
        else
        {
            chem_problems.append(problem);
        }
    }
    return chem_problems;
}


template<class ReactionThermo, class ThermoType>
Foam::chemistrySolution Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::callODE(Foam::chemistryProblem& prob)
{

    chemistrySolution soln;
    soln.cellid = prob.cellid;
    //- Initialize the time progress    
    scalar timeLeft = prob.deltaT;

    scalarList c0 = prob.c;
    // Calculate the chemical source terms
    while (timeLeft > small)
    {
        scalar dt = timeLeft;
       
        this->solve(prob.c, prob.Ti, prob.pi, dt, prob.deltaTChem);
        timeLeft -= dt;
    }
    scalarList RRtmp(this->nSpecie_,0.0);
    soln.c = prob.c;
    //RR_ should be weighted by density, due to NS formulation. rr = Y0-Y1 / dt, but RR_ = rho(Y0-Y1)/dt , the effect at Sh() function too!
    for (label i=0; i<this->nSpecie_; i++)
    {
        RRtmp[i] = prob.rhoi*(soln.c[i] - c0[i])/prob.deltaT; 
    }
    soln.RR = RRtmp;

    soln.deltaTChem = min(prob.deltaTChem, this->deltaTChemMax_);
    
    return soln;
}



template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }
    const scalar deltaT_ = deltaT[0];  // TODO: Check if this is true (all cells have same deltaT)

    // This assumes no refMapping, all the cells are in problems list
    DynamicList<chemistryProblem> problems = get_problems(Y_,deltaT_);

    List<chemistrySolution> solutions = solve_problems(problems);
    for (int i = 0; i < solutions.size(); i++)
    {
        label celli = solutions[i].cellid;
        for (label j=0; j<nSpecie_; j++)
        {
            RR_[j][celli] = solutions[i].RR[j];
        }
        deltaTMin = min(solutions[i].deltaTChem, deltaTMin);
        this->deltaTChem_[celli] = min(solutions[i].deltaTChem, this->deltaTChemMax_);   
    }

    return deltaTMin;


/*
    //- Create the reference solutions
    scalarList c_ref(nSpecie_,0.0);
    scalarList RR_ref(nSpecie_,0.0);

    bool refCellFound = false;

    //- Load Balancing
    if(load_balancer_->active())
    {
        //do stuff
    }
    //- CPU time analysis
    const clockTime clockTime_ = clockTime();
    clockTime_.timeIncrement();
    chemCPUT = 0.0;


    nActiveCells = 0;

    //- TODO: Call the loadcomputestats from load_balancer_
    forAll(rho, celli)
    {
        clockTime_.timeIncrement();

        scalar Ti = T[celli];

        if (Ti > Treact_)
        {

            const scalar rhoi = rho[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = Y_[i][celli];
                c0[i] = c_[i];
            }
            
            // Fill in the chemistry problem struct
            chemistryProblem prob;
            prob.pi = p[celli];
            prob.Ti = Ti;
            prob.c = c_;
            prob.rhoi = rhoi;
            prob.deltaTChem = this->deltaTChem_[celli];
            prob.cellid = celli;
            prob.deltaT = deltaT[celli];
            chemistrySolution soln;
            //- Refcell implementation
            if (refcell_mapper_->active())
            {
                //- First, try to find the first cell that can be calculated as a reference cell
                //- and map the solution to RR_ref and c_ref 

                if(!refCellFound)
                {
                    refCellFound = refcell_mapper_-> applyMapping(Y_,celli); 
                    if(refCellFound)
                    {
                        soln = callODE(prob);
                        for (label i=0; i<nSpecie_; i++)
                        {
                            RR_ref[i] = soln.RR[i];
                            c_ref[i] = soln.c[i];
                        }
                    }
                    soln = callODE(prob);
                    nActiveCells++;   

                }
                else if(refCellFound && refcell_mapper_-> applyMapping(Y_,celli))
                {
                    for (label i=0; i<nSpecie_; i++)
                    {
                        RR_[i][celli] = RR_ref[i];
                        c_[i] = c_ref[i];
                    }
                }
                else
                {
                    soln = callODE(prob);
                    nActiveCells++;
                }
            }
            else
            {
                soln = callODE(prob);
                nActiveCells++;

            }

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = soln.RR[i];
                c_[i] = soln.c[i];
            }
            deltaTMin = min(soln.deltaTChem, deltaTMin);
            this->deltaTChem_[celli] = min(soln.deltaTChem, this->deltaTChemMax_);
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }

        //- Update the CPU time spent on chemistry in the cell
        chemCPUT += clockTime_.timeIncrement();

    }
    Info<<"Number of active cells is : "<<nActiveCells<<endl;
*/
}



template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::solve
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


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::pyJacChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
