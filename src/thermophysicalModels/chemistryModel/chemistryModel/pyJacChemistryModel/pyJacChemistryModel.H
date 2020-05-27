/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::pyJacChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    pyJacChemistryModelI.H
    pyJacChemistryModel.C

Base OF-dev git commit : 09e6c93d4a7f2abb20243da013952c8a9369a9f2
Base OF-dev file path :
src/thermophysicalModels/chemistryModel/chemistryModel/StandardChemistryModel/StandardChemistryModel.H

\*---------------------------------------------------------------------------*/

#ifndef pyJacChemistryModel_H
#define pyJacChemistryModel_H

#ifdef DEBUG
#define OWN_DEBUG_VAR 1
#else
#define OWN_DEBUG_VAR 0
#endif


#include <memory> //std::unique_ptr
#include <vector> //std::vector

#include "BasicChemistryModel.H"
#include "ODESystem.H"
#include "Reaction.H"
#include "simpleMatrix.H"
#include "volFields.H"



#include "chemistryProblem.H"
#include "chemistrySolution.H"
#include "simpleLoadBalancing.H"
#include "bulutLoadBalancing.H"
#include "simpleRefMapping.H" //TODO: create a factory class that returns the correct ref cell mapper

// pyJac dependencies
extern "C" {
#include "chem_utils.h"
#include "dydt.h"
#include "jacob.h"
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class pyJacChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ReactionThermo, class ThermoType>
class pyJacChemistryModel : public BasicChemistryModel<ReactionThermo>, public ODESystem {
    // Private Member Functions

    //- Solve the reaction system for the given time step
    //  of given type and return the characteristic time
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);

    //- Disallow copy constructor
    pyJacChemistryModel(const pyJacChemistryModel&);

    //- Disallow default bitwise assignment
    void operator=(const pyJacChemistryModel&);

protected:
    typedef ThermoType thermoType; // what is this??

    // Ref cell mapper
    // autoPtr<chemistryRefMappingMethod<ReactionThermo, ThermoType>> refcell_mapper_;

    // TODO: DO NOT USE A RAW POINTER
    chemistryRefMappingMethod* refcell_mapper_;

    ///// Load balancing /////
    chemistryLoadBalancingMethod* load_balancer_;

    //- Number of active cells for load balance statistics
    label nActiveCells;

    //- Time consumed in chemistry evaluation
    scalar chemCPUT;

    // make_unique<Vec3>

    // Protected data

    //- Reference to the field of specie mass fractions
    PtrList<volScalarField>& Y_;

    //- Reactions
    const PtrList<Reaction<ThermoType>>& reactions_;

    //- Thermodynamic data of the species
    const PtrList<ThermoType>& specieThermo_;

    //- Number of species
    label nSpecie_;

    //- Number of reactions
    label nReaction_;

    //- Temperature below which the reaction rates are assumed 0
    scalar Treact_;

    //- List of reaction rate per specie [kg/m3/s]
    PtrList<volScalarField::Internal> RR_;

    //- Temporary concentration field
    mutable scalarField c_;

    //- Temporary rate-of-change of concentration field
    mutable scalarField dcdt_;

    //- List of enthalpies of formation
    scalarList sp_enth_form;

    // Protected Member Functions

    //- Write access to chemical source terms
    //  (e.g. for multi-chemistry model)
    inline PtrList<volScalarField::Internal>& RR();

public:
    //- Runtime type information
    TypeName("pyJac");

    // Constructors

    //- Construct from thermo
    pyJacChemistryModel(ReactionThermo& thermo);

    //- Destructor
    virtual ~pyJacChemistryModel();

    // Member Functions

    //- The reactions
    inline const PtrList<Reaction<ThermoType>>& reactions() const;

    //- Thermodynamic data of the species
    inline const PtrList<ThermoType>& specieThermo() const;

    //- The number of species
    virtual inline label nSpecie() const;

    //- The number of reactions
    virtual inline label nReaction() const;

    //- Temperature below which the reaction rates are assumed 0
    inline scalar Treact() const;

    //- Temperature below which the reaction rates are assumed 0
    inline scalar& Treact();

    //- dc/dt = omega, rate of change in concentration, for each species
    virtual void
    omega(const scalarField& c, const scalar T, const scalar p, scalarField& dcdt) const;

    //- Return the reaction rate for reaction r and the reference
    //  species and charateristic times
    virtual scalar omega(const Reaction<ThermoType>& r,
                         const scalarField&          c,
                         const scalar                T,
                         const scalar                p,
                         scalar&                     pf,
                         scalar&                     cf,
                         label&                      lRef,
                         scalar&                     pr,
                         scalar&                     cr,
                         label&                      rRef) const;

    //- Return the reaction rate for iReaction and the reference
    //  species and charateristic times
    virtual scalar omegaI(label              iReaction,
                          const scalarField& c,
                          const scalar       T,
                          const scalar       p,
                          scalar&            pf,
                          scalar&            cf,
                          label&             lRef,
                          scalar&            pr,
                          scalar&            cr,
                          label&             rRef) const;

    //- Calculates the reaction rates
    virtual void calculate();

    // Chemistry model functions (overriding abstract functions in
    // basicChemistryModel.H)

    //- Return const access to the chemical source terms for specie, i
    inline const volScalarField::Internal& RR(const label i) const;

    //- Return non const access to chemical source terms [kg/m3/s]
    virtual volScalarField::Internal& RR(const label i);

    //- Return reaction rate of the speciei in reactionI
    virtual tmp<volScalarField::Internal> calculateRR(const label reactionI,
                                                      const label speciei) const;

    //- Solve the reaction system for the given time step
    //  and return the characteristic time
    virtual scalar solve(const scalar deltaT);

    //- Solve the reaction system for the given time step
    //  and return the characteristic time
    virtual scalar solve(const scalarField& deltaT);

    //- Return the chemical time scale
    virtual tmp<volScalarField> tc() const;

    //- Return the heat release rate [kg/m/s3]
    virtual tmp<volScalarField> Qdot() const;

    // ODE functions (overriding abstract functions in ODE.H)

    //- Number of ODE's to solve
    inline virtual label nEqns() const;

    virtual void derivatives(const scalar t, const scalarField& c, scalarField& dcdt) const;

    virtual void
    jacobian(const scalar t, const scalarField& c, scalarField& dcdt, scalarSquareMatrix& J) const;

    virtual void
    solve(scalarField& c, scalar& T, scalar& p, scalar& deltaT, scalar& subDeltaT) const = 0;

private:
    chemistrySolution              solve_single(chemistryProblem& problem) const;
    DynamicList<chemistrySolution> solve_list(DynamicList<chemistryProblem>& problems) const;

    chemistryLoadBalancingMethod::buffer_t<chemistrySolution> solve_buffer(chemistryLoadBalancingMethod::buffer_t<chemistryProblem>& problems) const;






    DynamicList<chemistryProblem> get_problems(PtrList<volScalarField>& Y_, const scalar& deltaT);

    scalar update_reaction_rates(const DynamicList<chemistrySolution>& solutions);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pyJacChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "pyJacChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
