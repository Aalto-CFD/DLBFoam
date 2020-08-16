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

Class
    Foam::loadBalancedChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    loadBalancedChemistryModelI.H
    loadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalancedChemistryModel_H
#define loadBalancedChemistryModel_H

#include "OFstream.H"
#include "StandardChemistryModel.H"
#include "clockTime.H"

#include "bulutLoadBalancing.H"
#include "chemistryProblem.H"
#include "chemistrySolution.H"
#include "globalBalancingMethod.H"
#include "simpleRefMappingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class loadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ReactionThermo, class ThermoType>
class loadBalancedChemistryModel :
    // public aaltoChemistryModelBase<ReactionThermo, ThermoType>
    public StandardChemistryModel<ReactionThermo, ThermoType> {
    // Private Member Functions

    //- Solve the reaction system for the given time step
    //  of given type and return the characteristic time
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);

    //- Disallow copy constructor
    loadBalancedChemistryModel(const loadBalancedChemistryModel&);

    //- Disallow default bitwise assignment
    void operator=(const loadBalancedChemistryModel&);

public:
    //- Runtime type information
    TypeName("loadBalanced");

    // Constructors

    //- Construct from thermo
    loadBalancedChemistryModel(ReactionThermo& thermo);

    //- Destructor
    virtual ~loadBalancedChemistryModel();

    virtual scalar solve(const scalar deltaT) override;
    virtual scalar solve(const scalarField& deltaT) override;

    ///
    ///@brief
    ///
    ///@param prob
    ///@param soln
    ///
    void solve_single(chemistryProblem& prob, chemistrySolution& soln) const;

    // ODE functions (overriding abstract functions in ODE.H)
    // this function is somehow injected in the ODE
    virtual void
    solve(scalarField& c, scalar& T, scalar& p, scalar& deltaT, scalar& subDeltaT) const = 0;

    //- Log file for average time spent solving the chemistry
    autoPtr<OFstream> cpuSolveFile_;
    //- Create and return a loadbal log file of the given name
    inline autoPtr<OFstream> logFile(const word& name) const;

private:
    scalarField cpu_times_;

    chemistryLoadBalancingMethod* load_balancer_;

    chemistryRefMappingMethod* ref_mapper_;

    ///
    ///@brief A factory function for creating the balancer
    ///
    ///@return chemistryLoadBalancingMethod*
    ///
    static chemistryLoadBalancingMethod* create_balancer();

    ///
    ///@brief A factory function for creating the refmapper
    ///
    ///@return chemistryRefMappingMethod*
    ///
    chemistryRefMappingMethod* create_refmapper();

    ///
    ///@brief
    ///
    ///@param problems
    ///@return chemistryLoadBalancingMethod::buffer_t<chemistrySolution>
    ///
    chemistryLoadBalancingMethod::buffer_t<chemistrySolution>
    solve_buffer(chemistryLoadBalancingMethod::buffer_t<chemistryProblem>& problems) const;

    ///
    ///@brief Get the problems object
    ///
    ///@param Y_
    ///@param deltaT
    ///@return DynamicList<chemistryProblem>
    ///

    template <class DeltaTType>
    DynamicList<chemistryProblem> get_problems(PtrList<volScalarField>& Y_,
                                               const DeltaTType&        deltaT);

    ///
    ///@brief
    ///
    ///@param solutions
    ///@return scalar
    ///
    scalar update_reaction_rates(
        const chemistryLoadBalancingMethod::buffer_t<chemistrySolution>& solutions);

    ///
    ///@brief updates a single problem's reaction rate for index i
    ///
    ///
    virtual void update_reaction_rate(const chemistrySolution& solution, const label& i);

    ///
    ///@brief Compute c_ variable as concentration
    ///
    ///@param rho
    ///@param i
    ///@param celli
    ///@return double
    ///
    virtual double compute_c(const scalar& rho, const label& i, const label& celli) const;

    ///
    ///@brief Compute reaction rate for concentration
    ///
    ///@param j
    ///@param solution
    ///@return double
    ///
    virtual double compute_RR(const label& j, const chemistrySolution& solution) const;

    ///
    ///@brief Get the mass fraction field for refmapping
    ///
    ///@param problem
    ///@return scalarField
    ///
    virtual scalarField get_mass_fraction(const chemistryProblem& problem) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "loadBalancedChemistryModelI.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "loadBalancedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //