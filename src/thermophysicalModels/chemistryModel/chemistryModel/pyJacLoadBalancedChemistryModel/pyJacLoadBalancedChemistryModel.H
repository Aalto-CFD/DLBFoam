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
    Foam::pyJacLoadBalancedChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    pyJacLoadBalancedChemistryModelI.H
    pyJacLoadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef pyJacLoadBalancedChemistryModel_H
#define pyJacLoadBalancedChemistryModel_H

#include "clockTime.H"
#include "loadBalancedChemistryModel.H"

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
                      Class pyJacLoadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ReactionThermo, class ThermoType>
class pyJacLoadBalancedChemistryModel :
    public loadBalancedChemistryModel<ReactionThermo, ThermoType> {
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);

public:
    //- Runtime type information
    TypeName("loadBalanced_pyJac");

    // Constructors

    //- Construct from thermo
    pyJacLoadBalancedChemistryModel(const ReactionThermo& thermo);

    //- Destructor
    virtual ~pyJacLoadBalancedChemistryModel();

    //- Disallow default bitwise copy construction
    pyJacLoadBalancedChemistryModel(const pyJacLoadBalancedChemistryModel&);

    inline virtual label nEqns() const;

    //- Override jacobian and derivatives function from StandardChemistryModel

    virtual void
    jacobian(const scalar t, const scalarField& c, scalarField& dcdt, scalarSquareMatrix& J) const;

    virtual void derivatives(const scalar t, const scalarField& c, scalarField& dcdt) const;

    //- Override Qdot function from StandardChemistryModel for using pyJac's enthalpy of formation

    virtual tmp<volScalarField> Qdot() const;

private:
    ///
    ///@brief List of enthalpy of formation from pyJac
    ///
    ///
    scalarList sp_enth_form;

    ///
    ///@brief Compute c_ variable as mass fraction
    ///
    ///@param rho
    ///@param i
    ///@param celli
    ///@return scalar
    ///
    scalar compute_c(const scalar& rho, const label& i, const label& celli) const;

    ///
    ///@brief Compute reaction rate for mass fraction
    ///
    ///@param j
    ///@param solution
    ///@return scalar
    ///
    scalar compute_RR(const label& j, const chemistrySolution& solution) const;

    ///
    ///@brief Get the mass fraction field for refmapping
    ///
    ///@param problem
    ///@return scalarField
    ///
    scalarField get_mass_fraction(const chemistryProblem& problem) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "pyJacLoadBalancedChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "pyJacLoadBalancedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //