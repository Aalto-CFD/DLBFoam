/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           |
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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
    Foam::loadBalanced_pyJacChemistryModel

Description
    Extends loadBalancedChemistryModel model by adding analytical Jacobian.

SourceFiles
    loadBalanced_pyJacChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalanced_pyJacChemistryModel_H
#define loadBalanced_pyJacChemistryModel_H

#include "clockTime.H"
#include "loadBalancedChemistryModel.H"

// pyJac dependencies
extern "C" {
    #include "pyJacInclude/chem_utils.h"
    #include "pyJacInclude/dydt.h"
    #include "pyJacInclude/jacob.h"
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class loadBalanced_pyJacChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ThermoType>
class loadBalanced_pyJacChemistryModel
    : public loadBalancedChemistryModel<ThermoType> {
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);

public:
    //- Runtime type information
    TypeName("loadBalanced_pyJac");

    // Constructors

    //- Construct from thermo
    loadBalanced_pyJacChemistryModel(const fluidReactionThermo& thermo);

    //- Destructor
    virtual ~loadBalanced_pyJacChemistryModel();

    //- Disallow default bitwise copy construction
    loadBalanced_pyJacChemistryModel(const loadBalanced_pyJacChemistryModel&);

    inline virtual label nEqns() const;

    //- Override jacobian and derivatives function from chemistryModel
    virtual void jacobian(const scalar        t,
                          const scalarField&  c,
                          const label         li,
                          scalarField&        dcdt,
                          scalarSquareMatrix& J) const override;

    virtual void derivatives(const scalar       t,
                             const scalarField& c,
                             const label        li,
                             scalarField&       dcdt) const override;

    //- Override Qdot function from chemistryModel for using pyJac's enthalpy of formation
    virtual tmp<volScalarField> Qdot() const override;

    //- Override tc function from chemistryModel for using calculating chemical time scale
    //- via pyJac
    virtual tmp<volScalarField> tc() const override;


private:
    //- Temporary concentration field
    mutable scalarField c_;

    ///
    ///@brief List of enthalpy of formation from pyJac
    ///
    ///
    scalarList sp_enth_form;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "loadBalanced_pyJacChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "loadBalanced_pyJacChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
