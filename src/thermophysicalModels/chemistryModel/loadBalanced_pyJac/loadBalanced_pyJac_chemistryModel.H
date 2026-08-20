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
    Foam::chemistryModels::loadBalanced_pyJac

Description
    Extends chemistryModels::loadBalanced by adding an analytical Jacobian.

SourceFiles
    loadBalanced_pyJac_chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalanced_pyJac_chemistryModel_H
#define loadBalanced_pyJac_chemistryModel_H

#include "clockTime.H"
#include "loadBalanced_chemistryModel.H"

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

namespace chemistryModels
{

/*---------------------------------------------------------------------------*\
                     Class loadBalanced_pyJac Declaration
\*---------------------------------------------------------------------------*/

template <class ThermoType>
class loadBalanced_pyJac
    : public loadBalanced<ThermoType> {

public:
    //- Runtime type information
    TypeName("loadBalanced_pyJac");

    // Constructors

    //- Construct from thermo
    loadBalanced_pyJac(const fluidMulticomponentThermo& thermo);

    //- Destructor
    virtual ~loadBalanced_pyJac();

    //- Disallow default bitwise copy construction
    loadBalanced_pyJac(const loadBalanced_pyJac&);

    inline virtual label nEqns() const;

    //- Override jacobian and derivatives function from chemistryModel
    virtual void jacobian
    (
        const scalar t,
        const scalarField& TYp,
        const label li,
        scalarField& dTYpdt,
        scalarSquareMatrix& J
    ) const override;

    virtual void derivatives
    (
        const scalar t,
        const scalarField& TYp,
        const label li,
        scalarField& dTYpdt
    ) const override;

    //- Solve the ODE system
    virtual void solve
    (
        scalar& p,
        scalar& T,
        scalarField& c,
        const label li,
        scalar& deltaT,
        scalar& subDeltaT
    ) const override;

    //- Override Qdot function from chemistryModel for using pyJac's enthalpy of formation
    virtual tmp<volScalarField> Qdot() const override;

    //- Override tc function from chemistryModel for using calculating chemical time scale
    //- via pyJac
    virtual tmp<volScalarField> tc() const override;


private:
    //- Temporary concentration field
    mutable scalarField Y_;

    ///
    ///@brief List of enthalpy of formation from pyJac
    ///
    ///
    scalarList sp_enth_form;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "loadBalanced_pyJac_chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "loadBalanced_pyJac_chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Opts this model into the Standard reaction and reduction registrations of the
// chemistryModel dynamic code template, which tests
// "#if <type>_chemistryModel == Standard_chemistryModel"
#define loadBalanced_pyJac_chemistryModel 1

#endif

// ************************************************************************* //
