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

\*---------------------------------------------------------------------------*/

#include "ode_pyJac.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ode_pyJac<ChemistryModel>::ode_pyJac(const fluidMulticomponentThermo& thermo)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("odeCoeffs")),
    odeSolver_(ODESolver::New(*this, coeffsDict_)),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::ode_pyJac<ChemistryModel>::~ode_pyJac()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::ode_pyJac<ChemistryModel>::solve
(
    scalar& p,
    scalar& T,
    scalarField& c,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(cTp_);
    }

    const label nSpecie = this->nSpecie();
    const label nEqns = this->nEqns();

    if(mag(nEqns-nSpecie) == 2) // original implementation
    {
        // Copy the concentration, T and P to the total solve-vector
        for (label i=0; i<nSpecie; i++)
        {
            cTp_[i] = c[i];
        }
        cTp_[nSpecie] = T;
        cTp_[nSpecie+1] = p;

        odeSolver_->solve(0, deltaT, cTp_, li, subDeltaT);

        for (label i=0; i<nSpecie; i++)
        {
            c[i] = max(0.0, cTp_[i]);
        }
        T = cTp_[nSpecie];
        p = cTp_[nSpecie+1];
    }
    /* pyJac implementation:                                                            */
    /* Because pyJac does not consider the last specie in the system, it is not part of */
    /* the RHS and Jacobian and thus we need here the inert = 1.0-csum functionality    */
    else
    {
        // Copy the concentration, T and P to the total solve-vector
        cTp_[0] = T;
        cTp_[nSpecie] = p;

        for (label i=0; i<nSpecie-1; i++)
        {
            cTp_[i+1] = c[i];
        }

        odeSolver_->solve(0, deltaT, cTp_, li, subDeltaT);

        T = cTp_[0];
        p = cTp_[nSpecie];
        scalar csum = 0;

        for (label i=0; i<nSpecie-1; i++)
        {
            c[i] = max(0.0, cTp_[i+1]);
    	    csum += c[i];
        }
        //The last specie:
        c[nSpecie-1] = 1.0 - csum;
    }
}


// ************************************************************************* //