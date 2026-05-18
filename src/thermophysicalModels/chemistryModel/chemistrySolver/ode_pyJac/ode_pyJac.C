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
    odeSolver_(ODESolver::New(*this, this->typeDict("ode"))),
    YTp_(this->nEqns())
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
    scalarField& Y,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(YTp_);
    }

    const label nSpecie = this->nSpecie();
    const label nEqns = this->nEqns();

    if(mag(nEqns-nSpecie) == 2) // original implementation
    {
        // Copy the concentration, T and P to the total solve-vector
        for (label i=0; i<nSpecie; i++)
        {
            YTp_[i] = Y[i];
        }
        YTp_[nSpecie] = T;
        YTp_[nSpecie+1] = p;

        odeSolver_->solve(0, deltaT, YTp_, li, subDeltaT);

        for (label i=0; i<nSpecie; i++)
        {
            Y[i] = max(0.0, YTp_[i]);
        }
        T = YTp_[nSpecie];
        p = YTp_[nSpecie+1];
    }
    /* pyJac implementation:                                                            */
    /* Because pyJac does not consider the last specie in the system, it is not part of */
    /* the RHS and Jacobian and thus we need here the inert = 1.0-csum functionality    */
    else
    {
        // Copy the concentration, T and P to the total solve-vector
        YTp_[0] = T;
        YTp_[nSpecie] = p;

        for (label i=0; i<nSpecie-1; i++)
        {
            YTp_[i+1] = Y[i];
        }

        odeSolver_->solve(0, deltaT, YTp_, li, subDeltaT);

        T = YTp_[0];
        p = YTp_[nSpecie];
        scalar Ysum = 0;

        for (label i=0; i<nSpecie-1; i++)
        {
            Y[i] = max(0.0, YTp_[i+1]);
    	    Ysum += Y[i];
        }
        //The last specie:
        Y[nSpecie-1] = 1.0 - Ysum;
    }
}


// ************************************************************************* //
