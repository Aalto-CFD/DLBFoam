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
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::loadBalancedChemistryModel<ReactionThermo, ThermoType>::loadBalancedChemistryModel
(
    ReactionThermo& thermo
)
//: aaltoChemistryModelBase<ReactionThermo, ThermoType>(thermo)
: StandardChemistryModel<ReactionThermo, ThermoType>(thermo)
{
   // Info<< "loadBalancedChemistryModel: Number of species = " << nSpecie_
   //     << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::loadBalancedChemistryModel<ReactionThermo, ThermoType>::
~loadBalancedChemistryModel()
{}





template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
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



    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    auto& RR = this->RR();

    scalarField c0(this->nSpecie());

    
    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        
        if (Ti > this->Treact())
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];

            for (label i=0; i<this->nSpecie(); i++)
            {
                this->c_[i] = rhoi*this->Y_[i][celli]/this->specieThermo_[i].W();
                c0[i] = this->c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(this->c_, Ti, pi, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<this->nSpecie(); i++)
            {
                RR[i][celli] =
                    (this->c_[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
            }
        }
        else
        {
            for (label i=0; i<this->nSpecie(); i++)
            {
                RR[i][celli] = 0;
            }
        }
        
    }
    
    
    return deltaTMin;
}

template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::loadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}



// ************************************************************************* //