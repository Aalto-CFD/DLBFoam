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

#include "pyJacLoadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

template <class ReactionThermo, class ThermoType>
pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    pyJacLoadBalancedChemistryModel(const ReactionThermo& thermo)
    : LoadBalancedChemistryModel<ReactionThermo, ThermoType>(thermo),
      sp_enth_form(this->nSpecie_)
{

    if(this->chemistry_)
    {
        std::vector<scalar> sp_enth_form_(this->nSpecie_, 0.0);
        //- Enthalpy of formation is taken from pyJac at T-standard
        eval_h(298.15, sp_enth_form_.data());
        for(label i = 0; i < this->nSpecie_; i++)
        {
            sp_enth_form[i] = sp_enth_form_[i];
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    ~pyJacLoadBalancedChemistryModel()
{
}

template <class ReactionThermo, class ThermoType>
void pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::jacobian(
    const scalar        t,
    const scalarField&  c,
    const label         li,
    scalarField&        dcdt,
    scalarSquareMatrix& J) const
{
    std::vector<scalar> yToPyJac(this->nSpecie_ + 1, 0.0);
    std::vector<scalar> jac(this->nSpecie_ * this->nSpecie_, 0.0);

    J                 = Zero;
    dcdt              = Zero;
    const scalar T    = c[0];
    const scalar p    = c[this->nSpecie_];
    scalar       csum = 0.0;

    for(label i = 0; i < this->nSpecie_ - 1; i++)
    {
        this->c_[i] = max(c[i + 1], 0);
        csum += this->c_[i];
    }
    this->c_[this->nSpecie_ - 1] = 1.0 - csum; // The last specie
    yToPyJac[0]                  = T;
    // i=1->nSpecie are mass fractions
    for(label i = 1; i < this->nSpecie_; i++)
    {
        yToPyJac[i] = this->c_[i - 1];
    }
    // The last specie

    yToPyJac[this->nSpecie_] = this->c_[this->nSpecie_ - 1];
    // call pyJac Jacobian evaluation
    eval_jacob(0, p, yToPyJac.data(), jac.data());
    label k = 0;
    for(label j = 0; j < this->nSpecie_; j++)
    {
        for(label i = 0; i < this->nSpecie_; i++)
        {
            J[i][j] = jac[k + i];
        }
        k += this->nSpecie_;
    }

    // Last row and column to zero
    for(label j = 0; j < this->nSpecie_ + 1; j++)
    {
        J[this->nSpecie_][j] = 0.0;
        J[j][this->nSpecie_] = 0.0;
    }
}

template <class ReactionThermo, class ThermoType>
void pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::derivatives(
    const scalar       t,
    const scalarField& c,
    const label        li,
    scalarField&       dcdt) const
{

    std::vector<scalar> yToPyJac(this->nSpecie_ + 1, 0.0);
    std::vector<scalar> dy(this->nSpecie_, 0.0);

    const scalar T    = c[0];
    const scalar p    = c[this->nSpecie_];
    scalar       csum = 0.0;
    for(label i = 0; i < this->nSpecie_ - 1; i++)
    {
        this->c_[i] = max(c[i + 1], 0.0);
        csum += this->c_[i];
    }
    this->c_[this->nSpecie_ - 1] = 1.0 - csum; // The last specie

    yToPyJac[0] = T;
    // i=1->nSpecie are mass fractions
    for(label i = 1; i < this->nSpecie_; i++)
    {
        yToPyJac[i] = this->c_[i - 1];
    }
    // The last specie
    yToPyJac[this->nSpecie_] = this->c_[this->nSpecie_ - 1];

    // call pyJac RHS function
    dydt(0, p, yToPyJac.data(), dy.data());
    for(label i = 0; i < this->nSpecie_; i++)
    {
        dcdt[i] = dy[i];
    }
    // dp/dt = 0
    dcdt[this->nSpecie_] = 0.0;
}

template <class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::Qdot() const
{

    tmp<volScalarField> tQdot(new volScalarField(
        IOobject(
            "Qdot",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false),
        this->mesh_,
        dimensionedScalar("zero", dimEnergy / dimVolume / dimTime, 0)));

    if(this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(this->Y_, i)
        {
            forAll(Qdot, celli)
            {
                Qdot[celli] -= sp_enth_form[i] * this->RR_[i][celli];
            }
        }
    }

    return tQdot;
}

template <class ReactionThermo, class ThermoType>
scalar pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    computeConcentration(
        const scalar& rho, const label& i, const label& celli) const
{

    return (this->Y_[i][celli]);
}

template <class ReactionThermo, class ThermoType>
scalar pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    computeReactionRate(const label& j, const ChemistrySolution& solution) const
{

    return (solution.rhoi * solution.c_increment[j]);
}

template <class ReactionThermo, class ThermoType>
ChemistryProblem
pyJacLoadBalancedChemistryModel<ReactionThermo, ThermoType>::getMassFraction(
    const ChemistryProblem& problem) const
{

    return problem;
}

} // namespace Foam
