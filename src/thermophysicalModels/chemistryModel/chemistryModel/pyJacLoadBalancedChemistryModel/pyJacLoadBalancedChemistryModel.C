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

#include "pyJacLoadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {

template <class ThermoType>
pyJacLoadBalancedChemistryModel<ThermoType>::pyJacLoadBalancedChemistryModel(
    const fluidReactionThermo& thermo)
    : loadBalancedChemistryModel<ThermoType>(thermo)
    , sp_enth_form(this->nSpecie_) {

    if (this->chemistry_) {
        // TODO: prevent symbol look-up error in case of ill mechanism library compilation or wrong path
        // 1) Instead of providing libc_pyjac.so in controlDict, it should be given in chemistryProperties as a lib() argument (similar to functionObjects).
        // 2) Read the new lib() as dictionary path variable here in the constructor.
        // 3) Implement here "is_pyjac_lib_available(pyjac_lib_path)":
        //      - utilise dlopen for the test, see e.g. https://stackoverflow.com/questions/56747328/loading-shared-library-dynamically-using-dlopen
        // 4) If library is not available, safe exit and print out "check your libc_pyjac.so path in chemistryProperties."

        //- Enthalpy of formation is taken from pyJac at T-standard
        std::vector<scalar> sp_enth_form_(this->nSpecie_, 0.0);
        eval_h(298.15, sp_enth_form_.data());
        for (label i = 0; i < this->nSpecie_; i++) { sp_enth_form[i] = sp_enth_form_[i]; }
    }

    Info << "Overriding standardChemistryModel by pyJacLoadBalancedChemistryModel:" << endl;

    if (this->nSpecie_ == PYJAC_NSP())
    {
        Info << "pyJac mechanism information:" <<
                "\n\tNumber of species: " << PYJAC_NSP() <<
                "\n\tNumber of forward reactions: " << PYJAC_FWD_RATES() << "\n" << endl;
    }
    else
    {
        FatalErrorIn
        (
            "pyJacLoadBalancedChemistryModel::New"
        )   << "\nInconsistent definition of number of species between thermophysicalProperties (Nsp = " << this->nSpecie_ << ") and pyJac library (Nsp = " << PYJAC_NSP() << ")"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ThermoType>
pyJacLoadBalancedChemistryModel<ThermoType>::~pyJacLoadBalancedChemistryModel() {}

template <class ThermoType>
void pyJacLoadBalancedChemistryModel<ThermoType>::jacobian(
    const scalar t, const scalarField& c, const label li, scalarField& dcdt, scalarSquareMatrix& J)
    const {

    std::vector<scalar> yToPyJac(this->nSpecie_ + 1, 0.0);
    std::vector<scalar> jac(this->nSpecie_ * this->nSpecie_, 0.0);

    J = Zero;
    dcdt = Zero;
    const scalar T = c[0];
    const scalar p = c[this->nSpecie_];
    scalar csum = 0.0;

    for (label i = 0; i < this->nSpecie_ - 1; i++) {
        this->c_[i] = max(c[i + 1], 0);
        csum += this->c_[i];
    }

    this->c_[this->nSpecie_ - 1] = 1.0 - csum; // The last specie
    yToPyJac[0]                  = T;
    // i=1->nSpecie are mass fractions
    for (label i = 1; i < this->nSpecie_; i++) { yToPyJac[i] = this->c_[i - 1]; }
    // The last specie

    yToPyJac[this->nSpecie_] = this->c_[this->nSpecie_ - 1];
    // call pyJac Jacobian evaluation
    eval_jacob(0, p, yToPyJac.data(), jac.data());
    label k = 0;
    for (label j = 0; j < this->nSpecie_; j++) {
        for (label i = 0; i < this->nSpecie_; i++) { J[i][j] = jac[k + i]; }
        k += this->nSpecie_;
    }

    // Last row and column to zero
    for (label j = 0; j < this->nSpecie_ + 1; j++) {
        J[this->nSpecie_][j] = 0.0;
        J[j][this->nSpecie_] = 0.0;
    }

}

template <class ThermoType>
void pyJacLoadBalancedChemistryModel<ThermoType>::derivatives(
    const scalar t, const scalarField& c, const label li, scalarField& dcdt) const {

    std::vector<scalar> yToPyJac(this->nSpecie_ + 1, 0.0);
    std::vector<scalar> dy(this->nSpecie_, 0.0);

    const scalar T = c[0];
    const scalar p = c[this->nSpecie_];
    scalar csum = 0.0;
    for (label i = 0; i < this->nSpecie_ - 1; i++) {
        this->c_[i] = max(c[i + 1], 0.0);
        csum += this->c_[i];
    }
    this->c_[this->nSpecie_ - 1] = 1.0 - csum; // The last specie

    yToPyJac[0] = T;
    // i=1->nSpecie are mass fractions
    for (label i = 1; i < this->nSpecie_; i++) { yToPyJac[i] = this->c_[i - 1]; }
    // The last specie
    yToPyJac[this->nSpecie_] = this->c_[this->nSpecie_ - 1];

    // call pyJac RHS function
    dydt(0, p, yToPyJac.data(), dy.data());
    for (label i = 0; i < this->nSpecie_; i++) { dcdt[i] = dy[i]; }
    // dp/dt = 0
    dcdt[this->nSpecie_] = 0.0;
}

template <class ThermoType>
Foam::tmp<Foam::volScalarField>
pyJacLoadBalancedChemistryModel<ThermoType>::Qdot() const {

    tmp<volScalarField> tQdot(
        new volScalarField(IOobject("Qdot",
                                    this->mesh_.time().timeName(),
                                    this->mesh_,
                                    IOobject::NO_READ,
                                    IOobject::NO_WRITE,
                                    false),
                           this->mesh_,
                           dimensionedScalar("zero", dimEnergy / dimVolume / dimTime, 0)));

    if (this->chemistry_) {
        scalarField& Qdot = tQdot.ref();

        forAll(this->Y_, i) {
            forAll(Qdot, celli) { Qdot[celli] -= sp_enth_form[i] * this->RR_[i][celli]; }
        }
    }

    return tQdot;
}

// TODO: Make this work for pyjac
template <class ThermoType>
void Foam::pyJacLoadBalancedChemistryModel<ThermoType>::updateReactionRate
(
    const ChemistrySolution& solution, const label& i
)
{
    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][i] = solution.c_increment[j] * solution.rhoi;
    }
    this->deltaTChem_[i] = min(solution.deltaTChem, this->deltaTChemMax_);
}

template <class ThermoType>
Foam::scalarField Foam::pyJacLoadBalancedChemistryModel<ThermoType>::getVariable
(
    const scalarField& concentration, const scalarField& massFraction
)
{
    return massFraction;
}


} // namespace Foam
