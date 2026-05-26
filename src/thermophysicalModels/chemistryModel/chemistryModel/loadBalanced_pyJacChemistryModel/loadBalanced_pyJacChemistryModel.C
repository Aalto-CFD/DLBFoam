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

#include "loadBalanced_pyJacChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam {

template <class ThermoType>
loadBalanced_pyJacChemistryModel<ThermoType>::loadBalanced_pyJacChemistryModel(
    const fluidMulticomponentThermo& thermo)
    : loadBalancedChemistryModel<ThermoType>(thermo),
    Y_(this->nSpecie()),
    sp_enth_form(this->nSpecie()) {

    if (this->chemistry_) {
        // TODO: prevent symbol look-up error in case of ill mechanism library compilation or wrong path
        // 1) Instead of providing libc_pyjac.so in controlDict, it should be given in chemistryProperties as a lib() argument (similar to functionObjects).
        // 2) Read the new lib() as dictionary path variable here in the constructor.
        // 3) Implement here "is_pyjac_lib_available(pyjac_lib_path)":
        //      - utilise dlopen for the test, see e.g. https://stackoverflow.com/questions/56747328/loading-shared-library-dynamically-using-dlopen
        // 4) If library is not available, safe exit and print out "check your libc_pyjac.so path in chemistryProperties."

        //- Enthalpy of formation is taken from pyJac at T-standard
        std::vector<scalar> sp_enth_form_(this->nSpecie(), 0.0);
        eval_h(298.15, sp_enth_form_.data());
        for (label i = 0; i < this->nSpecie(); i++) { sp_enth_form[i] = sp_enth_form_[i]; }
    }

    Info << "Overriding chemistryModel by loadBalanced_pyJacChemistryModel:" << endl;

    if (this->nSpecie() == PYJAC_NSP())
    {
        Info << "pyJac mechanism information:" <<
                "\n\tNumber of species: " << PYJAC_NSP() <<
                "\n\tNumber of forward reactions: " << PYJAC_FWD_RATES() << "\n" << endl;
    }
    else
    {
        FatalErrorIn
        (
            "loadBalanced_pyJacChemistryModel::New"
        )   << "\nInconsistent definition of number of species between thermophysicalProperties (Nsp = " << this->nSpecie() << ") and pyJac library (Nsp = " << PYJAC_NSP() << ")"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ThermoType>
loadBalanced_pyJacChemistryModel<ThermoType>::~loadBalanced_pyJacChemistryModel() {}

template <class ThermoType>
void loadBalanced_pyJacChemistryModel<ThermoType>::jacobian(
    const scalar t, const scalarField& TYp, const label li, scalarField& dTYdt, scalarSquareMatrix& J)
    const {

    scalarField yToPyJac(this->nSpecie(), 0.0);
    scalarField jac(this->nSpecie() * this->nSpecie(), 0.0);

    J = Zero;
    dTYdt = Zero;
    const scalar p = TYp[this->nSpecie()];

    yToPyJac[0] = TYp[0];
    for (label i = 1; i < this->nSpecie(); i++) {
        yToPyJac[i] = max(TYp[i], 0);
    }

    eval_jacob(0, p, yToPyJac.begin(), jac.begin());

    for (label j = 0; j < this->nSpecie(); j++) {
        for (label i = 0; i < this->nSpecie(); i++) {
        J[i][j] = jac[i + j*this->nSpecie()];
        }
    }
}

template <class ThermoType>
void loadBalanced_pyJacChemistryModel<ThermoType>::derivatives(
    const scalar t, const scalarField& TYp, const label li, scalarField& dTYpdt) const {

    scalarField yToPyJac(this->nSpecie(), 0.0);

    const scalar p = TYp[this->nSpecie()];

    yToPyJac[0] = TYp[0];
    for (label i = 1; i < this->nSpecie(); i++) {
        yToPyJac[i] = max(TYp[i], 0);
    }

    // call pyJac RHS function
    dydt(0, p, yToPyJac.begin(), dTYpdt.begin());

    // dp/dt = 0
    dTYpdt[this->nSpecie()] = 0.0;
}

template <class ThermoType>
Foam::tmp<Foam::volScalarField>
loadBalanced_pyJacChemistryModel<ThermoType>::Qdot() const {

    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_) {
    scalarField& Qdot = tQdot.ref();

    forAll(this->Y(), i) {
        forAll(Qdot, celli) { Qdot[celli] -= sp_enth_form[i] * this->RR(i)[celli]; }
    }
    }

    return tQdot;
}

template <class ThermoType>
Foam::tmp<Foam::volScalarField>
loadBalanced_pyJacChemistryModel<ThermoType>::tc() const {

    if(PYJAC_FWD_RATES()!=this->nReaction()) {
    FatalErrorInFunction
        << "Number of reactions in pyJac (" << PYJAC_FWD_RATES() <<
        ") and in OpenFOAM's native chemistry model (" << this->nReaction() <<
        ") do not match.\nChemical time scale required in combustion models such "<<
        "as PaSR or EDC is calculated using standard chemistry model \n" <<
        "and thus reactions needs to be given for the native chemistry model as a list" << exit(FatalError);
    }
    return loadBalancedChemistryModel<ThermoType>::tc();
}



} // namespace Foam
