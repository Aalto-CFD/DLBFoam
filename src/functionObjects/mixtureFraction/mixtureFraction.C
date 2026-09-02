/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "mixtureFraction.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mixtureFraction, 0);
    addToRunTimeSelectionTable(functionObject, mixtureFraction, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mixtureFraction::mixtureFraction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_),
    alpha_(),
    beta_(2, 0.0),
    Z_st_(0.0),
    AFR_st_(0.0)
{
    read(dict);
    resetLocalObjectNames({"Z", "equivalenceRatio"});
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mixtureFraction::~mixtureFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mixtureFraction::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const fluidMulticomponentThermo& thermo =
        lookupObject<fluidMulticomponentThermo>(physicalProperties::typeName);

    initialise(dict, thermo);

    return true;
}


void Foam::functionObjects::mixtureFraction::initialise
(
    const dictionary& dict,
    const fluidMulticomponentThermo& thermo
)
{
    const speciesTable& species = thermo.species();
    const physicalProperties& thermoProperties =
        lookupObject<physicalProperties>(physicalProperties::typeName);

    alpha_.setSize(species.size(), 0.0);
    beta_ = 0.0;

    forAll(alpha_, i)
    {
        const dictionary& elements =
            thermoProperties.subDict(species[i]).subDict("elements");
        scalar a0(
            2.0 * elements.lookupOrDefault<label>("C", 0) / thermo.WiValue(i));
        scalar a1(
            0.5 * elements.lookupOrDefault<label>("H", 0) / thermo.WiValue(i));
        scalar a2(
            -elements.lookupOrDefault<label>("O", 0) / thermo.WiValue(i));
        alpha_[i] = a0 + a1 + a2;
    }

    List<List<scalar>> Yconst(2, List<scalar>(species.size(), 0.0));
    forAll(species, i)
    {
        Yconst[0][i] = dict.subDict("oxidizerMassFractions")
                           .lookupOrDefault<scalar>(species[i], 0.0);
        Yconst[1][i] = dict.subDict("fuelMassFractions")
                           .lookupOrDefault<scalar>(species[i], 0.0);
    }

    scalar YoxTot = 0.0;
    scalar YfuTot = 0.0;
    forAll(species, i)
    {
        YoxTot += Yconst[0][i];
        YfuTot += Yconst[1][i];
    }

    if (mag(1.0 - YoxTot) > SMALL || mag(1.0 - YfuTot) > SMALL)
    {
        FatalErrorIn("mixtureFraction.C :")
            << "oxidizerMassFractions or fuelMassFractions do not sum up to 1.0"
            << abort(FatalError);
    }

    forAll(species, i)
    {
        beta_[0] += alpha_[i] * Yconst[0][i]; // oxidizer
        beta_[1] += alpha_[i] * Yconst[1][i]; // fuel
    }

    // Stoichiometric mixture fraction
    Z_st_ = -beta_[0] / (beta_[1] - beta_[0]);
    AFR_st_ = beta_[1] / (-beta_[0]);

    Info << "Stoichiometric mixture fraction is: " << Z_st_ << endl;
    Info << "Stoichiometric air/fuel ratio is: " << AFR_st_ << endl;
    if (mag(1.0 - max(Z_st_ / (1.0 - Z_st_) * AFR_st_, 0.0)) > SMALL)
    {
        FatalErrorIn("mixtureFraction.C :")
            << "The stoichiometric mixture fraction does not convert to unity equivalence ratio"
            << abort(FatalError);
    }
}


Foam::wordList Foam::functionObjects::mixtureFraction::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::mixtureFraction::execute()
{
    const fluidMulticomponentThermo& thermo =
        lookupObject<fluidMulticomponentThermo>(physicalProperties::typeName);
    const PtrList<volScalarField>& Y = thermo.Y();

    tmp<volScalarField> tZ
    (
        volScalarField::New("Z", mesh_, dimensionedScalar(dimless, 0.0))
    );
    tmp<volScalarField> tEquivalenceRatio
    (
        volScalarField::New
        (
            "equivalenceRatio",
            mesh_,
            dimensionedScalar(dimless, 0.0)
        )
    );

    volScalarField& Z = tZ.ref();
    volScalarField& equivalenceRatio = tEquivalenceRatio.ref();

    forAll(Z, celli)
    {
        scalar beta = 0.0;
        forAll(Y, iField)
        {
            const scalarField& Yi = Y[iField];
            beta += alpha_[iField] * Yi[celli];
        }
        Z[celli] = min(max((beta - beta_[0]) / (beta_[1] - beta_[0]),0.0),1.0);

        equivalenceRatio[celli] =
            max(Z[celli] / ((1.0 - Z[celli]) + SMALL) * AFR_st_, 0.0);
    }

    Z.correctBoundaryConditions();
    equivalenceRatio.correctBoundaryConditions();

    store("Z", tZ);
    store("equivalenceRatio", tEquivalenceRatio);

    return true;
}


bool Foam::functionObjects::mixtureFraction::write()
{
    return writeLocalObjects::write();
}


// ************************************************************************* //