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
#include "fvcGrad.H"

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
    thermo(lookupObject<fluidMulticomponentThermo>(physicalProperties::typeName)),
    Y(thermo.Y()),
    species_(thermo.species()),
    alpha_(thermo.species().size(), 0.0),
    beta_(2, 0.0),
    Z
    (
        IOobject
        (
            "Z",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    equivalenceRatio
    (
        IOobject
        (
            "equivalenceRatio",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    mixFracDict_(dict)
{
    initialize(thermo);
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mixtureFraction::~mixtureFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::mixtureFraction::fields() const
{
    return wordList{};
}

void Foam::functionObjects::mixtureFraction::initialize(const fluidMulticomponentThermo& thermo)
{
    forAll(alpha_,i)
    {
        const dictionary& dict =
            mixFracDict_.subDict(species_[i]).subDict("elements");
        scalar a0(
            2.0 * dict.lookupOrDefault<label>("C", 0) / thermo.WiValue(i));
        scalar a1(
            0.5 * dict.lookupOrDefault<label>("H", 0) / thermo.WiValue(i));
        scalar a2(
            -1.0 * dict.lookupOrDefault<label>("O", 0) / thermo.WiValue(i));
        alpha_[i] = a0 + a1 + a2;

    }
    List<List<scalar>> Yconst(2, List<scalar>(species_.size(), 0.0));
    forAll(species_, i)
    {
        Yconst[0][i] = mixFracDict_.subDict("oxidizerMassFractions")
                           .lookupOrDefault<scalar>(species_[i], 0.0);
        Yconst[1][i] = mixFracDict_.subDict("fuelMassFractions")
                           .lookupOrDefault<scalar>(species_[i], 0.0);
    }
    scalar YoxTot = 0.0;
    scalar YfuTot = 0.0;
    forAll(species_, i)
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

    forAll(species_, i)
    {
        beta_[0] += alpha_[i] * Yconst[0][i]; // oxidizer
        beta_[1] += alpha_[i] * Yconst[1][i]; // fuel
    }
    // Stoichiometric mixture fraction
    Z_st = (0.0 - beta_[0]) / (beta_[1] - beta_[0]);
    AFR_st = beta_[1] / (-beta_[0]);

    Info << "Stoichiometric mixture fraction is: " << Z_st << endl;
    Info << "Stoichiometric air/fuel ratio is: " << AFR_st << endl;
    if (mag(1.0 - max(Z_st / (1.0 - Z_st) * AFR_st, 0.0)) > SMALL)
    {
        FatalErrorIn("mixtureFraction.C :")
            << "The stoichiometric mixture fraction does not convert to unity equivalence ratio"
            << abort(FatalError);
    }
}

bool Foam::functionObjects::mixtureFraction::execute()
{
    forAll(Z,celli)
    {

        scalar beta = 0.0;
        forAll(Y, iField)
        {
            const scalarField& Yi = Y[iField];
            beta += alpha_[iField] * Yi[celli];
        }
        Z[celli] = min(max((beta - beta_[0]) / (beta_[1] - beta_[0]),0.0),1.0);

        equivalenceRatio[celli] = max(Z[celli] / ( (1.0 - Z[celli]) + SMALL) * AFR_st, 0.0);

    }
    return true;
}


bool Foam::functionObjects::mixtureFraction::write()
{
    Z.correctBoundaryConditions();
    Z.write();

    equivalenceRatio.correctBoundaryConditions();
    equivalenceRatio.write();
    return true;
}


// ************************************************************************* //