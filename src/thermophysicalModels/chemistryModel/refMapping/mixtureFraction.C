/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

namespace Foam
{

// Constructor
mixtureFraction::mixtureFraction(
    const dictionary& mixFracDict, const wordList& species)
    : mixFracDict_(mixFracDict), species_(species), alpha_(species.size(), 0.0),
      beta_(2, 0.0)
{
}

void mixtureFraction::update(const basicSpecieMixture& composition)
{

    update_alpha(composition);
    auto Yconst = compute_yconst(composition);
    update_beta(Yconst, alpha_);
    print_information(composition, Yconst);
}

/// TODO: simplify this and consider if necessary
void mixtureFraction::print_information(
    const basicSpecieMixture& composition,
    const List<List<scalar>>& Yconst) const
{

    auto a   = compute_a(composition);
    auto Zox = compute_Zox(a, Yconst);
    auto Zfu = compute_Zfu(a, Yconst);

    //- stoichiometric mixture fraction
    scalar Z_st = (0.0 - beta_[0]) / (beta_[1] - beta_[0]);

    // Print out some information
    Info << "    Oxidizer:" << nl << "    speciesMassFractions:" << nl;
    forAll(species_, i)
    {
        if(Yconst[0][i] != 0.0)
        {
            Info << "        Y(" << species_[i] << ") = " << Yconst[0][i] << nl;
        }
    }
    Info << "    elementMassFractions:" << nl << "        Z(C) = " << Zox[0]
         << nl << "        Z(H) = " << Zox[1] << nl
         << "        Z(O) = " << Zox[2] << nl << endl;

    Info << "    Fuel:" << nl << "    speciesMassFractions:" << nl;
    forAll(species_, i)
    {
        if(Yconst[1][i] != 0.0)
        {
            Info << "        Y(" << species_[i] << ") = " << Yconst[1][i] << nl;
        }
    }
    Info << "    elementMassFractions:" << nl << "        Z(C) = " << Zfu[0]
         << nl << "        Z(H) = " << Zfu[1] << nl
         << "        Z(O) = " << Zfu[2] << nl << endl;

    // Stoichiometric mixture fraction
    Info << "    Stoichiomentric mixture fraction Zst = " << Z_st << endl;
}

void mixtureFraction::update_alpha(const basicSpecieMixture& composition)
{

    forAll(alpha_, i)
    {
        const dictionary& dict =
            mixFracDict_.subDict(species_[i]).subDict("elements");
        scalar a0(
            2.0 * dict.lookupOrDefault<label>("C", 0) / composition.Wi(i));
        scalar a1(
            0.5 * dict.lookupOrDefault<label>("H", 0) / composition.Wi(i));
        scalar a2(
            -1.0 * dict.lookupOrDefault<label>("O", 0) / composition.Wi(i));
        alpha_[i] = a0 + a1 + a2;
    }
}

void mixtureFraction::update_beta(
    const List<List<scalar>>& Yconst, const List<scalar>& alpha)
{

    scalar YoxTot = 0.0;
    scalar YfuTot = 0.0;
    forAll(species_, i)
    {
        YoxTot += Yconst[0][i];
        YfuTot += Yconst[1][i];
    }

    if(mag(1.0 - YoxTot) > SMALL || mag(1.0 - YfuTot) > SMALL)
    {
        FatalErrorIn("createMixtureFraction.H :")
            << "oxidizerMassFractions or fuelMassFractions do not sum up to 1.0"
            << abort(FatalError);
    }

    forAll(species_, i)
    {
        beta_[0] += alpha[i] * Yconst[0][i]; // oxidizer
        beta_[1] += alpha[i] * Yconst[1][i]; // fuel
    }
}

List<List<scalar>>
mixtureFraction::compute_yconst(const basicSpecieMixture& composition) const
{

    List<List<scalar>> Yconst(2, List<scalar>(species_.size(), 0.0));
    forAll(species_, i)
    {
        Yconst[0][i] = mixFracDict_.subDict("oxidizerMassFractions")
                           .lookupOrDefault<scalar>(species_[i], 0.0);
        Yconst[1][i] = mixFracDict_.subDict("fuelMassFractions")
                           .lookupOrDefault<scalar>(species_[i], 0.0);
    }
    return Yconst;
}

List<List<scalar>>
mixtureFraction::compute_a(const basicSpecieMixture& composition) const
{

    List<List<scalar>> a(species_.size(), List<scalar>(3, 0.0));
    forAll(species_, i)
    {
        const dictionary& dict =
            mixFracDict_.subDict(species_[i]).subDict("elements");
        a[i][0] = dict.lookupOrDefault<label>("C", 0) * atomicWeights["C"] /
                  composition.Wi(i);
        a[i][1] = dict.lookupOrDefault<label>("H", 0) * atomicWeights["H"] /
                  composition.Wi(i);
        a[i][2] = dict.lookupOrDefault<label>("O", 0) * atomicWeights["O"] /
                  composition.Wi(i);
    }
    return a;
}

List<scalar> mixtureFraction::compute_Zox(
    const List<List<scalar>>& a, const List<List<scalar>>& Yconst) const
{

    // Compute element mass fractions (only for printing)
    List<scalar> Zox(species_.size(), 0.0);

    forAll(Zox, j)
    {
        forAll(species_, i)
        {
            Zox[j] += a[i][j] * Yconst[0][i];
        }
    }
    return Zox;
}

List<scalar> mixtureFraction::compute_Zfu(
    const List<List<scalar>>& a, const List<List<scalar>>& Yconst) const
{

    List<scalar> Zfu(species_.size(), 0.0);

    forAll(Zfu, j)
    {
        forAll(species_, i)
        {
            Zfu[j] += a[i][j] * Yconst[1][i];
        }
    }
    return Zfu;
}

} // namespace Foam