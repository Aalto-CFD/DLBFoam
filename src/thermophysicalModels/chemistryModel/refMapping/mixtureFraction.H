/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

Class
    Foam::mixtureFraction

Description
    Calculates and outputs the mixture fraction calculated between the
    fuel and the oxidizer.
    \f
        Z = (beta - beta_0)/(beta_1 - beta_0)
    \f]

    with
    \f[
        beta = 2 \frac{Z_C}{W_C} + 0.5 \frac{Z_H}{W_H} - 1 \frac{Z_O}{W_O}
    \f]

    and
    \f[
        Z_k = \sum\limits_{i=1}^N(\frac{a_{k,i} W_k}{W_i} Y_i)
    \f]

    where
    \vartable
        Z       | Mixture fraction
        beta    | Conserved scalar
        beta_0  | Conserved scalar in the oxidiser stream
        beta_1  | Conserved scalar in the fuel stream
        Z_k     | Element mass fractions (k = [C, H, O])
        W_k     | Atomic weights (k = [C, H, O])
        a_{k,i} | Number of atoms of element 'k' in specie 'i'
        N       | Number of species
        Y_i     | Species mass fractions
    \endvartable

References:

    \verbatim
        Bilger, R.; Stårner, S. & Kee, R.
        On reduced mechanisms for methane-air combustion in nonpremixed flames
        Combustion and Flame, 1990, 80, 135 - 149

        J. C. Sutherland, P. J. Smith, And J. H. Chen
        Quantification of differential diffusion in nonpremixed systems
        Combustion Theory and Modelling, Vol. 9, No. 2, May 2005, 365–383

        "Large Eddy Simulation of Fuel Spray Combustion"
        Armin Wehrfritz,
        PhD thesis, Aalto University, Finland, 2016
    \endverbatim

SourceFiles
    mixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureFraction_H
#define mixtureFraction_H

#include "IOdictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// Mixture fraction headers
#include "atomicWeights.H"
#include "psiReactionThermo.H"

namespace Foam
{

class mixtureFraction
{

public:
    mixtureFraction() = default;

    mixtureFraction(const dictionary& mixFracDict, const wordList& species);

    const List<scalar>& get_alpha() const
    {
        return alpha_;
    }
    const List<scalar>& get_beta() const
    {
        return beta_;
    }

    ///
    ///@brief Updates the alpha and beta parameters based on the input
    /// composition
    ///
    ///@param composition
    ///
    void update(const basicSpecieMixture& composition);

    ///
    ///@brief Prints information on related to mixture fraction.
    ///
    ///@param composition
    ///@param Yconst
    ///
    void print_information(
        const basicSpecieMixture& composition,
        const List<List<scalar>>& Yconst) const;

private:
    dictionary   mixFracDict_;
    wordList     species_;
    List<scalar> alpha_;
    List<scalar> beta_;

    ///
    ///@brief updates alpha used in mixture fraction calculation.
    ///
    ///@param composition
    ///
    void update_alpha(const basicSpecieMixture& composition);

    ///
    ///@brief updates beta used in mixture fraction calculation.
    ///
    ///@param Yconst
    ///@param alpha
    ///
    void
    update_beta(const List<List<scalar>>& Yconst, const List<scalar>& alpha);

    ///
    ///@brief Sets the pure fuel and pure oxidizer stream conditions.
    ///
    ///@param composition
    ///@return List<List<scalar>>
    ///
    List<List<scalar>>
    compute_yconst(const basicSpecieMixture& composition) const;

    ///
    ///@brief Updates a_{k,i} for every specie in the system.
    ///
    ///@param composition
    ///@return List<List<scalar>>
    ///
    List<List<scalar>> compute_a(const basicSpecieMixture& composition) const;

    ///
    ///@brief Computes the normalized mixture fraction at oxidizer port.
    ///
    ///@param a
    ///@param Yconst
    ///@return List<scalar>
    ///
    List<scalar> compute_Zox(
        const List<List<scalar>>& a, const List<List<scalar>>& Yconst) const;

    ///
    ///@brief Computes the normalized mixture fraction at fuel port.
    ///
    ///@param a
    ///@param Yconst
    ///@return List<scalar>
    ///
    List<scalar> compute_Zfu(
        const List<List<scalar>>& a, const List<List<scalar>>& Yconst) const;
};

} // namespace Foam

#endif