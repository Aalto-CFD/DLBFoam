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

Class
    Foam::mixtureFractionRefMapper

Description
    Extends the base class RefMapperBase by implementing a
    reference mapping method using Bilger's mixture fraction definition,
    calculated between fuel and oxidizer.

References
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
    mixtureFractionRefMapper.C

See also
    mixtureFraction.H
\*---------------------------------------------------------------------------*/

#ifndef mixtureFractionRefMapper_H
#define mixtureFractionRefMapper_H

#include "IOdictionary.H"
#include "RefMapperBase.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class mixtureFractionRefMapper : public RefMapperBase
{

public:
    mixtureFractionRefMapper() = default;

    mixtureFractionRefMapper(
        const dictionary& dict, const basicSpecieMixture& composition)
        : RefMapperBase(), dict_(dict), coeffsDict_(dict.subDict("refmapping")),
          active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
          tolerance_(coeffsDict_.lookupOrDefault<scalar>("tolerance", 1e-4)),
          temperatureDict_(coeffsDict_.subDict("temperatureChecker")),
          temp_active_(
              temperatureDict_.lookupOrDefault<Switch>("active", false)),
          temperature_tolerance_(
              temperatureDict_.lookupOrDefault<scalar>("tolerance", 0.1)),
          refCellFound_(false), Tref_(0.0),
          mixture_fraction_(
              coeffsDict_.subDict("mixtureFractionProperties"),
              composition.species())
    {
        // TODO: move to the mixture fraction constructor and dont call anything
        // here
        mixture_fraction_.update(composition);
    }

    bool shouldMap(const ChemistryProblem& problem);

    bool active() const
    {
        return active_;
    }

protected:
    // None

private:
    // TODO: for runtime switching, these need to be references
    const dictionary dict_;
    const dictionary coeffsDict_;

    //- Is reference mapping active?
    Switch active_;

    //- Tolerance for reference mapping
    scalar           tolerance_;

    const dictionary temperatureDict_;

    //- Is temperature checking active?
    Switch temp_active_;

    //- Tolerance for temperature checking (in Kelvins)
    scalar temperature_tolerance_;

    bool   refCellFound_;
    scalar Tref_;

    mixtureFraction mixture_fraction_;

    bool check_if_refcell(const ChemistryProblem& mass_fraction);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
