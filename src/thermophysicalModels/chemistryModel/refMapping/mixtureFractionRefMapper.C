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

#include "mixtureFractionRefMapper.H"

namespace Foam
{

bool mixtureFractionRefMapper::check_if_refcell(
    const scalarField& mass_fraction) const
{
    // Note this assumes that mixture_fraction.update() has been called!
    auto beta_of = mixture_fraction_.get_beta();
    auto alpha   = mixture_fraction_.get_alpha();

    scalar beta = 0.0; // TODO: rename!
    forAll(mass_fraction, iField)
    {
        beta += alpha[iField] * mass_fraction[iField];
    }

    scalar Z = (beta - beta_of[0]) / (beta_of[1] - beta_of[0]);
    if(Z > tolerance_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool mixtureFractionRefMapper::shouldMap(const scalarField& mass_fraction) const
{
    return check_if_refcell(mass_fraction);
}

} // namespace Foam