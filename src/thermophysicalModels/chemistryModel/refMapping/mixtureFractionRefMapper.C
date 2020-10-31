/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM-Aalto library, derived from OpenFOAM.

    https://github.com/blttkgl/OpenFOAM-Aalto

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

bool Foam::mixtureFractionRefMapper::shouldMap(const ChemistryProblem& problem)
{
    scalar Z = mixture_fraction_.getZ(problem);

    if(!refCellFound_)
    {
        if(Z > tolerance_)
        {
            return false;
        }
        else
        {
            Tref_         = problem.Ti;
            refCellFound_ = true;
            return true;
        }
    }
    else
    {
        if((Z > tolerance_) || (abs(problem.Ti - Tref_) > deltaT_))
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

bool Foam::mixtureFractionRefMapper::shouldMap(const scalarField& concentration, scalar T)
{


    scalar Z = mixture_fraction_.getZ(concentration);

    if(!refCellFound_)
    {
        if(Z > tolerance_)
        {
            return false;
        }
        else
        {
            Tref_         = T; 
            refCellFound_ = true;
            return true;
        }
    }
    else
    {
        if((Z > tolerance_) || (abs(T - Tref_) > deltaT_))
        {
            return false;
        }
        else
        {
            return true;
        }
    }


}