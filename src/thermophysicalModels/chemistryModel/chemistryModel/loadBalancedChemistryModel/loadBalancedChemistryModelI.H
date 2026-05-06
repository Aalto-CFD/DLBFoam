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

#include "OSspecific.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::loadBalancedChemistryModel<ThermoType>::logFile(
    const word& name) const
{
    mkDir(this->mesh().time().path() / "loadBal" / this->group());
    return autoPtr<OFstream>(new OFstream(
        this->mesh().time().path() / "loadBal" / this->group() / name));
}

template <class ThermoType>
inline const Foam::volScalarField::Internal&
Foam::loadBalancedChemistryModel<ThermoType>::RR(const label i) const{
    const PtrList<volScalarField::Internal>& rr = chemistryModel<ThermoType>::RR();
    return rr[i];
}

// Non const access to chemical source terms was removed in OF 11
// and thus casting is required
template <class ThermoType>
inline Foam::volScalarField::Internal&
Foam::loadBalancedChemistryModel<ThermoType>::RR(const label i) {
    PtrList<volScalarField::Internal>& rr =
    	const_cast<PtrList<volScalarField::Internal>&>(chemistryModel<ThermoType>::RR());
    return rr[i];
}

template <class ThermoType>
inline Foam::scalar Foam::loadBalancedChemistryModel<ThermoType>::timeValue() const
{
    if (repeat_ > 0)
    {
        return begin_ + fmod(this->time().userTimeValue() - begin_, repeat_);
    }
    else
    {
        return this->time().userTimeValue();
    }
}

template <class ThermoType>
inline Foam::Switch Foam::loadBalancedChemistryModel<ThermoType>::chemistry() const
{
    Foam::scalar currentTime = timeValue();

    // Check if endTime_ is over the repeat period (wrap-around case)
    if (startTime_ < endTime_)
    {
        // Normal case: check if currentTime is within startTime_ and endTime_
        if (currentTime >= startTime_ && currentTime <= endTime_)
        {
            return this->chemistry_;
        }
        else
        {
            return false;
        }
    }
    else
    {
        // Wrap-around case: active period spans from startTime_ to end of repeat, then from 0 to endTime_
        if (currentTime >= startTime_ || currentTime <= endTime_)
        {
            return this->chemistry_;
        }
        else
        {
            return false;
        }
    }

    return false;
}

// ************************************************************************* //
