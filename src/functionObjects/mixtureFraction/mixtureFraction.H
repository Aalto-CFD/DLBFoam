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

Class
    Foam::functionObjects::mixtureFraction

Description
    Calculates and outputs the mixture fraction calculated between the
    fuel and the oxidizer. It also outputs the equivalence ratio, calculated
    from the mixture fraction.

    \f
        Z = (beta - beta_0)/(beta_1 - beta_0)
    \f]
    and
    \f
        phi = max(Z / (1.0 - Z) * AFR_st, 0.0);
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

    Example of function object specification:
    \verbatim
    mixtureFraction
    {
        type        mixtureFraction;
        libs        ("libFunctionObjects_DLB.so");

        oxidizerMassFractions
        {
            N2              0.77;
            O2              0.23;
        }
        fuelMassFractions
        {
            CH4           1;
        }
    }
    \endverbatim

Usage
    \table
        Property  | Description                 | Required  | Default value
        type      | type name: mixtureFraction   | yes       |
    \endtable

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
        PhD thesis, 2020, Aalto University, Finland, 2016
    \endverbatim

SourceFiles
    mixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureFraction_functionObject_H
#define mixtureFraction_functionObject_H

#include "fvMeshFunctionObject.H"
#include "writeLocalObjects.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "volFieldsFwd.H"
#include "psiMulticomponentThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class mixtureFraction Declaration
\*---------------------------------------------------------------------------*/

class mixtureFraction
:
    public fvMeshFunctionObject,
    public writeLocalObjects
{
    // Private Data

    List<scalar> alpha_;
    List<scalar> beta_;
    scalar Z_st_;
    scalar AFR_st_;


    // Private Member Functions

        //- Initialize the mixture fraction coefficients
        void initialise
        (
            const dictionary& dict,
            const fluidMulticomponentThermo& thermo
        );

public:

    //- Runtime type information
    TypeName("mixtureFraction");


    // Constructors

        //- Construct from Time and dictionary
        mixtureFraction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        mixtureFraction(const mixtureFraction&) = delete;

    //- Destructor
    virtual ~mixtureFraction();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate the mixture fraction
        virtual bool execute();

        //- Write the mixture fraction fields
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mixtureFraction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //