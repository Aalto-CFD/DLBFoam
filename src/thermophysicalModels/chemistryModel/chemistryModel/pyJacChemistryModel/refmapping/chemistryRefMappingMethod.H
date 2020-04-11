/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::chemistryRefMappingMethod

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryRefMappingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryRefMappingMethod_H
#define chemistryRefMappingMethod_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"



// Mixture fraction headers
#include "atomicWeights.H"
#include "psiReactionThermo.H"
//#include "BasicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class chemistryRefMappingMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryRefMappingMethod
{

public:


    

    

    // Constructors

    //- Construct from components
    chemistryRefMappingMethod
    (
        const dictionary& dict
    )
    :
    dict_(dict),
    coeffsDict_(dict.subDict("refmapping")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false))
    {


        if (this->active_)
        {
        // scalar test = readScalar(mixFracDict.lookup("test"));
        }
    }


    ///
    ///@brief THIS IS THE BREAD AND BUTTER WHICH SHOULD BE DEFINED IN ALL CLASSES!
    ///       THIS IS NOT PURE AT THIS POINT 
    ///
    // TODO add params
    virtual void applyMapping() const = 0;


    


    void init_mixture_fraction(basicSpecieMixture& composition){
      
        dictionary mixFracDict(this->coeffsDict_.subDict("mixtureFractionProperties"));

        const wordList& species = composition.species();

        alpha = compute_alpha(composition);
        List<List<scalar>> Yconst = compute_yconst(composition);
        beta_of = compute_beta_of(species, Yconst, alpha);
        

        // Compute element mass fraction coefficients (only for printing information)
        List<List<scalar> > a = compute_a(composition);
        

        // Compute element mass fractions (only for printing)
        List<scalar> Zox = compute_Zox(a, Yconst, species);
        List<scalar> Zfu = compute_Zfu(a, Yconst, species);
         
        print_information(Zox, Zfu, Yconst, beta_of, species);

    };


    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli)
    {
        scalar beta = 0.0;
        scalar Z;
        forAll(Y, iField)
        {
            const scalarField& Yi = Y[iField];
            beta += alpha[iField]*Yi[celli];
        }
        Z = (beta - beta_of[0])/(beta_of[1] - beta_of[0]);
        if (Z>1e-5)
        {
            return false;
        }
        else
        {
            return true;
        }

    }


    //- Destructor
    virtual ~chemistryRefMappingMethod() {}


    // Member Functions

    inline bool active()
    {
        return active_;
    }
    


protected:

    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;


    int mytest;
    List<scalar> alpha; //what is this?
    List<scalar> beta_of; //what is this?


private:

    ///
    ///@brief Computes alpha. What is alpha?
    ///
    ///@param composition 
    ///@return List<scalar> 
    ///
    List<scalar> compute_alpha(basicSpecieMixture& composition) const{

        dictionary mixFracDict(this->coeffsDict_.subDict("mixtureFractionProperties"));
        const wordList& species = composition.species();
        List<scalar> ret(species.size(), 0.0);
        forAll(ret, i)
        {
            const dictionary& dict = mixFracDict.subDict(species[i]).subDict("elements");
            scalar a0
            (
                2.0
            *dict.lookupOrDefault<label>("C",0)
            /composition.W(i)
            );
            scalar a1
            (
                0.5
            *dict.lookupOrDefault<label>("H",0)
            /composition.W(i)
            );
            scalar a2
            (
            -1.0
            *dict.lookupOrDefault<label>("O",0)
            /composition.W(i)
            );
            ret[i] = a0 + a1 + a2;
        }
        return ret;

    }

    ///
    ///@brief Computes yconst. What is yconst?
    ///
    ///@param composition 
    ///@return List<List<scalar>> 
    ///
    List<List<scalar>> compute_yconst(basicSpecieMixture& composition) const {

        dictionary mixFracDict(this->coeffsDict_.subDict("mixtureFractionProperties"));
        const wordList& species = composition.species();
        List<List<scalar> > Yconst( 2, List<scalar>(species.size(),0.0) );
        forAll(species, i)
        {
            Yconst[0][i] =
                mixFracDict.subDict
                (
                    "oxidizerMassFractions"
                ).lookupOrDefault<scalar>(species[i], 0.0);
            Yconst[1][i] =
                mixFracDict.subDict
                (
                    "fuelMassFractions"
                ).lookupOrDefault<scalar>(species[i], 0.0);
        }
        return Yconst;

    }


    List<scalar> compute_beta_of(const wordList& species, const List<List<scalar>>& Yconst, const List<scalar>& alphaparam) const {


        scalar YoxTot = 0.0;
        scalar YfuTot = 0.0;
        forAll(species, i)
        {
            YoxTot += Yconst[0][i];
            YfuTot += Yconst[1][i];
        }

        if (mag(1.0 - YoxTot) > SMALL || mag(1.0 - YfuTot) > SMALL)
        {
            FatalErrorIn
            (
                "createMixtureFraction.H :"
            )
            << "oxidizerMassFractions or fuelMassFractions do not sum up to 1.0"
            << abort(FatalError);
        }

        // constructing the beta function at the oxidizer and fuel streams
        List<scalar> beta_of_(2, 0.0);
        forAll(species, i)
        {
            beta_of_[0] += alphaparam[i]*Yconst[0][i]; // oxidizer
            beta_of_[1] += alphaparam[i]*Yconst[1][i]; // fuel
        }

        return beta_of_;

    }

    ///
    ///@brief Computes a
    ///
    ///@param composition 
    ///@return List<List<scalar>> 
    ///
    List<List<scalar>> compute_a(basicSpecieMixture& composition) const{

        dictionary mixFracDict(this->coeffsDict_.subDict("mixtureFractionProperties"));
        const wordList& species = composition.species();

        List<List<scalar> > a(species.size(), List<scalar>(3,0.0));
        forAll(species, i)
        {
            const dictionary& dict = mixFracDict.subDict(species[i]).subDict("elements");
            a[i][0] =
                dict.lookupOrDefault<label>("C",0)
            *atomicWeights["C"]
            /composition.W(i) ;
            a[i][1] =
                dict.lookupOrDefault<label>("H",0)
            *atomicWeights["H"]
            /composition.W(i) ;
            a[i][2] =
                dict.lookupOrDefault<label>("O",0)
            *atomicWeights["O"]
            /composition.W(i) ;
        }
        return a;
    }
    
    ///
    ///@brief Computes Zox
    ///
    ///
    List<scalar> compute_Zox(const List<List<scalar>>& a, const List<List<scalar>>& Yconst, const wordList& species ) const{
    
        // Compute element mass fractions (only for printing)
        List<scalar> Zox(species.size(),0.0);
        
        forAll(Zox, j)
        {
            forAll(species, i)
            {
                Zox[j] += a[i][j]*Yconst[0][i];
            }
        }
        return Zox;
    }

    List<scalar> compute_Zfu(const List<List<scalar>>& a, const List<List<scalar>>& Yconst, const wordList& species ) const{

        List<scalar> Zfu(species.size(),0.0);
         
        forAll(Zfu, j)
        {
            forAll(species, i)
            {
                Zfu[j] += a[i][j]*Yconst[1][i];
            }
        }
        return Zfu;
    }


    void print_information(const List<scalar>& Zox, 
                           const List<scalar>& Zfu, 
                           const List<List<scalar>>& Yconst, 
                           const List<scalar>& beta_of_,
                           const wordList& species ) const{
                           

        //- stoichiometric mixture fraction
        scalar Z_st =  (0.0 - beta_of_[0])/(beta_of_[1] - beta_of_[0]);

        // Print out some information
        Info<< "    Oxidizer:" << nl
            << "    speciesMassFractions:" << nl;
        forAll(species, i)
        {
            if (Yconst[0][i] != 0.0)
            {
                Info<< "        Y("<<species[i]<<") = " << Yconst[0][i] << nl;
            }
        }
        Info<< "    elementMassFractions:" << nl
            << "        Z(C) = " << Zox[0] << nl
            << "        Z(H) = " << Zox[1] << nl
            << "        Z(O) = " << Zox[2] << nl
            << endl;

        Info<< "    Fuel:" << nl
            << "    speciesMassFractions:" << nl;
        forAll(species, i)
        {
            if (Yconst[1][i] != 0.0)
            {
                Info<< "        Y("<<species[i]<<") = " << Yconst[1][i] << nl;
            }
        }
        Info<< "    elementMassFractions:" << nl
            << "        Z(C) = " << Zfu[0] << nl
            << "        Z(H) = " << Zfu[1] << nl
            << "        Z(O) = " << Zfu[2] << nl
            << endl;

        // Stoichiometric mixture fraction
        Info<<"    Stoichiomentric mixture fraction Zst = " << Z_st << endl;


    }    

        




};



} //namespace Foam


#endif

// ************************************************************************* //
