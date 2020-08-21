/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::loadBalancedChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    loadBalancedChemistryModelI.H
    loadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalancedChemistryModel_H
#define loadBalancedChemistryModel_H

#include "StandardChemistryModel.H"
#include "clockTime.H"
#include "OFstream.H"


#include "globalBalancingMethod.H"
#include "bulutLoadBalancing.H"
#include "chemistryProblem.H"
#include "chemistrySolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class loadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template <class ReactionThermo, class ThermoType>
class loadBalancedChemistryModel : public StandardChemistryModel<ReactionThermo, ThermoType> 
{
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);


    autoPtr<OFstream> cpuSolveFile_;

    inline autoPtr<OFstream> logFile(const word& name) const;

    //using balancer_ptr = typename std::unique_ptr<chemistryLoadBalancingMethod>;
    using balancer_ptr = autoPtr<chemistryLoadBalancingMethod>;


public:
    //- Runtime type information
    TypeName("loadBalanced");

    // Constructors

    //- Construct from thermo
    loadBalancedChemistryModel(const ReactionThermo& thermo);

    //- Destructor
    virtual ~loadBalancedChemistryModel();

    //- Disallow default bitwise copy construction
    loadBalancedChemistryModel(const loadBalancedChemistryModel&);


    virtual scalar solve(const scalar deltaT) override;
    virtual scalar solve(const scalarField& deltaT) override;

    virtual void solve(scalar&      p,
                       scalar&      T,
                       scalarField& c,
                       const label  li,
                       scalar&      deltaT,
                       scalar&      subDeltaT) const = 0;

    ///
    ///@brief 
    ///
    ///@param prob 
    ///@param soln 
    ///
    void solve_single(chemistryProblem& prob, chemistrySolution& soln) const;

    //- Disallow default bitwise assignment
        void operator=(const loadBalancedChemistryModel&) = delete;



private:

    scalarField                   cpu_times_;

    balancer_ptr load_balancer_;

    ///
    ///@brief A factory function for creating the balancer
    ///
    ///@return chemistryLoadBalancingMethod* 
    ///
    static balancer_ptr create_balancer();


    ///
    ///@brief 
    ///
    ///@param problems 
    ///@return chemistryLoadBalancingMethod::buffer_t<chemistrySolution> 
    ///
    chemistryLoadBalancingMethod::buffer_t<chemistrySolution>
    solve_buffer(chemistryLoadBalancingMethod::buffer_t<chemistryProblem>& problems) const;


    ///
    ///@brief Get the problems object
    ///
    ///@param Y_ 
    ///@param deltaT 
    ///@return DynamicList<chemistryProblem> 
    ///
    
    template<class DeltaTType>
    DynamicList<chemistryProblem> get_problems(const DeltaTType& deltaT);


    ///
    ///@brief 
    ///
    ///@param solutions 
    ///@return scalar 
    ///
    scalar update_reaction_rates(const  chemistryLoadBalancingMethod::buffer_t<chemistrySolution>& solutions);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "loadBalancedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
