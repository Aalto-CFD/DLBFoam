#ifndef mixtureFraction_H
#define mixtureFraction_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"



// Mixture fraction headers
#include "atomicWeights.H"
#include "psiReactionThermo.H"

namespace Foam{


class mixtureFraction{

public:

    mixtureFraction() = default;

    mixtureFraction(const dictionary& mixFracDict, const wordList& species);


    const List<scalar>& get_alpha() const {return alpha_;}
    const List<scalar>& get_beta()  const {return beta_;}


    ///
    ///@brief Updates the alpha and beta parameters based on the input composition
    ///
    ///@param composition NO IDEA
    ///
    void update(basicSpecieMixture& composition);



    ///
    ///@brief Prints information on related to mixture fraction.
    ///
    ///@param composition 
    ///@param Yconst 
    ///
    void print_information(basicSpecieMixture& composition, const List<List<scalar>>& Yconst) const;
    


protected:

private:

    dictionary mixFracDict_;
    wordList species_;
    List<scalar> alpha_;
    List<scalar> beta_;

    ///
    ///@brief NO IDEA WHAT THIS DOES
    ///
    ///@param composition 
    ///
    void update_alpha(basicSpecieMixture& composition);


    ///
    ///@brief NO IDEA WHAT THIS DOES
    ///
    ///@param Yconst
    ///@param alpha 
    ///
    void update_beta(const List<List<scalar>>& Yconst, const List<scalar>& alpha);

    ///
    ///@brief NO IDEA WHAT THIS DOES
    ///
    ///@param composition 
    ///@return List<List<scalar>> 
    ///
    List<List<scalar>> compute_yconst(basicSpecieMixture& composition) const;

    ///
    ///@brief NO IDEA WHAT THIS DOES
    ///
    ///@param composition 
    ///@return List<List<scalar>> 
    ///
    //TODO: RENAME!
    List<List<scalar>> compute_a(basicSpecieMixture& composition) const;
    
    ///
    ///@brief NO IDEA WHAT THIS DOES
    ///
    ///@param a 
    ///@param Yconst 
    ///@return List<scalar> 
    ///
    List<scalar> compute_Zox(const List<List<scalar>>& a, const List<List<scalar>>& Yconst ) const;

    ///
    ///@brief NO IDEA WHAT THIS DOES
    ///
    ///@param a 
    ///@param Yconst 
    ///@return List<scalar> 
    ///
    List<scalar> compute_Zfu(const List<List<scalar>>& a, const List<List<scalar>>& Yconst ) const;
    

};


} //namespace Foam


#endif