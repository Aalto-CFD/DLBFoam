#include "../third_party/catch.hpp"
#include "mixtureFraction.H"


TEST_CASE("mixtureFraction constructors"){


    using namespace Foam;

    REQUIRE_NOTHROW(mixtureFraction());
    
    dictionary mix_frac_dict;
    mix_frac_dict.add("asd", 3);
    wordList species;

    REQUIRE_NOTHROW(mixtureFraction(mix_frac_dict, species));



    


    //const dictionary& mixFracDict, const wordList& species

}