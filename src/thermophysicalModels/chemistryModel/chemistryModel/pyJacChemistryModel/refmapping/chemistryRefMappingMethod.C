#include "chemistryRefMappingMethod.H"


namespace Foam{


//- Construct from dict and a list of species
chemistryRefMappingMethod::chemistryRefMappingMethod(const dictionary& dict, const wordList& species)
:
dict_(dict),
coeffsDict_(dict.subDict("refmapping")),
active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
mixture_fraction_(coeffsDict_.subDict("mixtureFractionProperties"), species)
{
}

}//namespace Foam

