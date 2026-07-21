#ifndef loadBalancedChemistryModel_chemistryModel_H
#define loadBalancedChemistryModel_chemistryModel_H

#include "loadBalancedChemistryModel.H"

namespace Foam
{
namespace chemistryModels
{
template<class ThermoType>
using loadBalancedChemistryModel = Foam::loadBalancedChemistryModel<ThermoType>;
}
}

// Keep Standard reaction/reduction registrations enabled for this derived model.
#define loadBalancedChemistryModel_chemistryModel 1
#define DLB_CHEM_MODEL_LOADBALANCED 1

#endif
