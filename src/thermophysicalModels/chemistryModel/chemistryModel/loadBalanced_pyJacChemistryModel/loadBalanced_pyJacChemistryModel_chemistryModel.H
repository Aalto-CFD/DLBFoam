#ifndef loadBalanced_pyJacChemistryModel_chemistryModel_H
#define loadBalanced_pyJacChemistryModel_chemistryModel_H

#include "loadBalanced_pyJacChemistryModel.H"

namespace Foam
{
namespace chemistryModels
{
template<class ThermoType>
using loadBalanced_pyJacChemistryModel = Foam::loadBalanced_pyJacChemistryModel<ThermoType>;
}
}

// Keep Standard reaction/reduction registrations enabled for this derived model.
#define loadBalanced_pyJacChemistryModel_chemistryModel 1
#define DLB_CHEM_MODEL_LOADBALANCED_PYJAC 1

#endif
