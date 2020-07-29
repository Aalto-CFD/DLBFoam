#include "simpleRefMappingMethod.H"

namespace Foam {


bool simpleRefMappingMethod::check_if_refcell(const scalarField& mass_fraction) const{
    // Note this assumes that mixture_fraction.update() has been called!
    auto beta_of = mixture_fraction_.get_beta();
    auto alpha   = mixture_fraction_.get_alpha();

    scalar beta = 0.0; // TODO: rename!
    forAll(mass_fraction, iField) { beta += alpha[iField] * mass_fraction[iField]; }

    scalar Z = (beta - beta_of[0]) / (beta_of[1] - beta_of[0]);
    if (Z > tolerance_) {
        return false;
    } else {
        return true;
    }
}

bool simpleRefMappingMethod::shouldMap(const scalarField& mass_fraction) const{ return check_if_refcell(mass_fraction); }

} // namespace Foam