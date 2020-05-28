#include "simpleRefMapping.H"

namespace Foam {


bool simpleRefMapping::check_if_refcell(const chemistryProblem& problem) const{
    // Note this assumes that mixture_fraction.update() has been called!
    auto beta_of = mixture_fraction_.get_beta();
    auto alpha   = mixture_fraction_.get_alpha();

    scalar beta = 0.0; // TODO: rename!
    forAll(problem.c, iField) { beta += alpha[iField] * problem.c[iField]; }

    scalar Z = (beta - beta_of[0]) / (beta_of[1] - beta_of[0]);
    if (Z > tolerance_) {
        return false;
    } else {
        return true;
    }
}

bool simpleRefMapping::shouldMap(const chemistryProblem& problem) const{ return check_if_refcell(problem); }

} // namespace Foam