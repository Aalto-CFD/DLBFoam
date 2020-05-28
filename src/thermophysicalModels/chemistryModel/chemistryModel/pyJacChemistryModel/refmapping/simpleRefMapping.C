#include "simpleRefMapping.H"

namespace Foam{

int simpleRefMapping::compute_active_cells(PtrList<volScalarField>& Y) const{

    /* please someone implement, I have no idea what the other functions do */

    return 54;
}




bool simpleRefMapping::check_if_refcell(chemistryProblem& problem)
{
    //Note this assumes that mixture_fraction.update() has been called!
    auto beta_of = mixture_fraction_.get_beta();
    auto alpha = mixture_fraction_.get_alpha();

    scalar beta = 0.0; //TODO: rename!
    scalar Z;
    forAll(problem.c, iField)
    {
        const scalar& Yi = problem.c[iField];
        beta += alpha[iField]*Yi;
    }
    Z = (beta - beta_of[0])/(beta_of[1] - beta_of[0]);
    if (Z>tolerance_)
    {
        return false;
    }
    else
    {
        return true;
    }

}

bool simpleRefMapping::shouldMap(chemistryProblem& problem) {
    return check_if_refcell(problem);
}

} //namespace Foam