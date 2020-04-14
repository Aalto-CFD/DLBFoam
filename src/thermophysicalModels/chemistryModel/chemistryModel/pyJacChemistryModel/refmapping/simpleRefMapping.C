#include "simpleRefMapping.H"

namespace Foam{

void simpleRefMapping::init_mixture_fraction(basicSpecieMixture& composition){
    mixture_fraction_.update(composition);
};

bool simpleRefMapping::check_if_refcell(PtrList<volScalarField>& Y, const label celli)
{
    //Note this assumes that mixture_fraction.update() has been called!
    auto beta_of = mixture_fraction_.get_beta();
    auto alpha = mixture_fraction_.get_alpha();

    scalar beta = 0.0; //TODO: rename!
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

bool simpleRefMapping::applyMapping(PtrList<volScalarField>& Y, const label celli) {
    return check_if_refcell(Y,celli);
}

} //namespace Foam