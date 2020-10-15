#pragma once

namespace Foam {


void sanity_check(volScalarField& p, volScalarField& rho, PtrList<volScalarField>& Y, psiReactionThermo& thermo) {

    basicChemistryModel* model1 = new ode<StandardChemistryModel<psiReactionThermo, ThermoType>>(thermo);
    basicChemistryModel* model2 = new ode<LoadBalancedChemistryModel<psiReactionThermo, ThermoType>>(thermo);

    set_master_heavy(rho, thermo);

    for (size_t i = 0; i < 2; ++i){
        model1->solve(detail::FLOWTIMESTEP);
        model2->solve(detail::FLOWTIMESTEP);
    }



    for (label i = 0; i < model1->nSpecie(); ++i){


        const scalarField& r1 = model1->RR(i);
        const scalarField& r2 = model2->RR(i);


        scalar sum = 0.0;
        forAll(r1, celli) {

            scalar diff = r1[celli] - r2[celli];
            sum += (diff * diff);
        }

        scalar norm =std::sqrt(sum);


        if (norm > 1E-10){
            Info << "Load balanced model and standard model give different answers! Be careful!" << endl;
        }

    }


    delete model1; delete model2;

}

}