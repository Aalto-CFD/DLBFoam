### Test module for validating pyJac chemistry models

- Covered tests:
    - Reference data obtained by Cantera-2.5.1-dev (c.f. XX.py)

- pyjacTestMechanism:
    - Available pyJac generated C-code for a particularly chosen chemical mechanism (modified GRI30).
    - Due to modifications discussed below, this mechanism SHOULD NOT BE UTLIZED IN REAL SIMULATIONS.
    - Modifications are made to minimize error introduced by inconsistencies between pyJac/Cantera/OpenFOAM thermodynamics.
    - Agreement in chemistry/thermodynamics computations  between cantera/pyjac/OpenFOAM can be achieved to only certain accuracy due to the following reasons:
        - Different Gas constant definition: (R_ct - R_pyjac)/R_ct = 9e-07. In addition, gas constant in OpenFOAM differs from these.  
        - eval_h and cantera enthalpies reveal difference of 1.0000000623195142*sp_enth_form[i] --> which is mainly explained by differing gas constant.
        - OpenFoam utilises a different averaging for thermodynamic quantities evaluated per given species mixture compared to Cantera. See e.g. janafThermoI.H
        - While the ODE solution vector (T,Yi) is highly consistent between references, the final temperature value on OpenFOAM side is obtained via Newton iteration 
        from enthalpy w.r.t. utilized NASA (janaf) polynomials. OpenFOAM requires that the NASA polynomials share a unique Tcommon value over all species, hence typically a refitting is required, altering results from a reference. This is avoided in the presented test mechanism. 
    - Note that often users define molecular weights according to their tabulated values while Cantera defines molecular weights as a sum of element masses.