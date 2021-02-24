
if (constProp != "temperature")
{
    volScalarField& h = thermo.he();
    
    if (constProp != "pressure")
    {
        //h[0] = u0 + p[0]/rho[0] + integratedHeat;
        FatalErrorIn
        ("hEqn.H:") 
            << "Only constant pressure is supported by pyJacChemistryModel"  
            << exit(FatalError);
    }
    else
    {
        h[0] = h0 + integratedHeat;
    }

    thermo.correct();
}
