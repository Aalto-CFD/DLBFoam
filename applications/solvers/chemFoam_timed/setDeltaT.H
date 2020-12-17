if (adjustTimeStep)
{
    runTime.setDeltaT(min(max(dtChem, minDeltaT), maxDeltaT));
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}
