runTime.controlDict().lookup("adjustTimeStep") >> adjustTimeStep;
maxDeltaT = runTime.controlDict().lookup<scalar>("maxDeltaT");
minDeltaT = runTime.controlDict().lookup<scalar>("minDeltaT");

