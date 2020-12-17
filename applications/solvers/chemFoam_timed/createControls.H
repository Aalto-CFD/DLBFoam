Switch adjustTimeStep(runTime.controlDict().lookup("adjustTimeStep"));
scalar maxDeltaT(runTime.controlDict().lookup<scalar>("maxDeltaT"));
scalar minDeltaT(runTime.controlDict().lookup<scalar>("minDeltaT"));

