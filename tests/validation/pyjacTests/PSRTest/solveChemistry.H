StopWatch   chemistryTime;  
chemistryTime.start();
dtChem = chemistry.solve(endTime);
chemistryTime.stop();
scalar Qdot = chemistry.Qdot()()[0]/rho[0];
integratedHeat += Qdot*runTime.deltaT().value();
(void)dtChem; // avoid unused error
#include "YEqn.H"
#include "hEqn.H"
#include "pEqn.H"