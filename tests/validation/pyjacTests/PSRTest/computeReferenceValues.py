# Note: 
# - This python script computes reference results used in the OpenFOAM validation tests.
# - Cantera environment must be created to run this python script. 
# - See https://cantera.org/index.html for further information
# - The utilised mechanism is a modified GRI30 to achieve thermodynamic consistency with openfoam
 
import time
import cantera as ct
import numpy as np

mechRelPath = "../pyjacTestMechanism/mechanism.cti"

# Enthalpy at standard conditions
gasStd = ct.Solution(mechRelPath)
gasStd.TPX = 298.15, ct.one_atm, 'CH4:0.5,O2:1,N2:3.76'
r = ct.IdealGasConstPressureReactor(gasStd)
dh0 = np.sum(gasStd.standard_enthalpies_RT*gasStd.T*ct.gas_constant*(1/gasStd.molecular_weights)*gasStd.Y)
print("\nsum(Hf*Yi): " + repr(dh0) + "\n")


gas = ct.Solution(mechRelPath)
gas.TPX = 1000.0, 1.36789e+06, 'CH4:0.5,O2:1,N2:3.76'
r = ct.IdealGasConstPressureReactor(gas)
dh0 = np.sum(gas.standard_enthalpies_RT*gas.T*ct.gas_constant*(1/gasStd.molecular_weights)*gas.Y)
print("\nsum(Hf*Yi): " + repr(dh0) + "\n")

sim = ct.ReactorNet([r])
sim.verbose = False

# limit advance when temperature difference is exceeded
delta_T_max = 20.
r.set_advance_limit('temperature', delta_T_max)

states = ct.SolutionArray(gas, extra=['t'])
print('{:10s} {:10s} {:10s} {:14s}'.format(
    't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))

sim.rtol = 1e-12
sim.atol = 1e-12
tEnd = 0.07
t = time.time()
sim.advance(tEnd)
elapsed = time.time() - t

states.append(r.thermo.state, t=sim.time*1e3)
print('{:10.3e} {:10.6f} {:10.3f} {:14.10f}'.format(
        sim.time, r.T, r.thermo.P, gas.Y[gas.species_index('CH4')]))#r.thermo.u))
        
print("\n Wall clock time: " + repr(elapsed))
