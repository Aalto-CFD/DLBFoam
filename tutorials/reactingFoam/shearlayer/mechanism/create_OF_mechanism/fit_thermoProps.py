#####################################################################################################################################
## Copyright <2015-2019> <Aalto University, Kahila & Vuorinen>                                                                     ##
##                                                                                                                                 ##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated                    ##
## documentation files (the "Software"), to deal in the Software without restriction, including without limitation the             ##
## rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit          ##
## persons to whom the Software is furnished to do so, subject to the following conditions:                                        ##
##                                                                                                                                 ## 
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  ##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES ##
## OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE ##
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR  ##
## IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                                   ##
#####################################################################################################################################

#from IPython import get_ipython
#get_ipython().magic('reset -sf')

import cantera as ct
import numpy as np
import os
import shutil
import writing_routines as wr
import fitNASApol
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import time
import sys



def sutherlandFormula(x,T,y):
    return  x[0]*np.sqrt(T)/(1.0 + x[1]/T)  - y

# The Euken formulation used in Openfoam
def sutherlandKappa(As,Ts,mu,cv_mass,cv_mole,R):
    return mu*cv_mass*(1.32 + 1.77*R/cv_mole)



# -------------------------- #
# -- Choose the mechanism -- #
# -------------------------- #
# --------------------------------------------- #
gas = ct.Solution('../chemkin/chem.cti')

# --------------------------------------------- #

#Fit parameters
Tlow = 300.
Thigh = 3000.
nT = round((Thigh-Tlow)/25)
Trange = np.linspace(Tlow,Thigh,nT+1)

Tcommon = 1000.0
Tcommon_i = np.where(Trange==Tcommon)[0][0]

logTypeFitForTransport=1
degree = 3 #CK type fit

errLim = 1e-2; # maximum total L2 error in the cp/R,h/RT,s/R fitting



# --------------------------------------------- #
# Set up the 1d-flame IdealGasMix object used to compute mixture properties
p0 = ct.one_atm
width=0.03
#gas.TPX = 280.0, p0, reactants0
gas.TP = 280.0, p0

f = ct.FreeFlame(gas, width=width)

# Variables describing the original data obtained from Cantera
mu = np.zeros((len(Trange)))
kappa = np.zeros((len(Trange)))
cp_over_R = np.zeros((len(Trange)))
cv_mass = np.zeros((len(Trange)))
cv_mole = np.zeros((len(Trange)))
h_over_RT = np.zeros((len(Trange)))
s_over_R = np.zeros((len(Trange)))
# --------------------------------------------- #


# --------------------------------------------- #
try:  
    os.makedirs('foam_mech_files/',exist_ok=True)
except OSError:  
    None    
# File name you want to give for the openfoam style thermodynamics properties
thermoFN= 'foam_mech_files/thermo.foam'
try:
    os.remove(thermoFN) 
except OSError:
    None
try:  
    shutil.rmtree('fit_validation_figs/')
except OSError:  
    None
try:  
    os.makedirs('fit_validation_figs/')
except OSError:  
    None

# Add the chem.foam include to the file first
with open(thermoFN,'a') as output:
    
    output.write('#include "chem.foam"')
    output.write('\n')

#This would write the species names list in OF style. However, we want to include that into the reaction file in pyJac style. 
#See makeMechFileForOF.sh at pyjac folder
#wr.write_sp_list(gas.species_names,thermoFN)
# --------------------------------------------- #

for sp_i in gas.species_names:
    print("\n\n\n")
    print("########################")
    print(' Species: %s' % sp_i) 
    print("########################")
   
    reactants = sp_i + ':1.0'
    
    j=0
    for Tin in Trange:

        gas.TPX = Tin, p0, reactants
        f = ct.FreeFlame(gas, width=width)
        r = ct.IdealGasConstPressureReactor(gas)

        mu[j] = gas.viscosity
        kappa[j] = gas.thermal_conductivity
        # divide by gas constant according to NASA JANAF definitions
        cp_over_R[j] = gas.cp_mole/(ct.gas_constant)

        cv_mass[j] = gas.cv_mass
        cv_mole[j] = gas.cv_mole

        h_over_RT[j]  = gas.enthalpy_mole/(ct.gas_constant*Tin)
        s_over_R[j]  = gas.entropy_mole/(ct.gas_constant)
        j=j+1

    if( (gas.species(gas.species_index(sp_i)).thermo.min_temp > Tlow) or (gas.species(gas.species_index(sp_i)).thermo.max_temp < Thigh) ):
        print("\n\nWARNING: Given Tmin/Tmax are out of bounds!!!!!\n\n")
        print( repr(gas.species(gas.species_index(sp_i)).thermo.min_temp) + ' vs. ' + repr(Tlow) ) 
        print( repr(gas.species(gas.species_index(sp_i)).thermo.max_temp) + ' vs. ' + repr(Thigh) ) 
        #sys.exit(0)





    # -------------------------------------- #
    #  --  Fit the transport properties  --  #
    # -------------------------------------- #
    print("\nFitting transport properties:")
    print("-----------------------------")
    if(logTypeFitForTransport):
        mu_forFit = np.log(mu)
        kappa_forFit = np.log(kappa)
        Trange_forFit = np.log(Trange)
    else:
        mu_forFit = mu
        kappa_forFit = kappa
        Trange_forFit = Trange
    
    mu_fit_coeffs = np.polyfit(Trange_forFit, mu_forFit, degree)
    kappa_fit_coeffs = np.polyfit(Trange_forFit, kappa_forFit, degree)
    
    #report the error
    if(logTypeFitForTransport):
        fitFunc_mu = np.poly1d(mu_fit_coeffs)    
        fitted_mu = np.exp(fitFunc_mu(np.log(Trange)))
        fitFunc_kappa = np.poly1d(kappa_fit_coeffs)    
        fitted_kappa = np.exp(fitFunc_kappa(np.log(Trange)))
    else:
        fitFunc_mu = np.poly1d(mu_fit_coeffs)
        fitted_mu = fitFunc_mu(Trange)
        fitFunc_kappa = np.poly1d(kappa_fit_coeffs)
        fitted_kappa = fitFunc_kappa(Trange)
    
    #Compute L2 norm
    L2e_mu = sum(np.power((fitted_mu - mu),2))/(sum(np.power(mu,2)))
    L2e_kappa = sum(np.power((fitted_kappa - kappa),2))/(sum(np.power(kappa,2)))
    
    #Make fir for sutherland's coefficients
    x0 = np.array([1.0, 1.0]) #initial guess
    res_lsq = least_squares(sutherlandFormula, x0, args=(Trange, mu))
    #Compute L2 norm for sutherland
    fitted_mu_sutherland = res_lsq.x[0]*np.sqrt(Trange)/(1.0 + res_lsq.x[1]/Trange) 
    L2e_mu_sutherland = sum(np.power((fitted_mu_sutherland - mu),2))/(sum(np.power(mu,2)))
    fitted_kappa_sutherland = sutherlandKappa(res_lsq.x[0],res_lsq.x[1],fitted_mu,cv_mass,cv_mole,ct.gas_constant)
    L2e_kappa_sutherland = sum(np.power((fitted_kappa_sutherland - kappa),2))/(sum(np.power(kappa,2)))

    print('\nL2 errors of the transport fits:')
    print('mu: ' + repr(np.sqrt(L2e_mu)))
    print('Sutherland mu: ' + repr(np.sqrt(L2e_mu_sutherland)))
    print('kappa: ' + repr(np.sqrt(L2e_kappa)))
    print('Sutherland (=Euken) kappa: ' + repr(np.sqrt(L2e_kappa_sutherland)))

    # plot fits
    fig1=plt.figure(num=1,figsize=(7.5,10))
    ax4 = plt.subplot(411)
    plt.plot(Trange,mu,'-',color='r',label='Orig.')
    plt.plot(Trange,fitted_mu,'--',color='b', label='Fit')
    ax4.set_ylabel('mu')
    plt.legend(loc=4)
    ax5 = plt.subplot(412)
    plt.plot(Trange,mu,'-',color='r')
    plt.plot(Trange,fitted_mu_sutherland,'--',color='b')
    ax5.set_ylabel('mu Suth.')
    ax6 = plt.subplot(413)
    plt.plot(Trange,kappa,'-',color='r')
    plt.plot(Trange,fitted_kappa,'--',color='b')
    ax6.set_ylabel('kappa')
    ax7 = plt.subplot(414)
    plt.plot(Trange,kappa,'-',color='r')
    plt.plot(Trange,fitted_kappa_sutherland,'--',color='b')
    ax7.set_ylabel('kappa Euken')

    fig1.savefig('fit_validation_figs/'+sp_i+'_trans.png', bbox_inches='tight')
    fig1.clf()


    # --------------------------------------------------------- #
    #  --  Fit the thermodynamical properties: cp, h and s  --  #
    # --------------------------------------------------------- #

    coeffs=np.zeros(14)
    # Define functions to be used later
    cpFunc_low = lambda T: coeffs[0] + coeffs[1]*T + coeffs[2]*pow(T,2) + coeffs[3]*pow(T,3) + coeffs[4]*pow(T,4) 
    cpFunc_high = lambda T: coeffs[7] + coeffs[8]*T + coeffs[9]*pow(T,2) + coeffs[10]*pow(T,3) + coeffs[11]*pow(T,4) 
    cpdTFunc_low = lambda T: coeffs[1] + 2.*coeffs[2]*T + 3.*coeffs[3]*pow(T,2) + 4.*coeffs[4]*pow(T,3) 
    cpdTFunc_high = lambda T: coeffs[8] + 2.*coeffs[9]*T + 3.*coeffs[10]*pow(T,2) + 4.*coeffs[11]*pow(T,3) 
    hFunc_low = lambda T: coeffs[0] + coeffs[1]*T/2 + coeffs[2]*pow(T,2)/3 + coeffs[3]*pow(T,3)/4 + coeffs[4]*pow(T,4)/5 + coeffs[5]/T
    hFunc_high = lambda T: coeffs[7] + coeffs[8]*T/2 + coeffs[9]*pow(T,2)/3 + coeffs[10]*pow(T,3)/4 + coeffs[11]*pow(T,4)/5 + coeffs[12]/T
    sFunc_low = lambda T: coeffs[0]*np.log(T) + coeffs[1]*T + coeffs[2]*pow(T,2)/2 + coeffs[3]*pow(T,3)/3 + coeffs[4]*pow(T,4)/4 + coeffs[6]
    sFunc_high = lambda T: coeffs[7]*np.log(T) + coeffs[8]*T + coeffs[9]*pow(T,2)/2 + coeffs[10]*pow(T,3)/3 + coeffs[11]*pow(T,4)/4 + coeffs[13]

    # Check the quality of the existing fit (w.r.t Tcommon)
    low_T_coeffs = gas.species(gas.species_index(sp_i)).thermo.coeffs[8:15]
    high_T_coeffs = gas.species(gas.species_index(sp_i)).thermo.coeffs[1:8]
    Tcommon_orig = gas.species(gas.species_index(sp_i)).thermo.coeffs[0]
    coeffs[:7] = low_T_coeffs
    coeffs[7:14] = high_T_coeffs

    desired_Tcommon = gas.species(gas.species_index(sp_i)).thermo.coeffs[0] == Tcommon

    # relative error limits for consistency check are chosen by experience: 
    # - C0 continuity: 1e-6 error is acceptable 
    # - C1 continuityi (not mandatory but recommended): For some mechanism fits C1 continuity is not guaranteed --> a loose tolerance is chosen
    C0_at_Tc_cp = np.abs(cpFunc_low(Tcommon_orig) - cpFunc_high(Tcommon_orig))/np.abs(cpFunc_low(Tcommon_orig)) < 1e-6
    C1_at_Tc_cp = np.abs(cpdTFunc_low(Tcommon_orig) - cpdTFunc_high(Tcommon_orig))/np.abs(cpdTFunc_low(Tcommon_orig)) < 0.01 
    C0_at_Tc_h = np.abs(hFunc_low(Tcommon_orig) - hFunc_high(Tcommon_orig))/np.abs(hFunc_low(Tcommon_orig)) < 1e-6
    C0_at_Tc_s = np.abs(sFunc_low(Tcommon_orig) - sFunc_high(Tcommon_orig))/np.abs(sFunc_low(Tcommon_orig)) < 1e-6
    
    good_fit = C0_at_Tc_cp and C0_at_Tc_h and C0_at_Tc_s and C1_at_Tc_cp
    
    # Make a new fit if Tcommon is not the one wanted or if the original fit is poor
    if( not (desired_Tcommon and good_fit) ):
        print("\nRefitting thermodynamics.")
        print("------------------------")
        #Fitting via solving a constrained linear least squares solution yielding the NASA polynomials   
        # - for definition of simultaneous_fit, see fitNasapol.py
        if(good_fit):
            simultaneous_fit = False
            print('Refitting due to different Tcommon.')
            print("Tcommon: " + repr(Tcommon_orig) + ' vs. ' + repr(Tcommon))
        else:
            simultaneous_fit = True
            print('Original polynomials are not continuous at Tcommon :')
            print("Tcommon: " + repr(Tcommon_orig) + ' vs. ' + repr(Tcommon))
            print("C0 err of cp/R: "+repr(np.abs(cpFunc_low(Tcommon_orig) - cpFunc_high(Tcommon_orig))/np.abs(cpFunc_low(Tcommon_orig))))
            print("C1 err of cp/R: "+repr(np.abs(cpdTFunc_low(Tcommon_orig) - cpdTFunc_high(Tcommon_orig))/np.abs(cpdTFunc_low(Tcommon_orig))))
            print("C0 err of h/RT: "+repr(np.abs(hFunc_low(Tcommon_orig) - hFunc_high(Tcommon_orig))/np.abs(hFunc_low(Tcommon_orig))))
            print("C0 err of s/R: "+repr(np.abs(sFunc_low(Tcommon_orig) - sFunc_high(Tcommon_orig))/np.abs(sFunc_low(Tcommon_orig))))

        # - enthalpy of formation and entropy at standard conditions are required by definition
        cp0_over_R =  gas.species(gas.species_index(sp_i)).thermo.cp(298.15)/ct.gas_constant # this calls molar cp <==> consistent
        dhf_by_R =  gas.species(gas.species_index(sp_i)).thermo.h(298.15)/ct.gas_constant
        s0_by_R =  gas.species(gas.species_index(sp_i)).thermo.s(298.15)/ct.gas_constant

        coeffs = fitNASApol.fitNASAp(Trange,Tcommon_i,cp_over_R,h_over_RT,s_over_R,cp0_over_R,dhf_by_R,s0_by_R,simultaneous_fit) 
 
        low_T_coeffs = coeffs[:7]
        high_T_coeffs = coeffs[7:]
            
        # Check error
        T_l = Trange[0:Tcommon_i+1]
        T_h = Trange[Tcommon_i:]
        cp_over_R_tot = np.concatenate((cp_over_R[0:Tcommon_i+1],cp_over_R[Tcommon_i:]))
        h_over_RT_tot = np.concatenate((h_over_RT[0:Tcommon_i+1],h_over_RT[Tcommon_i:]))
        s_over_R_tot = np.concatenate((s_over_R[0:Tcommon_i+1],s_over_R[Tcommon_i:]))
        
        fitted_cp_over_R = np.concatenate((cpFunc_low(T_l),cpFunc_high(T_h)))
        fitted_h_over_RT = np.concatenate((hFunc_low(T_l),hFunc_high(T_h)))
        fitted_s_over_R = np.concatenate((sFunc_low(T_l),sFunc_high(T_h)))

        cp_over_R_error_Tc = np.abs(cpFunc_high(Tcommon) - cpFunc_low(Tcommon))
        cpdT_over_R_error_Tc = np.abs(cpdTFunc_high(Tcommon) - cpdTFunc_low(Tcommon))
        h_over_RT_error_Tc = np.abs(hFunc_high(Tcommon) - hFunc_low(Tcommon))
        s_over_R_error_Tc = np.abs(sFunc_high(Tcommon) - sFunc_low(Tcommon))
        h_over_RT_error_Tstd = np.abs( hFunc_low(298.15) - gas.species(gas.species_index(sp_i)).thermo.h(298.15)/(ct.gas_constant * 298.15) )
        s_over_R_error_Tstd = np.abs(  sFunc_low(298.15) - gas.species(gas.species_index(sp_i)).thermo.s(298.15)/ct.gas_constant  )

        L2e_cp_over_R = np.sqrt(sum(np.power((fitted_cp_over_R - cp_over_R_tot),2))/(sum(np.power(cp_over_R_tot,2))))
        L2e_h_over_RT = np.sqrt(sum(np.power((fitted_h_over_RT - h_over_RT_tot),2))/(sum(np.power(h_over_RT_tot,2))))
        L2e_s_over_R = np.sqrt(sum(np.power((fitted_s_over_R - s_over_R_tot),2)))/np.sqrt((sum(np.power(s_over_R_tot,2))))
        
        print('\nChecking the quality of the new fit:')
        print('------------------------------------')
        print( 'Thermodynamics consistency errors at Tcommon / T=298.15K' )
        print( 'cp/R: '+repr( cp_over_R_error_Tc ) )
        print( 'dcpdT/R: '+repr( cpdT_over_R_error_Tc ) )
        print( 'h/RT: '+repr( h_over_RT_error_Tc ) + ' / ' +repr( h_over_RT_error_Tstd ) ) 
        print( 's/R: '+repr( s_over_R_error_Tc ) + ' / ' + repr( s_over_R_error_Tstd ) )
        print('\nL2 errors of thermodynamics')
        print('cp/R: ' + repr(L2e_cp_over_R))
        print('h/RT: ' + repr(L2e_h_over_RT))
        print('s/R: ' + repr(L2e_s_over_R))

        if(L2e_cp_over_R>errLim or L2e_h_over_RT>errLim or L2e_s_over_R>errLim):
            print("\nERROR: Too large error in thermodynamical fit (>"+repr(errLim)+")")
            sys.exit(0)
        
        '''
        print '\ncp_over_R error at Tmin: '+repr( np.abs(cp_over_R[0] - fitted_cp_over_R[0] ))
        print 'h error at Tmin: '+repr(np.abs(  h_over_RT[0] - fitted_h_over_RT[0]  ))
        print 's error at Tmin: '+repr(np.abs(  s_over_R[0] - fitted_s_over_R[0]  ))
        print 'cp_over_R error at Tmax: '+repr( np.abs(  cp_over_R[-1] - fitted_cp_over_R[-1]  ))
        print 'h error at Tmax: '+repr(np.abs(  h_over_RT[-1] - fitted_h_over_RT[-1]  ))
        print 's error at Tmax: '+repr(np.abs(  s_over_R[-1] - fitted_s_over_R[-1]  ))
        '''
    
        fig2=plt.figure(num=2,figsize=(7.5,10))
        Tt = np.concatenate((T_l,T_h))
        ax1 = plt.subplot(311)
        plt.plot(Tt,cp_over_R_tot,'-',color='r',label='Orig.')
        plt.plot(Tt,fitted_cp_over_R,'--',color='b',label='Fit')
        ax1.set_ylabel('cp/R')
        plt.legend(loc=4)
        ax2 = plt.subplot(312)
        plt.plot(Tt,h_over_RT_tot,'-',color='r')
        plt.plot(Tt,fitted_h_over_RT,'--',color='b')
        ax2.set_ylabel('h/RT')
        ax3 = plt.subplot(313)
        plt.plot(Tt,s_over_R_tot,'-',color='r')
        plt.plot(Tt,fitted_s_over_R,'--',color='b')
        ax3.set_ylabel('s/R')
        
        fig2.savefig('fit_validation_figs/'+sp_i+'_thermo.png', bbox_inches='tight')
        fig2.clf()

    #plt.show()

    #Write the thermo output in openfoam format
    with open(thermoFN,'a') as output:
        
        output.write(sp_i+'\n{\n')
        output.write('\t')
        ############################################################################# 
        output.write('specie\n\t{\n')
        MW = gas.molecular_weights[gas.species_index(sp_i)]
        output.write('\t\tnMoles \t 1;\n')
        output.write('\t\tmolWeight \t'+repr(MW)+';')
        output.write('\n\t}\n\n')
        #############################################################################

        ############################################################################# 
        output.write('\tthermodynamics')
        #############################################################################         
        output.write('\n\t{\n')
        output.write('\t\tTlow\t\t'+repr(Tlow)+';\n')
        output.write('\t\tThigh\t\t'+repr(Thigh)+';\n')       
        output.write('\t\tTcommon\t\t'+repr(Tcommon)+';\n')
        output.write('\t\tlowCpCoeffs\t(\t' )
        for wi in range(0,7): #NASA pol has 7 coeffs
            output.write(repr(low_T_coeffs[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t\thighCpCoeffs\t(\t' )
        for wi in range(0,7):#NASA pol has 7 coeffs
            output.write(repr(high_T_coeffs[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')
        
        ############################################################################# 
        output.write('\ttransport \n\t{\n')
        #############################################################################
        output.write('\t\tAs\t'+repr(res_lsq.x[0])+';\n' )
        output.write('\t\tTs\t'+repr(res_lsq.x[1])+';\n' )
        if(logTypeFitForTransport):
            output.write('\t\tmuLogCoeffs<'+repr(degree+1)+'>\t(\t' )
        else:
            output.write('\t\tmuCoeffs<'+repr(degree+1)+'>\t(\t' )
        for wi in range(degree,-1,-1):
            output.write(repr(mu_fit_coeffs[wi])) 
            output.write(' ')
        output.write(' );\n')
        #############################################################################
        if(logTypeFitForTransport):
            output.write('\t\tkappaLogCoeffs<'+repr(degree+1)+'>\t(\t' )
        else:
            output.write('\t\tkappaCoeffs<'+repr(degree+1)+'>\t(\t' )
        for wi in range(degree,-1,-1):
            output.write(repr(kappa_fit_coeffs[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')
        ############################################################################# 

        ############################################################################# 
        output.write('\telements')
        #############################################################################      
        output.write('\n\t{\n')   
        # Variables for the elemental composition entry:
        for elem_i in gas.element_names:
            na = gas.n_atoms(sp_i,elem_i)
            if(na > 0):
                output.write("\t\t" + elem_i+"\t"+repr(int(na))+';\n')                
        #############################################################################
        output.write('\t}\n}\n\n')

