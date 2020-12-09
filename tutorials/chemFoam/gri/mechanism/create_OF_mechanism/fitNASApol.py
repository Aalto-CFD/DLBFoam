import numpy as np
import lsqlin as fitter #3rd party matlab like lsqlin, much better than any scipy shit

# lsqlin follows the matlab style https://se.mathworks.com/help/optim/ug/lsqlin.html

# See the following links for further information on the NASA polynomial fitting procedure:
# - https://pdfs.semanticscholar.org/4920/6eb96b41fcbc8b19526b2cf2a3b10e02e0b1.pdf
# - http://shepherd.caltech.edu/EDL/publications/reprints/RefittingThermoDataNew.pdf
# - https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930003779.pdf

'''
Description below, follows http://shepherd.caltech.edu/EDL/PublicResources/sdt/thermo.html:

This function fits polynomials to specific heat, using NASA 7 term format from SP-272. Gordon and McBride 1971.
- Eq. 90: cp0/R=a1+a2T+a3*T^2+a4*T^3+a5*T^4 
- Eq. 91: H0_T/(RT)=a1+a2/2*T+a3/3*T^2+a4/4*T^3+a5/5*T^4+a6/T 
- Eq. 92: S0_T/R=a1*ln(T)+a2*T+a3/2*T^2+a4/3*T^3+a5/4*T^4+a7

Depending on the initial data, the fitting is made by either of the following methods:

1) Original data is from high-quality fit or from partition functions.
- Fit the specific heat only and use the analytical integration to obtain the additional coefficients for entropy and enthalpy.  
The specific heat fit constrains the derivative to be continuous at Tcommon and the value of 
specific heat at 298.15 K is constrained to be the specified value of the standard state (enthalpy of formation). 
The entropy and enthalpy values at the midpoints are forced to be continuous by analytical constaining.

- This procedure produces high quality fits as long as the input entropy and enthalpy are thermodynamically 
consistent with the specific heat data.  This will always be the case if the thermodynamic data have been 
derived from partition functions and presumably when refitting existing fit functions. Be cautious in refitting!

2) Original data from low-quality fit or from experiments
- Create a simultaneous least squares fit considering the whole system with prescribed constraints.


Fitting procedure: 

Solve linear constrained l2-regularized least squares by a Matlab's lsqlin equivalent algorithm. It is actually wrapper around CVXOPT QP solver.
    min_x ||C*x  - d||^2_2 + reg * ||x||^2_2
    s.t.  A * x <= b
            Aeq * x = beq
            lb <= x <= ub

NOTES:
- Note that multiplying with a fraction prior to the exponent ((1./2.)**(1./3.)*T)**3 ensures higher numerical arithmetic accuracy
'''

def fitNASAp(Trange0, Tc_i, cp_over_R, h_over_RT, s_over_R, cp0_by_R, dhf_by_R, s0_by_R, simultaneous_fit):
 
    Tcommon = Trange0[Tc_i]
    Tstd = 298.15

    T_l = Trange0[0:Tc_i+1]  
    T_h = Trange0[Tc_i:]
    T = np.concatenate((T_l,T_h))
    
    cp_over_R_L = cp_over_R[0:Tc_i+1]
    cp_over_R_H = cp_over_R[Tc_i:]

    h_over_RT_L = h_over_RT[0:Tc_i+1]
    h_over_RT_H = h_over_RT[Tc_i:]

    s_over_R_L = s_over_R[0:Tc_i+1]
    s_over_R_H = s_over_R[Tc_i:]

    Nl=len(T_l)
    Nh=len(T_h)
    N=len(T)
       
    if(not simultaneous_fit):

        print( '\nFitting cp/R only and derive the remaining coeffcients analytically:\n' )

        M = 5
        C = np.zeros((N,2*M))
        d = np.zeros(N)

        # cp equation
        for i in range(0,5):
            C[:Nl,i] = pow(T_l,i)
            C[Nl:Nl+Nh,i+M] = pow(T_h,i)       

        #The right hand side
        d[:Nl]=cp_over_R_L
        d[Nl:Nl+Nh]=cp_over_R_H

        # Constraints
        Aeq=np.zeros((5,2*M))
        beq=np.zeros(5)
        
        # cp_over_R : C0 continuity at Tcommon and equal first and last element values
        for i in range(0,5):
            Aeq[0,i] = Tcommon**i
            Aeq[0,i+M] = -Tcommon**i
            Aeq[1,i] = Tcommon**i 
            Aeq[2,i] = Tstd**i
            Aeq[3,i+M] = T[-1]**i
        beq[0] = 0.0
        beq[1] = cp_over_R[Tc_i]
        beq[2] = cp0_by_R
        beq[3] = cp_over_R_H[-1]
        
        # cp_over_R : C1 continuity at Tcommon
        Aeq[4,1] = 1.
        Aeq[4,2] = 2.*Tcommon
        Aeq[4,3] = ((3.**(1./2.))*Tcommon)**2
        Aeq[4,4] = ((4.**(1./3.))*Tcommon)**3
        Aeq[4,1+M] = -1.
        Aeq[4,2+M] = -2.*Tcommon
        Aeq[4,3+M] = -((3.**(1./2.))*Tcommon)**2
        Aeq[4,4+M] = -((4.**(1./3.))*Tcommon)**3
        
        # Find the least-squares solution
        sol=fitter.lsqlin(C, d, 0, None, None, Aeq, beq,-1e9,1e9,None,{'show_progress': False, 'abstol': 1e-12, 'reltol': 1e-8})
        coeffs_tmp=sol['x']

        # fill the coefficient matrix
        coeffs=np.zeros(14) 
        for i in range(0,5):
            coeffs[i] = coeffs_tmp[i]
            coeffs[i+7] = coeffs_tmp[i+M]

        
        # solving additional constant by means of conservation of enthalpy and entropy + ensuring C0 continuity at T=Tcommon 
        
        # Define coeff[5] in terms of enthalpy of formation at standard conditions
        coeffs[5] = dhf_by_R - ( coeffs[0]*Tstd + coeffs[1]*((1./2.)**(1./2.)*Tstd)**2 \
                + coeffs[2]*((1./3.)**(1./3.)*Tstd)**3 + coeffs[3]*((1./4.)**(1./4.)*Tstd)**4 + coeffs[4]*((1./5.)**(1./5.)*Tstd)**5 )
        # Define coeffs[5+7] by ensuring continuity at Tcommon
        h_sum_term_l = lambda T:  coeffs[0] + coeffs[1]*(1./2.)*T + coeffs[2]*((1./3.)**(1./2.)*T)**2 + coeffs[3]*((1./4.)**(1./3.)*T)**3 + coeffs[4]*((1./5.)**(1./4.)*T)**4 
        h_sum_term_h = lambda T:  coeffs[0+7] + coeffs[1+7]*(1./2.)*T + coeffs[2+7]*((1./3.)**(1./2.)*T)**2 + coeffs[3+7]*((1./4.)**(1./3.)*T)**3 + coeffs[4+7]*((1./5.)**(1./4.)*T)**4 
        coeffs[5+7] = Tcommon*( h_sum_term_l(Tcommon) + coeffs[5]/Tcommon - h_sum_term_h(Tcommon) )
        
        # Define coeff[6] in terms of entropy at standard conditions
        coeffs[6] = s0_by_R - ( coeffs[0]*np.log(Tstd) + coeffs[1]*Tstd + coeffs[2]*((1./2.)**(1./2.)*Tstd)**2 \
                + coeffs[3]*((1./3.)**(1./3.)*Tstd)**3 + coeffs[4]*((1./4.)**(1./4.)*Tstd)**4  )
        # Define coeffs[6+7] by ensuring continuity at Tcommon
        s_sum_term_l = lambda T:  coeffs[0]*np.log(T) + coeffs[1]*T + coeffs[2]*((1./2.)**(1./2.)*T)**2 + coeffs[3]*((1./3.)**(1./3.)*T)**3 + coeffs[4]*((1./4.)**(1./4.)*T)**4 
        s_sum_term_h = lambda T:  coeffs[0+7]*np.log(T) + coeffs[1+7]*T + coeffs[2+7]*((1./2.)**(1./2.)*T)**2 + coeffs[3+7]*((1./3.)**(1./3.)*T)**3 + coeffs[4+7]*((1./4.)**(1./4.)*T)**4 
        coeffs[6+7] = s_sum_term_l(Tcommon) + coeffs[6] - s_sum_term_h(Tcommon) 
        
        





    if(simultaneous_fit):

        print( '\nFitting cp/R, h/RT and s/R simultaneously:\n' )

        M = 5
        C = np.zeros((3*N,2*M))
        d = np.zeros(3*N)

        # cp equation
        for i in range(0,5):
            C[:Nl,i] = pow(T_l,i)
            C[Nl:Nl+Nh,i+M] = pow(T_h,i)           
        
        # h equation
        c1 = 1./2.
        c2 = (1./3.)**(1./2.)
        c3 = (1./4.)**(1./3.)
        c4 = (1./5.)**(1./4.)
        
        i1=Nl+Nh
        i2=i1+Nl
        C[i1:i2,0] = np.ones((1,Nl)) - Tstd/T_l
        C[i1:i2,1] = (c1*T_l)    - (Tstd/T_l)*(c1*Tstd)   
        C[i1:i2,2] = (c2*T_l)**2 - (Tstd/T_l)*(c2*Tstd)**2
        C[i1:i2,3] = (c3*T_l)**3 - (Tstd/T_l)*(c3*Tstd)**3
        C[i1:i2,4] = (c4*T_l)**4 - (Tstd/T_l)*(c4*Tstd)**4

        i3=i2
        i4=i2+Nh
        C[i3:i4,0+M] = np.ones((1,Nh)) - Tstd/T_h
        C[i3:i4,1+M] = (c1*T_h)    - (Tstd/T_h)*(c1*Tstd)
        C[i3:i4,2+M] = (c2*T_h)**2 - (Tstd/T_h)*(c2*Tstd)**2
        C[i3:i4,3+M] = (c3*T_h)**3 - (Tstd/T_h)*(c3*Tstd)**3
        C[i3:i4,4+M] = (c4*T_h)**4 - (Tstd/T_h)*(c4*Tstd)**4

        
        # s equation
        c2 = (1./2.)**(1./2.)
        c3 = (1./3.)**(1./3.)
        c4 = (1./4.)**(1./4.)

        i5=i4
        i6=i4+Nl
        C[i5:i6,0] = np.log(T_l/Tstd)
        C[i5:i6,1] = T_l - Tstd
        C[i5:i6,2] = (c2*T_l)**2 - (c2*Tstd)**2
        C[i5:i6,3] = (c3*T_l)**3 - (c3*Tstd)**3
        C[i5:i6,4] = (c4*T_l)**4 - (c4*Tstd)**4
        
        i7=i6
        i8=i6+Nh
        C[i7:i8,0+M] = np.log(T_h/Tstd)
        C[i7:i8,1+M] = T_h - Tstd
        C[i7:i8,2+M] = (c2*T_h)**2 - (c2*Tstd)**2
        C[i7:i8,3+M] = (c3*T_h)**3 - (c3*Tstd)**3
        C[i7:i8,4+M] = (c4*T_h)**4 - (c4*Tstd)**4
        
        #The right hand side
        d[:Nl]=cp_over_R_L
        d[Nl:Nl+Nh]=cp_over_R_H
        
        d[i1:i2] = h_over_RT_L - dhf_by_R/T_l
        d[i3:i4] = h_over_RT_H - dhf_by_R/T_h
        
        d[i5:i6] = s_over_R_L - s0_by_R
        d[i7:i8] = s_over_R_H - s0_by_R




        # Constraints
        Aeq=np.zeros((2,2*M))
        beq=np.zeros(2)
        
        # Note that you can add here constraints if e.g. maximum values are known etc.

        # cp_over_R : C0 continuity at Tcommon (presumably no knowledge on min/max/std values)
        for i in range(0,5):
            Aeq[0,i] = Tcommon**i
            Aeq[0,i+M] = -Tcommon**i

        # cp_over_R : C1 continuity at Tcommon
        Aeq[1,1] = 1.
        Aeq[1,2] = 2.*Tcommon
        Aeq[1,3] = ((3.**(1./2.))*Tcommon)**2
        Aeq[1,4] = ((4.**(1./3.))*Tcommon)**3
        Aeq[1,1+M] = -1.
        Aeq[1,2+M] = -2.*Tcommon
        Aeq[1,3+M] = -((3.**(1./2.))*Tcommon)**2
        Aeq[1,4+M] = -((4.**(1./3.))*Tcommon)**3
        
        # - h and C0 continuity is guaranteed after solution by analytical consideration
        # - if forcing the constraints to the linear solution (below), result are typically less good
        #   and the other coefficients do the job anyways as suggested by pen and paper

        '''       
        # h : C0 continuity 
        Aeq[2,0] = 1. - Tstd/Tcommon
        Aeq[2,M] = -(1. - Tstd/Tcommon)
        for i in range(1,5):
            cf = (1./(i+1))**(1./i)
            Aeq[2,i] = (cf*Tcommon)**i - (Tstd/Tcommon)*(cf*Tstd)**i  
            Aeq[2,i+M] = -((cf*Tcommon)**i - (Tstd/Tcommon)*(cf*Tstd)**i)  
        
        # s : C0 continuity 
        
        Aeq[3,0] = np.log(Tcommon/Tstd)
        Aeq[3,M] = -np.log(Tcommon/Tstd)
        for i in range(1,5):
            cf = (1./i)**(1./i)
            Aeq[3,i] = (cf*Tcommon)**i - (cf*Tstd)**i
            Aeq[3,i+M] = -(cf*Tcommon)**i + (cf*Tstd)**i
        '''


        # Find the least-squares solution
        sol=fitter.lsqlin(C, d, 0, None, None, Aeq, beq,-1e9,1e9,None,{'show_progress': False, 'abstol': 1e-12, 'reltol': 1e-8})
        coeffs_tmp=sol['x']

        # fill the coefficient matrix
        coeffs=np.zeros(14) 
        for i in range(0,5):
            coeffs[i] = coeffs_tmp[i]
            coeffs[i+7] = coeffs_tmp[i+M]
        
        # solving additional constant by means of conservation of enthalpy and entropy + ensuring C0 continuity at T=Tcommon 
        
        # Define coeff[5] in terms of the enthalpy of formation at standard conditions
        coeffs[5] = dhf_by_R - ( coeffs[0]*Tstd + coeffs[1]*((1./2.)**(1./2.)*Tstd)**2 \
                + coeffs[2]*((1./3.)**(1./3.)*Tstd)**3 + coeffs[3]*((1./4.)**(1./4.)*Tstd)**4 + coeffs[4]*((1./5.)**(1./5.)*Tstd)**5 )
        
        # Define coeffs[5+7] by ensuring continuity at Tcommon
        h_sum_term_l = lambda T:  coeffs[0] + coeffs[1]*(1./2.)*T + coeffs[2]*((1./3.)**(1./2.)*T)**2 + coeffs[3]*((1./4.)**(1./3.)*T)**3 + coeffs[4]*((1./5.)**(1./4.)*T)**4 
        h_sum_term_h = lambda T:  coeffs[0+7] + coeffs[1+7]*(1./2.)*T + coeffs[2+7]*((1./3.)**(1./2.)*T)**2 + coeffs[3+7]*((1./4.)**(1./3.)*T)**3 + coeffs[4+7]*((1./5.)**(1./4.)*T)**4 
        coeffs[5+7] = Tcommon*( h_sum_term_l(Tcommon) + coeffs[5]/Tcommon - h_sum_term_h(Tcommon) )

        # Define coeff[6] in terms of entropy at standardconditions
        coeffs[6] = s0_by_R - ( coeffs[0]*np.log(Tstd) + coeffs[1]*Tstd + coeffs[2]*((1./2.)**(1./2.)*Tstd)**2 \
                + coeffs[3]*((1./3.)**(1./3.)*Tstd)**3 + coeffs[4]*((1./4.)**(1./4.)*Tstd)**4  )
        # Define coeffs[6+7] by ensuring continuity at Tcommon
        s_sum_term_l = lambda T:  coeffs[0]*np.log(T) + coeffs[1]*T + coeffs[2]*((1./2.)**(1./2.)*T)**2 + coeffs[3]*((1./3.)**(1./3.)*T)**3 + coeffs[4]*((1./4.)**(1./4.)*T)**4 
        s_sum_term_h = lambda T:  coeffs[0+7]*np.log(T) + coeffs[1+7]*T + coeffs[2+7]*((1./2.)**(1./2.)*T)**2 + coeffs[3+7]*((1./3.)**(1./3.)*T)**3 + coeffs[4+7]*((1./4.)**(1./4.)*T)**4 

        coeffs[6+7] = s_sum_term_l(Tcommon) + coeffs[6] - s_sum_term_h(Tcommon) 
        

    return coeffs
