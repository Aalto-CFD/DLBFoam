import os
import numpy as np
import pdb
import subprocess
odes = ['pyJac','standard']
abstol = [1e-4,1e-6]
Nsample = 2
for ode in odes:
        os.system('mkdir {}'.format(ode))
        os.system('cd {}'.format(ode))
        for tolerance in abstol:
                    case_id = ode + '/' + str(tolerance) 
                    subprocess.call('foamCloneCase base {}'.format(case_id),shell=True)
                    if(ode=='pyJac'):
                        subprocess.call('foamDictionary  -entry chemistryType/solver {}/constant/chemistryProperties -set {}'.format(case_id,'ode_pyJac'),shell=True)
                        subprocess.call('foamDictionary  -entry chemistryType/method {}/constant/chemistryProperties -set {}'.format(case_id,'loadBalanced_pyJac'),shell=True)
                        subprocess.call('foamDictionary  -entry odeCoeffs/solver {}/constant/chemistryProperties -set {}'.format(case_id,'seulex_LAPACK'),shell=True)
                    else:
                        subprocess.call('foamDictionary  -entry chemistryType/solver {}/constant/chemistryProperties -set {}'.format(case_id,'ode'),shell=True)
                        subprocess.call('foamDictionary  -entry chemistryType/method {}/constant/chemistryProperties -set {}'.format(case_id,'loadBalanced'),shell=True)
                        subprocess.call('foamDictionary  -entry odeCoeffs/solver {}/constant/chemistryProperties -set {}'.format(case_id,'seulex'),shell=True)
                    subprocess.call('foamDictionary  -entry odeCoeffs/absTol {}/constant/chemistryProperties -set {}'.format(case_id,str(tolerance)),shell=True)
                    os.chdir(case_id)
                    meantime = 0
                    for i in range(Nsample):
                        subprocess.call('chemFoam > log.chemFoam',shell=True)
                        line = subprocess.check_output(['tail', '-5', 'log.chemFoam'])
                        time = float(line.split()[6]) 
                        print("Time: {}".format(time))
                        meantime += time
                    meantime/=Nsample
                    file = open('time.dat', 'w')
                    file.write(str(meantime))
                    file.close()

                    os.chdir('../../')
        os.system('cd ../')

