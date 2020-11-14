import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib

fig,ax = plt.subplots()
mean = 0
for rankid,rank in enumerate(sorted(glob.glob('processor*'))):
    # Plot add, retrieve, grow and solve

    get_problem = []
    update_state = []
    balance    = []
    solve_buffer = []
    unbalance = []

    time = []
    cpu_solve = open(rank+'/loadBal/cpu_solve.out')
    lines_solve = cpu_solve.readlines()[1:]
    for x in lines_solve:
        time.append(x.split()[0])
        get_problem.append(x.split()[1])
        update_state.append(x.split()[2])
        balance.append(x.split()[3])
        solve_buffer.append(x.split()[4])
        unbalance.append(x.split()[5])

    if(rankid==0):
        size = np.size(time)-2
    time = np.array([float(i) for i in time])
    get_problem = np.array([float(i) for i in get_problem])
    update_state = np.array([float(i) for i in update_state])
    balance = np.array([float(i) for i in balance])
    solve_buffer = np.array([float(i) for i in solve_buffer])
    unbalance = np.array([float(i) for i in unbalance])         
    total = get_problem[:size]  + update_state[:size] + balance[:size] + solve_buffer[:size] + unbalance[:size]
    mean += solve_buffer[:size]
    ax.plot(solve_buffer,linewidth=0.6,label='Processor ID' if rankid==0 else "")

ax.plot(mean/np.size(glob.glob('processor*')),'k--',label='Mean',linewidth=0.9)
ax.tick_params(bottom=True,top=False,left=True,right=True,labeltop=False,labelright=False,length=2,direction='in')
ax.set_ylabel('Chemistry CPU time [s]')


ax.set_xlim([0,None])
ax.legend(loc=1,frameon=False)
fig.text(0.5, -0.01, 'Number of iterations', ha='center')
fig.tight_layout()
fig.savefig('rankbased_solve.png',bbox_inches='tight',pad_inches=0,dpi=600)
plt.show()

