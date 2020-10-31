import numpy as np
import pylab as pl

def main():

    fig, ax = pl.subplots()
    plot(ax, "standard",'-')
    plot(ax, "balanced",'--')
    
    ax.legend(loc="best")
    #ax.set_ylim(1900, 2000)
    pl.savefig("comparison.png")
    pl.show()
def plot(ax, case,line):

    path = "./{0}/postProcessing/fieldMinMax1/0/fieldMinMax.dat".format(case)
    data = np.loadtxt(path)
    t = data[:,0]
    mmin = data[:,1]
    mmax = data[:,2]
    ax.plot(t, mmax, label=case,linestyle=line)


main()
