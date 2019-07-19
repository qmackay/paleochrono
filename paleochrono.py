"""
TODO: what about symbolic links in github?
TODO: extend the chronology down to the bedrock by extrapolating the accumulation
TODO: optinally use a restart file to have a bootstrap method
TODO: is there an elegant way to unpack the variables vector in the model function?
TODO: allow to save the correction vector to be able to restart while changing the resolution
TODO: include some checks for when ddelta_depth/dz>1
TODO: Delta-depth observations should be lognormal?
TODO: we should superpose two charts for ice and air ages, one for the age and
    one for the uncertainty, since the min age is not always near 0.
TODO: also compute the prior uncertainties and show them in the figures.
TODO: is there really a computation gain with the change of variable for the
    correction functions? Avoiding this change of variables would make the code
    easier to understand. I think there is no gain since solving A^-1 b when we
    have the LU factorisation of A does not cost more than computing A^-1 * b
    when we have computed A^-1.
"""

import imp
import sys
import time
import multiprocessing
import math as m
import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import leastsq, minimize
import pccfg
from pcsite import Site
from pcsitepair import SitePair

#Reload pccfg
imp.reload(pccfg)

###Registration of start time
START_TIME = time.perf_counter()


###Opening of output.txt file
OUTPUT_FILE = open(pccfg.datadir+'output.txt', 'a')

##Global
VARIABLES = np.array([])
D = {}
DC = {}



def residuals(var):
    """Calculate the residuals."""
    resi = np.array([])
    index = 0
    for i, dlab in enumerate(pccfg.list_sites):
        D[dlab].variables = var[index:index+np.size(D[dlab].variables)]
        index = index+np.size(D[dlab].variables)
        resi = np.concatenate((resi, D[dlab].residuals(D[dlab].variables)))
        for j, dlab2 in enumerate(pccfg.list_sites):
            if j < i:
                resi = np.concatenate((resi, DC[dlab2+'-'+dlab].residuals()))
    return resi

def cost_function(var):
    """Calculate the cost function terms related to a pair of sites."""
    cost = np.dot(residuals(var), np.transpose(residuals(var)))
    return cost


def deriv_res(var):
    """Calculate derivatives for each parameter using pool."""
    zeropred = residuals(var)
    derivparams = []
    results = []
    delta = m.sqrt(np.finfo(float).eps) #Stolen from the leastsq code
    #fixme: This loop is probably sub-optimal. Have a look at what does leastsq to improve this.
    for i in range(len(var)):
        copy = np.array(var)
        copy[i] += delta
        derivparams.append(copy)
#        results.append(residuals(derivparams))
    if __name__ == "__main__":
        pool = multiprocessing.Pool(pccfg.nb_nodes)
    results = pool.map(residuals, derivparams)
    derivs = [(r - zeropred)/delta for r in results]
    return derivs

##MAIN


##Initialisation
for di, dlabel in enumerate(pccfg.list_sites):

    print('Initialization of site '+dlabel)

    D[dlabel] = Site(dlabel)
    D[dlabel].model(D[dlabel].variables)
#    D[dlabel].a_init=D[dlabel].a
#    D[dlabel].lid_init=D[dlabel].lid
    D[dlabel].write_init()
#    D[dlabel].display_init()
    VARIABLES = np.concatenate((VARIABLES, D[dlabel].variables))

for di, dlabel in enumerate(pccfg.list_sites):
    for dj, dlabel2 in enumerate(pccfg.list_sites):
        if dj < di:
            print('Initialization of site pair '+dlabel2+'-'+dlabel)
            DC[dlabel2+'-'+dlabel] = SitePair(D[dlabel2], D[dlabel])
#            DC[dlabel2+'-'+dlabel].display_init()


##Optimization
START_TIME_OPT = time.perf_counter()
print('cost function: ', cost_function(VARIABLES))
if pccfg.opt_method == 'leastsq':
    print('Optimization by leastsq')
    VARIABLES, HESS, INFODICT, MESG, LER = leastsq(residuals, VARIABLES, full_output=1)
elif pccfg.opt_method == 'leastsq-parallel':
    print('Optimization by leastsq-parallel')
    VARIABLES, HESS, INFODICT, MESG, LER = leastsq(residuals, VARIABLES, Dfun=deriv_res,
                                                   col_deriv=1, full_output=1)
elif pccfg.opt_method == "L-BFGS-B":
    print('Optimization by L-BFGS-B')
    RESULT = minimize(cost_function, VARIABLES, method='L-BFGS-B', jac=False)
    VARIABLES = RESULT.x
    print('number of iterations: ', RESULT.nit)
    HESS = np.zeros((np.size(VARIABLES), np.size(VARIABLES)))
    print('Message: ', RESULT.message)
#    cost=cost_function(VARIABLES)
elif pccfg.opt_method == 'none':
    print('No optimization')
#    HESS=np.zeros((np.size(VARIABLES),np.size(VARIABLES)))
else:
    print(pccfg.opt_method, ': Optimization method not recognized.')
    sys.exit()
print('Optimization execution time: ', time.perf_counter() - START_TIME_OPT, 'seconds')
#print 'solution: ',VARIABLES
print('cost function: ', cost_function(VARIABLES))
if pccfg.opt_method != 'none' and np.size(HESS) == 1 and HESS is None:
    print('singular matrix encountered (flat curvature in some direction)')
    sys.exit()
print('Calculation of confidence intervals')
INDEXSITE = 0
for dlabel in pccfg.list_sites:
    if pccfg.opt_method == 'none':
        D[dlabel].sigma_zero()
    else:
        D[dlabel].variables = VARIABLES[INDEXSITE:INDEXSITE+np.size(D[dlabel].variables)]
        D[dlabel].hess = HESS[INDEXSITE:INDEXSITE+np.size(D[dlabel].variables),\
            INDEXSITE:INDEXSITE+np.size(D[dlabel].variables)]
        INDEXSITE = INDEXSITE+np.size(D[dlabel].variables)
        D[dlabel].sigma()

###Final display and output
print('Display of results')
for di, dlabel in enumerate(pccfg.list_sites):
#    print dlabel+'\n'
    D[dlabel].save()
    D[dlabel].figures()
    for dj, dlabel2 in enumerate(pccfg.list_sites):
        if dj < di:
#            print dlabel2+'-'+dlabel+'\n'
            DC[dlabel2+'-'+dlabel].figures()

###Program execution time
MESSAGE = 'Program execution time: '+str(time.perf_counter()-START_TIME)+' seconds.'
print(MESSAGE)
OUTPUT_FILE.write(MESSAGE)

if pccfg.show_figures:
    mpl.show()

###Closing output file
OUTPUT_FILE.close()
