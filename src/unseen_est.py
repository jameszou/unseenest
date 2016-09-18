from __future__ import division
import os, sys, random, cPickle, sets, subprocess, pandas, gzip, math, itertools
import numpy as np
from operator import itemgetter
from scipy.stats import binom_test
from scipy.stats import chi2_contingency
from scipy.stats import entropy
from scipy.stats import poisson
from scipy.stats import binom
from scipy.stats import hypergeom
import statsmodels.api as sm
import networkx as nx
import matplotlib.pyplot as plt
import scipy.io
from cvxopt import matrix, solvers

# f is a list of fingerprint values
# n_samples is the number of alleles in the cohort
def unseen_est(filename, n_samples):
    file = open(filename,'r')
    f = []
    for line in file:
        f.append(int(line.strip()))
    
    ########### BASIC CONSTANTS ###################
    gridFactor = 1.05
    maxLPIters = 1000
    xLPmax = len(f)/n_samples
    xLPmin = 1./(n_samples*100)
    N_max = 65000000
    #N_max = 650000000
    
    ########### SETTING UP THE LP ###################
    fLP = f + [0]*int(np.ceil(np.sqrt(len(f))))
    szLPf = len(fLP)
    xLP = xLPmin*np.power(gridFactor, np.arange(0, np.ceil(np.log(xLPmax/xLPmin)/np.log(gridFactor))+1))
    szLPx = np.max(xLP.shape)
    
    ## set up the objective function
    objf = np.zeros((1, szLPx + 2*szLPf))
    objf[0, np.arange(szLPx, szLPx + 2*szLPf, 2)] = 1./np.sqrt(np.array(fLP) + 1)
    objf[0, np.arange(szLPx+1, szLPx + 2*szLPf, 2)] = 1./np.sqrt(np.array(fLP) + 1)
    
    ## set up the inequality constraints corresponding to the moment matching
    ## first 2*szLPf are for moment matching, next szLPx+2*szLPf for >=0, last for <= N_max
    A = np.zeros((2*szLPf+szLPx+2*szLPf+1, szLPx+2*szLPf))  
    b = np.zeros((2*szLPf+szLPx+2*szLPf+1, 1))
    
    rv_list = [binom(n_samples, x) for x in xLP]
    # moment matching constraints
    for i in range(szLPf):
        A[2*i, np.arange(szLPx)] = [rv.pmf(i+1) for rv in rv_list]
        A[2*i+1, np.arange(szLPx)] = -A[2*i, np.arange(szLPx)]
        A[2*i, szLPx+2*i] = -1
        A[2*i+1, szLPx+2*i+1] = -1
        b[2*i, 0] = fLP[i]
        b[2*i+1, 0] = -fLP[i]
    
    # >= 0 constraints
    for i in range(szLPx+2*szLPf):
        A[i+2*szLPf,i] = -1
        b[i+2*szLPf,0] = 0
        
    # <= N_max constraint
    A[-1,range(szLPx)] = 1
    b[-1,0] = N_max
    
        
    ## set up the equality constraints
    Aeq = np.zeros((1, szLPx+2*szLPf))
    Aeq[0, range(szLPx)] = xLP
    beq = np.sum(np.array(f)*(1+np.arange(len(f))))/n_samples
    
    ########### RUNNING THE LP ###################
    
    solvers.options['show_progress'] = False
    
    ## rescaling for better conditioning
    for j in range(np.max(xLP.shape)):
        A[:,j] = A[:,j]/xLP[j]
        Aeq[0,j] = Aeq[0,j]/xLP[j]
    
    #return objf, A, b, szLPf, szLPx, xLP
    sol = solvers.lp(matrix(objf.T), matrix(A), matrix(b), matrix(Aeq), matrix(beq))    
    #res = linprog(list(objf[0]), A_ub = A, b_ub = list(b.T[0]), A_eq = Aeq, b_eq = [beq] , options = {'maxiter': maxLPIters})
    
    ## remove the scaling
    histx = np.array(sol['x'])[0:szLPx]
    histx = [histx[i]/xLP[i] for i in range(szLPx)]
    
    return np.array(histx), xLP

def write_output(histx, xLP, outname):
    out = open(outname, 'w')
    out.write('\t'.join(['frequency', '# of variants'])+'\n')
    for i in range(len(xLP)):
        out.write('\t'.join([str(xLP[i]), str(histx[i,0])])+'\n')
    out.close()
    
if __name__ == '__main__':
    filename = sys.argv[1]
    n_alleles = int(sys.argv[2])
    outname = sys.argv[3]
    histx, xLP = unseen_est(filename, n_alleles)
    write_output(histx, xLP, outname)



    
    