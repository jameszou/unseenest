from __future__ import division
import os, sys, random, subprocess, pandas, gzip, math, itertools
try:
    import cPickle
except ImportError:
    import pickle as cPickle
import argparse
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
def unseen_est(filename, n_samples, gridFactor, low_percentage_bound, N_max, samples_allelle_ratio):
    file = open(filename,'r')
    f = []
    for line in file:
        f.append(int(line.strip()))
    
    ########### BASIC CONSTANTS ###################
    if samples_allelle_ratio:
        n_samples = np.sum(np.array(f)*(1+np.arange(len(f)))) * samples_allelle_ratio
    maxLPIters = 1000
    xLPmax = len(f)/n_samples
    xLPmin = low_percentage_bound*1./(n_samples*100)
    
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
    print(np.sum(np.array(f)*(1+np.arange(len(f)))),n_samples)
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

def write_output(histx, xLP, outname, discrete=False, low_percentage_bound = 0.01):
    out = open(outname, 'w')
    out.write('\t'.join(['frequency', '# of variants'])+'\n')
    for i in range(len(xLP)):
        if not discrete or histx[i,0] >= low_percentage_bound:
            out.write('\t'.join([str(xLP[i]), str(histx[i,0])])+'\n')
    out.close()
    
if __name__ == '__main__':
    # python3 unseen_est.py testing_input.txt 4 tesint_output.txt -l 50 -s 1
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input file name, read the README for the correct format.')
    parser.add_argument("alleles_numbers", help="k - The total number of alleles in the dataset. This should be the number of sequenced individuals times 2",
                    type=int)
    parser.add_argument('outname', help='Output file name, read the README for the correct format.')
    parser.add_argument("-g", "--gridFactor", help="How thin should the grid be",
                    type=float, default=1.05)
    parser.add_argument("-l", "--low_percentage_bound", help="The percentage below which a frequency should be considered low",
                    type=float, default=1)
    parser.add_argument("-n", "--N_max", help="The maximum variants number",
                    type=int, default=65000000)
    parser.add_argument("-s", "--samples_allelle_ratio", help="Replace alleles_numbers by a ratio between samples and alleles",
                    type=float, default=0)
    parser.add_argument("-d", "--discrete", help="should only output numbers higher than low_percentage_bound", action="store_true")
    args = parser.parse_args()
    filename = args.filename
    n_alleles = args.alleles_numbers
    outname = args.outname
    gridFactor = args.gridFactor
    low_percentage_bound = args.low_percentage_bound
    N_max = args.N_max
    samples_allelle_ratio = args.samples_allelle_ratio
    discrete = args.discrete
    histx, xLP = unseen_est(filename, n_alleles, gridFactor, low_percentage_bound, N_max, samples_allelle_ratio)
    write_output(histx, xLP, outname, discrete, low_percentage_bound)
