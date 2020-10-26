import sys
import os
import numpy as np
from sklearn import metrics
sys.path.append("..")
from SINC_functions_tau import *

def TPR_FPR_B(phi,B):
    B_est = np.absolute(phi) > 0.050
    B_true = B != 0

    TP = np.sum(((phi > 0.50) == 1) & ((B != 0) == 1))  * 1.0
    FP = np.sum(((phi > 0.50) == 1) & ((B != 0) == 0))  * 1.0
    FN = np.sum(((phi > 0.50) == 0) & ((B != 0) == 1))  * 1.0
    TN = np.sum(((phi > 0.50) == 0) & ((B != 0) == 0))  * 1.0
    NN = np.sum(B_true == 0) * 1.0
    NP = np.sum(B_true == 1) * 1.0

    prec = TP / (TP + FP)
    recall = TP / (TP + FN)

    F1 = 2*(prec*recall) / (prec + recall)

    MCC = ((TP*TN) - (FP*FN))/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

    TPR = TP/NP
    FPR = FP/(NN)
    return(TPR,FPR,F1,MCC)

## get shapes for simulation
B_true = np.genfromtxt("./SimulatedData/Hub/B1.csv")
X_true = np.genfromtxt("./SimulatedData/Hub/x1.csv")

N = X_true.shape[0]
P = X_true.shape[1]
Q = B_true.shape[1]


reps = 40
goal = .10
adj_est_all = np.zeros((P, P, reps))
EZ_all = np.zeros((P, P, reps))
B_est_all = np.zeros((P,Q,reps))
phi_all = np.copy(B_est_all)
spars = np.zeros(reps)
v0s = np.linspace(.00001, .01, reps)
TPR = np.zeros(reps)
FPR = np.zeros(reps)
from itertools import chain

results = np.zeros((25, 5))
results_B = np.zeros((25,4))
total_iters = np.zeros(reps)

Z_track_all = []
omega_track_all = []
B_track_all = []

for sim in range(0,25):
    TPR = np.zeros(reps)
    FPR = np.zeros(reps)
    iters = np.zeros(reps)
    spars = np.zeros(reps)
    
    v1 = 10
    lamb = 150
    vB = 1
    a_pi = 2
    b_pi = 2
    a_gamma = 2
    b_gamma = 2
    max_iters = 100
    tol_prec = 0.01
    tol_elbo = 50.0
    cpus = 24

    print("**********************start iteration " + str(sim) + "*************************")
    x = np.genfromtxt("./SimulatedData/Hub/x" + str(1 + sim) + ".csv")
    m = np.genfromtxt("./SimulatedData/Hub/m" + str(1 + sim) + ".csv")
    adj_true = np.genfromtxt("./SimulatedData/Hub/adj" + str(1 + sim) + ".csv")
    B_true = np.genfromtxt("./SimulatedData/Hub/B" + str(1 + sim) + ".csv")

    for v in range(reps):
        
        v0 = v0s[v]

        phi,B,iters_total, elbo, elbo_score  = SINC_fixed_Omega(x, m, v0, v1, lamb, vB,a_gamma,b_gamma,a_pi,b_pi, max_iters, tol_prec, tol_elbo, cpus)
        
        B_est_all[:,:,v] = B
        phi_all[:,:,v] = phi

    B_edges_true = (B_true != 0)*1.0
    B_edges_est = (phi_ind > 0.50) *1
    results_B[sim,] = Performance_B(B_edges_true,B_edges_est)

    np.savetxt("./Results/CoefficientResults_Hub_Fixed_Omega.csv",results_B, delimiter=",")

