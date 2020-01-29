import numpy as np
import scipy.optimize
from numpy import linalg as la
from scipy.special import gammaln as lgamma
from scipy.special import digamma as dgamma
from multiprocess import Pool
import time
import multiprocessing


def E_step(omega, v0, v1, theta):
    temp_dens = np.zeros((P, P, 2));

    dens1 = - 0.5 * np.log(v1 * v1 * 2 * np.pi) - (omega ** 2) / (2 * v1 * v1) + np.log(theta);
    temp_dens[:, :, 1] = dens1;
    dens0 = - 0.5 * np.log(v0 * v0 * 2 * np.pi) - (omega ** 2) / (2 * v0 * v0) + np.log(1 - theta);
    temp_dens[:, :, 0] = dens0;
    abmax = np.amax(temp_dens, axis=2);
    EZ_new = np.exp(dens1 - abmax) / (np.exp(dens0 - abmax) + np.exp(dens1 - abmax));
    np.fill_diagonal(EZ_new, 0)

    Ed_new = EZ_new / v1 / v1 + (1 - EZ_new) / v0 / v0;
    return (EZ_new, Ed_new)

def M_step_prec(omega, Ed, S, lamb):
    out = omega
    pseq = range(P)
    for p in range(P):
        remove_i = np.delete(pseq, p);
        out_mi_mi = out[remove_i];
        out_mi_mi = out_mi_mi[:, remove_i];
        Ed_i_mi = Ed[remove_i];
        Ed_i_mi = Ed_i_mi[:, p]
        Ed_i_mi = np.diag(Ed_i_mi);
        S_i_mi = S[remove_i];
        S_i_mi = S_i_mi[:, p]
        v = N / (lamb + S[p, p]);
        invsub = la.inv(out_mi_mi);
        u = -np.matmul(la.inv((lamb + S[p, p]) * invsub + Ed_i_mi), S_i_mi)
        out[p, remove_i] = u.T;
        out[remove_i, p] = u;
        univu = np.matmul(np.matmul(u.T, invsub), u);
        out[p, p] = v + univu;

    return ([out, omega])

def M_step_theta(EZ, a, b):
    EZ_sum = np.sum(EZ, axis=(0, 1)) / 2
    PP = (P * (P - 1)) / 2
    theta = (a - 1 + EZ_sum) / (a + b + PP - 2);
    return (theta)


def s(x):
    return (sum(x))

def TPR_FPR(adj_true, adj_est):
    from sklearn.metrics import confusion_matrix
    import numpy as np

    ind = np.tril_indices(P, -1)

    tn, fp, fn, tp = (confusion_matrix(adj_true[ind], adj_est[ind]).ravel()) * 1.0
    np = (np.sum(adj_true[ind])) * 1.00
    nn = sum((adj_true[ind] == 0) * 1) * 1.00
    TPR = tp / np
    FPR = fp / (fp + tn)
    return ([TPR, FPR])


def F1(adj_true, adj_est):
    from sklearn.metrics import confusion_matrix
    import numpy as np

    ind = np.tril_indices(P, -1)

    tn, fp, fn, tp = (confusion_matrix(adj_true[ind], adj_est[ind]).ravel()) * 1.0
    prec = tp / (tp + fp)
    recall = tp / (tp + fn)
    F1 = 2 * (prec * recall) / (prec + recall)
    MCC = ((tp * tn) - (fp * fn)) / (np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return ([F1, MCC])

def Performance_B(B_true,B_est):
    TP = np.sum((B_true == 1) & (B_est == 1)) * 1.0
    FP = np.sum((B_true == 0) & (B_est == 1)) * 1.0
    FN = np.sum((B_true == 1) & (B_est == 0)) * 1.0
    TN = np.sum((B_true == 0) & (B_est == 0)) * 1.0
    NN = np.sum((B_true == 0)) * 1.0
    NP = np.sum((B_true == 1)) * 1.0
    
    prec = TP / (TP + FP)
    recall = TP / (TP + FN)
    
    F1 = 2*(prec*recall) / (prec + recall)
    MCC = ((TP*TN) - (FP*FN)) / np.sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
    
    TPR = TP / (TP + FP)
    FPR = FP / (FP + NN)
    return([TPR,FPR,MCC,F1])

def Performance_Omega(adj_true,adj_est):
    TP = np.sum((adj_true == 1) & (adj_est == 1)) / 2.0
    FP = np.sum((adj_true == 0) & (adj_est == 1)) / 2.0
    FN = np.sum((adj_true == 1) & (adj_est == 0)) / 2.0
    TN = np.sum((adj_true == 0) & (adj_est == 0)) / 2.0
    NN = np.sum((adj_true == 0)) / 2.0
    NP = np.sum((adj_true == 1)) / 2.0
    
    prec = TP / (TP + FP)
    recall = TP / (TP + FN)
    
    F1 = 2*(prec*recall) / (prec + recall)
    MCC = ((TP*TN) - (FP*FN)) / np.sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
    
    TPR = TP / (TP + FP)
    FPR = FP / (FP + NN)
    return([TPR,FPR,MCC,F1])

def objective_z(z, *args):

    x = args[1]
    mu_b = args[2]
    prec = args[3]
    alpha = np.exp(z)
    A = z - mu_b
    A = np.matrix(A)
    re = - np.sum(lgamma(alpha + x) - lgamma(alpha)) + (lgamma(np.sum(alpha + x)) - lgamma(np.sum(alpha))) + np.matmul(
        np.matmul(A, prec), A.T)
    return (float(re))


def gradient_z(z, *args):
    x = args[1]
    mu_b = args[2]
    prec = args[3]
    alpha = np.exp(z)
    A = z - mu_b
    A = np.matrix(A)
    re = - ((dgamma(alpha + x) - dgamma(alpha)) + (dgamma(np.sum(alpha + x)) - dgamma(np.sum(alpha)))) * alpha + prec.dot(A.T).T
    return (re)


def M_step_Z(alpha, x, mu, mu_b, prec):
    Z_est = np.zeros((N, P))
    mybounds = []
    for i in range(P):
        mybounds.append((None, None))
    for i in range(N):
        init = Z[i,]
        arg = (alpha[i,], x[i,], mu, mu_b[i,], prec, i)
        z_est = scipy.optimize.fmin_l_bfgs_b(objective_z, fprime=gradient_z, x0=init, args=arg, approx_grad=1,
                                             bounds=mybounds)
        Z_est[i,] = z_est[0]
    return (Z_est)


def inverse_logit(x):
    return(np.exp(x) / (1 + np.exp(x)))     



def M_step_Z_parallel(args):
    mybounds = []
    for i in range(P):
        mybounds.append((None, None))
    alpha, x, mu_b, prec, init = args
    arg = (alpha, x, mu_b, prec)
    z_est = scipy.optimize.fmin_l_bfgs_b(objective_z, fprime=gradient_z, x0=init, args=arg, approx_grad=1,
                                         bounds=mybounds)
    return (z_est[0])

def VI_VS_parallel(args):
    Y, X, v1, sig_hat, mu, phi = args
    
    sig = np.ones((Q,1)) / 2500
    theta = 0.01
    xtx = X.T.dot(X)

    mu_change = 100
    while mu_change > .01:
        mu_old = np.copy(mu)
        for j in range(Q):
            X_temp = np.copy(X)
            X_l = np.delete(X_temp,j,axis = 1)
            X_j = X_temp[:,j]
            XtX = xtx[:,j]
            XtX = np.delete(XtX,j,axis = 0)
            phi_temp = np.copy(phi)
            phi_l = np.delete(phi_temp,j,axis = 0)
            mu_temp = np.copy(mu)
            mu_l = np.delete(mu_temp,j,axis = 0)

            mu[j] = (sig[j]/sig_hat) * (X.T.dot(Y)[j] - np.sum(XtX.dot(phi_l * mu_l)))
            sig[j] = sig_hat / (xtx[j,j] + 1/v1)
            phi_j = np.log(theta/(1 - theta)) + np.log(sig[j]/(v1*sig_hat)) + mu[j]**2/(2*sig[j])
            phi[j] = inverse_logit(phi_j)
            if phi_j > 709:
                phi[j] = 1.0

        mu_change = np.max(abs(mu_old - mu))
        theta = np.sum(phi) / Q
    mu = np.squeeze(mu)
    PHI = np.squeeze(phi)
    return(mu,PHI)

def SINC(x, m, v0, v1, lamb, vB, max_iters, tol_prec, tol_B, cpus):
    #######
    # x are OTU measurements where rows are the samples
    # M are the environmental factors
    # v0 is the v0 value in the precision prior
    # v1 is the v1 value in the precision prior
    # vB is the vB value in the B coefficient prior
    # max_iters
    # tol_prec is how much change precision matrix can change each iteration before stopping
    # tol_B is how much change B matrix can change each iteration before stopping
    # lamb is the prior hyperparameter on EMGS part
    ########

    # get dimensions
    global P, N, Q
    P = x.shape[1]
    N = x.shape[0]
    Q = m.shape[1]

    ## initialize values
    global Z
    Z = np.log(x + 1)
    B0 = np.mean(Z,axis = 0)
    alpha = np.exp(Z)
    B = np.zeros((P, Q))
    phi = np.ones((P,Q)) / 4
    mu_b = m.dot(B.T)
    centr = np.mean(Z, axis=0)
    S = ((Z - centr).T).dot((Z - centr))
    sig = S / N
    omega = la.inv(sig)
    sigs_j = np.zeros(P)


    ## v0 values
    pii = .5
    
    change_z = 10000
    change_prec = 10000
    change_B = 10000

    iters_total = 0

    while (change_B > tol_B or change_prec > tol_prec):
        if (iters_total >= max_iters):
            break
        iters_total += 1
        
        ## get variance of each row
        pseq = range(P)
        SIG = np.linalg.inv(omega)
        for p in range(P):
            remove_i = np.delete(pseq,p)
            r11 = SIG[remove_i]
            r11 = r11[:,remove_i]
            r12 = SIG[p,remove_i]
            sigs_j[p] = SIG[p,p] - r12.dot(np.linalg.inv(r11)).dot(r12)

        #### update B with VI
        B_old = np.copy(B)            
        pool = Pool(cpus)
        args = [(Z[:,i] - B0[i],m,vB,sigs_j[i],B[i,],phi[i,]) for i in range(P)]
        VI_update = np.asarray(pool.map(VI_VS_parallel,args))
        B = np.squeeze(VI_update[:,0,:])
        phi = np.squeeze(VI_update[:,1,:])
        pool.close()
        pool.join

        B_mult = B * (phi > .50)
        change_B = np.max(np.absolute(B_old - B))
        mu_b = m.dot(B_mult.T)

        ## E Step
        EZ, Ed = E_step(omega, v0, v1, pii)

        ## M Step

        for p in range(P):
            B0[p] = np.mean(Z[:,p] - mu_b[:,p])
            mu_b[:,p] = mu_b[:,p] + B0[p]

        ## update centr
        S = ((Z - mu_b).T).dot((Z - mu_b))

        ## update prec
        change_p = 1
        prec_old = np.copy(omega)

        for i in range(25):
            out, omega = M_step_prec(omega, Ed, S, lamb)
            pii = M_step_theta(EZ, 2, 2)
            change_prec = np.amax(np.absolute(prec_old - omega))

        ## update Z
        pool = Pool(cpus)
        Z_old = np.copy(Z)
        args = [(alpha[i,], x[i,], mu_b[i,],omega, Z[i,]) for i in range(N)]
        Z = np.asarray(pool.map(M_step_Z_parallel, args))
        change_z = np.amax(np.absolute(Z_old - Z))
        pool.close()
        pool.join
        
        print("Finished Iteration " + str(iters_total) + ": ","change in Omega,B,Z",change_prec,change_B,change_z)
        
        adj_est = EZ > .5
        spars = np.mean(adj_est)/2

    print("v0 = ", v0, "Sparsity = ", spars)
    return (omega, EZ, phi,B,iters_total)

