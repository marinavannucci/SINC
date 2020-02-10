library(Rcpp,quietly = TRUE);
suppressMessages(library(RcppEigen));
library(lbfgs,quietly = TRUE);
library(QUIC,quietly = TRUE);
suppressMessages(sourceCpp("/rsrch3/home/biostatistics/ntosborne/EM-MN/seadragon/mLDM-master/mLDM_NO.cpp",verbose = FALSE));
suppressMessages(sourceCpp("/rsrch3/home/biostatistics/ntosborne/EM-MN/seadragon/mLDM-master/mLDM.cpp",verbose = FALSE))
suppressMessages(library(foreach,quietly = TRUE))
suppressMessages(library(doParallel))


TPR_FPR = function(adj_true,adj_est) {
  adj_true_vec = adj_true[lower.tri(adj_true)]
  adj_est_vec = adj_est[lower.tri(adj_est)]
  NP = sum(adj_true_vec)
  NN = sum(adj_true_vec == 0)

  TP = sum(adj_est_vec == 1 & adj_true_vec == 1)
  FN = sum(adj_est_vec == 0 & adj_true_vec == 1)
  FP = sum(adj_est_vec == 1 & adj_true_vec == 0)
  TN = sum(adj_est_vec == 0 & adj_true_vec == 0)

  prec = TP / (TP + FP)
  recall = TP / (TP + FN)

  F1 = 2*(prec*recall) / (prec + recall)

  MCC = ((TP*TN) - (FP*FN))/sqrt((TP+FP)*(TP+FN))
  MCC = MCC * (1/sqrt((TN+FP)*(TN+FN)))

  TPR = TP/NP
  FPR = FP/(FP + TN)
  return(c(TPR,FPR,F1,MCC))
}



registerDoParallel()

## mLDM
mLDM_results = c()
mLDM_results_B = c()

## what happens when you try and do parrallel with two rbinds?

mLDM.loop = foreach(sim = 1:25, .combine = rbind) %dopar% {
  X = as.matrix(read.table(paste("./SimulatedData/Random/x",i,".csv",sep = "")))
  m = as.matrix(read.table(paste("./SimulatedData/Random/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Random/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Random/B",i,".csv",sep = "")))
  
  ### mLDM

  ### get AUC of mLDM
  output = mLDM(x,m,max_iteration = 500,model_selection_num = 5)
  output_optimal = output
  lambda1_optimal = output$optimal[[6]]
  rm(output)
  output = mLDM_NO(x,m,max_iteration = 500,model_selection_num = 5,ratio2 = .999,ratio1 = .00001,lam1 = lambda1_optimal,lam2 = 10)

  results = c()

  for(i in 1:(1*10)) {
    oa = output$all[[i]]
    if(length(oa) != 0) {
      # print(c(sum(oa[[3]] != 0),oa[[5]],oa[[6]]))
      results = rbind(results,TPR_FPR(adj_true,oa[[3]] != 0))
      print(c(oa[[5]],oa[[6]]))
    }
  }

  TPR = results[,1]
  ord = order(TPR)
  TPR = c(0,results[ord,1],1)
  FPR = c(0,results[ord,2],1)

  FPR_keep = sapply(1:length(TPR),function(x) min(FPR[TPR == TPR[x]]))
  TPR_keep = sapply(1:length(FPR_keep),function(x) max(TPR[FPR == FPR[x]]))
  ord = order(FPR_keep)
  TPR_keep = TPR_keep[ord]
  FPR_keep = FPR_keep[ord]

  ind = c(TRUE,sapply(2:length(TPR_keep),function(x) TPR_keep[x] > TPR_keep[x-1]))

  TPR_keep = TPR_keep[ind]
  FPR_keep = FPR_keep[ind]

  ind = c(TRUE,sapply(2:length(TPR_keep),function(x) TPR_keep[x] > TPR_keep[x-1]))

  TPR_keep = TPR_keep[ind]
  FPR_keep = FPR_keep[ind]

  ind = c(TRUE,sapply(2:length(FPR_keep),function(x) FPR_keep[x] > FPR_keep[x-1]))

  TPR_keep = TPR_keep[ind]
  FPR_keep = FPR_keep[ind]

  auc = 0
  for(i in 1:(length(TPR_keep) - 1)) {
    b1 = TPR_keep[i]
    b2 = TPR_keep[i+1]
    h = FPR_keep[i+1] - FPR_keep[i]
    auc = auc + (b1 + b2)*h/2
  }

  adj_est = (output_optimal$optimal[[3]] != 0)
  B_est = (output_optimal$optimal[[1]] != 0)

  mLDM_results = rbind(mLDM_results,c(TPR_FPR(adj_true,adj_est),auc,error_B(B_est, B_true)))
  #  mLDM_results_B = rbind(mLDM_results_B,c(TPR_FPR(B_true,B_est)))

}


write.csv(mLDM.loop,"./Results/mLDM_results_Random.csv")
