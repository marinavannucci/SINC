library(SpiecEasi)
library(SSLASSO)
library(MASS)
library(mSSL)

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}


se_results = c()
se_auc = c()
for(i in 1:25) {
  
  X = as.matrix(read.table(paste("./SimulatedData/Random/x",i,".csv",sep = "")))
  m = as.matrix(read.table(paste("./SimulatedData/Random/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Random/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Random/B",i,".csv",sep = "")))
  
  se = spiec.easi(as.matrix(X),nlambda = 100)
  se$refit$stars
  se_results = rbind(se_results,round(error_Omega(se$refit$stars,adj_true),4))
  
  
  TPR = FPR = 0
  for(p in 1:100) {
    adj_est = se$est$path[[p]]
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  se_auc = c(se_auc,simple_auc(TPR,FPR))

}

write.csv(cbind(se_results,mSSL_auc),"./Results/SE_results_Random.csv")

############
############

prec_mSSL = c()
B_mSSL = c()
mSSL_auc = c()

for(i in 1:25) {
  X = as.matrix(read.table(paste("./SimulatedData/Random/x",i,".csv",sep = "")))
  X = clr(X + 1)
  m = as.matrix(read.table(paste("./SimulatedData/Random/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Random/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Random/B",i,".csv",sep = "")))
  B_true = t(B_true)
  fit_mSSL_dpe <- mSSL_dpe(m,X)
  
  
  B_mSSL = rbind(B_mSSL,error_B(fit_mSSL_dpe$B, B_true))
  prec_mSSL = rbind(prec_mSSL,error_Omega(fit_mSSL_dpe$Omega, adj_true))

  TPR = FPR = 1
  for(p in 1:100) {
    adj_est = fit_mSSL_dcpe$Omega_path[,,p]
    diag(adj_est) = 0
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  mSSL_auc = c(mSSL_auc,simple_auc(sort(TPR),sort(FPR)))
}

write.csv(cbind(prec_mSSL,mSSL_auc,B_mSSL),"./Results/mSSL_results_Random.csv")

############

se_results = c()
se_auc = c()
for(i in 1:25) {
  
  X = as.matrix(read.table(paste("./SimulatedData/Cluster/x",i,".csv",sep = "")))
  m = as.matrix(read.table(paste("./SimulatedData/Cluster/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Cluster/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Cluster/B",i,".csv",sep = "")))
  
  se = spiec.easi(as.matrix(X),nlambda = 100)
  se$refit$stars
  se_results = rbind(se_results,round(error_Omega(se$refit$stars,adj_true),4))
  
  
  TPR = FPR = 0
  for(p in 1:100) {
    adj_est = se$est$path[[p]]
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  se_auc = c(se_auc,simple_auc(TPR,FPR))
  
}

write.csv(cbind(se_results,mSSL_auc),"./Results/SE_results_Cluster.csv")

############
############

prec_mSSL = c()
B_mSSL = c()
mSSL_auc = c()

for(i in 1:25) {
  X = as.matrix(read.table(paste("./SimulatedData/Cluster/x",i,".csv",sep = "")))
  X = clr(X + 1)
  m = as.matrix(read.table(paste("./SimulatedData/Cluster/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Cluster/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Cluster/B",i,".csv",sep = "")))
  B_true = t(B_true)
  fit_mSSL_dpe <- mSSL_dpe(m,X)
  
  
  B_mSSL = rbind(B_mSSL,error_B(fit_mSSL_dpe$B, B_true))
  prec_mSSL = rbind(prec_mSSL,error_Omega(fit_mSSL_dpe$Omega, adj_true))
  
  TPR = FPR = 1
  for(p in 1:100) {
    adj_est = fit_mSSL_dcpe$Omega_path[,,p]
    diag(adj_est) = 0
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  mSSL_auc = c(mSSL_auc,simple_auc(sort(TPR),sort(FPR)))
}

write.csv(cbind(prec_mSSL,mSSL_auc,B_mSSL),"./Results/mSSL_results_Cluster.csv")

#####################
se_results = c()
se_auc = c()
for(i in 1:25) {
  
  X = as.matrix(read.table(paste("./SimulatedData/Hub/x",i,".csv",sep = "")))
  m = as.matrix(read.table(paste("./SimulatedData/Hub/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Hub/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Hub/B",i,".csv",sep = "")))
  
  se = spiec.easi(as.matrix(X),nlambda = 100)
  se$refit$stars
  se_results = rbind(se_results,round(error_Omega(se$refit$stars,adj_true),4))
  
  
  TPR = FPR = 0
  for(p in 1:100) {
    adj_est = se$est$path[[p]]
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  se_auc = c(se_auc,simple_auc(TPR,FPR))
  
}

write.csv(cbind(se_results,mSSL_auc),"./Results/SE_results_Hub.csv")

############
############

prec_mSSL = c()
B_mSSL = c()
mSSL_auc = c()

for(i in 1:25) {
  X = as.matrix(read.table(paste("./SimulatedData/Hub/x",i,".csv",sep = "")))
  X = clr(X + 1)
  m = as.matrix(read.table(paste("./SimulatedData/Hub/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Hub/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Hub/B",i,".csv",sep = "")))
  B_true = t(B_true)
  fit_mSSL_dpe <- mSSL_dpe(m,X)
  
  
  B_mSSL = rbind(B_mSSL,error_B(fit_mSSL_dpe$B, B_true))
  prec_mSSL = rbind(prec_mSSL,error_Omega(fit_mSSL_dpe$Omega, adj_true))
  
  TPR = FPR = 1
  for(p in 1:100) {
    adj_est = fit_mSSL_dcpe$Omega_path[,,p]
    diag(adj_est) = 0
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  mSSL_auc = c(mSSL_auc,simple_auc(sort(TPR),sort(FPR)))
}

write.csv(cbind(prec_mSSL,mSSL_auc,B_mSSL),"./Results/mSSL_results_Hub.csv")

##################

se_results = c()
se_auc = c()
for(i in 1:25) {
  
  X = as.matrix(read.table(paste("./SimulatedData/Band/x",i,".csv",sep = "")))
  m = as.matrix(read.table(paste("./SimulatedData/Band/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Band/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Band/B",i,".csv",sep = "")))
  
  se = spiec.easi(as.matrix(X),nlambda = 100)
  se$refit$stars
  se_results = rbind(se_results,round(error_Omega(se$refit$stars,adj_true),4))
  
  
  TPR = FPR = 0
  for(p in 1:100) {
    adj_est = se$est$path[[p]]
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  se_auc = c(se_auc,simple_auc(TPR,FPR))
  
}

write.csv(cbind(se_results,mSSL_auc),"./Results/SE_results_Band.csv")

############
############

prec_mSSL = c()
B_mSSL = c()
mSSL_auc = c()

for(i in 1:25) {
  X = as.matrix(read.table(paste("./SimulatedData/Band/x",i,".csv",sep = "")))
  X = clr(X + 1)
  m = as.matrix(read.table(paste("./SimulatedData/Band/m",i,".csv",sep = "")))
  adj_true = as.matrix(read.table(paste("./SimulatedData/Band/adj",i,".csv",sep = "")))
  B_true = as.matrix(read.table(paste("./SimulatedData/Band/B",i,".csv",sep = "")))
  B_true = t(B_true)
  fit_mSSL_dpe <- mSSL_dpe(m,X)
  
  
  B_mSSL = rbind(B_mSSL,error_B(fit_mSSL_dpe$B, B_true))
  prec_mSSL = rbind(prec_mSSL,error_Omega(fit_mSSL_dpe$Omega, adj_true))
  
  TPR = FPR = 1
  for(p in 1:100) {
    adj_est = fit_mSSL_dcpe$Omega_path[,,p]
    diag(adj_est) = 0
    res = round(error_Omega(adj_est,adj_true),4)
    TPR = c(TPR,res[1] / (res[1] + res[4]))
    FPR = c(FPR,1 - res[6])
  }
  mSSL_auc = c(mSSL_auc,simple_auc(sort(TPR),sort(FPR)))
}

write.csv(cbind(prec_mSSL,mSSL_auc,B_mSSL),"./Results/mSSL_results_Band.csv")
