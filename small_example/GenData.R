library(mvtnorm)
library(huge)
library(MCMCpack)

P = 25
N = 350
Q = 15
theta = huge.generator(d = P,graph = "random",prob = 0.05)
prec = theta$omega  

m = matrix(rnorm(N*Q),nrow = N)
m = scale(m)

B = matrix(0,nrow = P,ncol = Q)
for(p in 1:P) {
  for(q in 1:Q) {
    if(runif(1) > .85) {
      B[p,q] = runif(1,.5,1)
      if(runif(1) > .5) {
        B[p,q] = B[p,q]*(-1)
      }
    }
  }
}
mu = m %*% t(B)
B0 = sample(c(runif(160,2,4),runif(40,5,6)),P)
# B0 = sample(c(rep(1,1000),rep(8,100)),P)

# Z = apply(mu,1,function(x) mvtnorm::rmvnorm(1,x,solve(prec)))
# Z = t(Z)
Z = mvtnorm::rmvnorm(N,B0,solve(prec))
alpha = exp(Z + mu)
h = apply(alpha,1,function(x) MCMCpack::rdirichlet(1,x))
h = t(h)
# x = apply(h,1,function(x) rmultinom(1,sample(pops,1),x))
x = apply(h,1,function(x) rmultinom(1,round(rnorm(1,3000,500)),x))
x = t(round(x/1))

# se = spiec.easi(as.matrix(x), nlambda=100)
adj_true = as.matrix(theta$theta)

###################################################################
###################################################################
write.table(x,"/home/nathan/SINC/small_example/x.csv",col.name = F,row.name = F)
write.table(m,"/home/nathan/SINC/small_example/m.csv",col.name = F,row.name = F)
write.table(adj_true,"/home/nathan/SINC/small_example/adj_true.csv",col.name = F,row.name = F)
write.table(B,"/home/nathan/SINC/small_example/B_true.csv",col.name = F,row.name = F)
write.table(prec,"/home/nathan/SINC/small_example/prec_true.csv",col.name = F,row.name = F)
###################################################################
###################################################################
