rm(list=ls())
library(MASS)
# read source files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # setting the workspace to source file location
source("joint_selection_new_GitHub.R")

####################################################
# Set dimension
####################################################
n = 100
p = 150




####################################################
# Generation of covariate matrix X
####################################################
set.seed(12)
Omega0 = matrix(0, p,p)

#Senario 1
Omega0[1:10,1] = Omega0[1,1:10] = 0.3



diag(Omega0) = 1
library(MASS)
X = mvrnorm(n,rep(0, p),Sigma = solve(Omega0))
G0 = 1*(Omega0 != 0)
diag(G0) = 0




##################################################
# Generation of coefficient vector \beta
##################################################

beta0 = matrix(0, nrow=p, ncol=1)
t = 10
beta0[,1] = c(rep(3,t), rep(0,p-t))
gam0.loc = 1:t
gam0.len = t
gamma0 = c(rep(1,gam0.len), rep(0, p-gam0.len))

####################################################
# Generate a data vector Y : n X 1
####################################################

set.seed(12)

#X = scale(X)
Y = rlogis(n, location = X %*% beta0)
Y = ifelse(Y > 0, 1, 0) 

####################################################
# Hyperparameters
####################################################
a.hp = 2.5 
b.hp = 0.5
qn = 0.005

tau2 = 1


r = 1e-4
s = 1e-8



niter = 2000
nburn = 2000


init.gam = rep(0, p)

sd = 12

set.seed(sd)

##################################################
# Start joint selection using JSSL
##################################################
time = Sys.time()
source("joint_selection_new.R")
res = BSSC(Y, X, qn, init.A=NULL, init.gam=NULL, Rj=NULL, a.hp, b.hp, tau2, r, s, niter, nburn)
elapsed = Sys.time() - time
elapsed



################################################## 
# variable selection result
################################################## 
post.ind = as.numeric(colMeans(res$gamma.mat[2000:4000,])>0.5)
summary.joint(gamma0, post.ind)

################################################## 
# inverse covariance estimation result
################################################## 
post.Omega = apply(res$Omega.mat[,,-(1:nburn)], c(1, 2), mean, na.rm = TRUE)
post.Omega[which(abs(post.Omega) < 10^(-4))] = 0
EvaluationNorm(Omega0, post.Omega)

################################################## 
# graph selection result
##################################################
post.G = apply(res$G.mat[,,-(1:nburn)], c(1, 2), mean, na.rm = TRUE)
post.G[which(post.G > 0.5)] = 1
post.G[which(post.G < 0.5)] = 0
post.G[upper.tri(post.G)] = 0
Evaluation.G(G0, post.G)