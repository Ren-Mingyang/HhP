################################################################################
# Codes for conducting simulation studies
################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(Matrix)
library(MASS)
library(fmrs)
source("function.R")
source("sim-func.R")


n <- 120                         # The sample size
q <- 5                           # The dimension of x
p <- 10                          # The dimension of z
pp <- q+p                        # The dimension of each sample
S_E <- 2                         # The dimension of non-zero coefficients of x
S_GE <- 4                        # The dimension of non-zero coefficients of z
E.cat <- 2                       # The dimension of non-zero coefficients of z
beta00 <- 2                      # The signal intensity of non-zero coefficients
balance <- F                     # Balanced/Imbalanced subgroup
gr.main <- 2                     # The number of rough subgroup 
gr.GE <- 4                       # The number of refined subgroup 
num_main <- list(c(1,2),c(3,4))  # Hierarchical subgroup structure
sd.epsi<-0.5                     # The residual error
corZX <- F                       # Correlated/Uncorrelated X and Z
intercorr0<-"AR1"                # The correlation structure of Z


# ------------ Necessary parameters to support algorithm implementation --------
sample.index <- combn(n,2)
a.matrix <- kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(pp)))
a.matrix.main <- kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(q)))
a.matrix.GE <- kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(p)))
one.matrix <- bdiag(rep(list(rep(1,pp)),(n*(n-1)/2)))
one.matrix.main <- bdiag(rep(list(rep(1,q)),(n*(n-1)/2)))
one.matrix.GE <- bdiag(rep(list(rep(1,p)),(n*(n-1)/2)))

# ------------ Generate the regression coefficients -----------
beta_main.vec <- c(-beta00,beta00)
beta_GE.vec <- c(-1.5*beta00,-0.5*beta00,1.5*beta00,0.5*beta00)
if(balance){pr_GE<-rep(1/gr.GE,gr.GE)}else{pr_GE<-c(1/6,1/6,1/3,1/3)}
group.num <- group_GE(n, pr_GE=pr_GE, num_main=num_main)
group.num_main <- group.num$group.num_main
group.num_GE <- group.num$group.num_GE
Beta0 <- Generatebeta(n, q, p, group.num_main, group.num_GE, beta_main.vec, beta_GE.vec, sparse=T, S_GE=S_GE, S_E=S_E)
beta.true.list <- Beta0$beta.true.list


# ------------ Generate data ------------ 
whole.data <- GenerateData(n, p, q, beta.true.list, maincorr="AR_E", intercorr=intercorr0,
                        sd.epsi=sd.epsi, corZX = corZX, E.cat = E.cat)

# ------------ Apply the proposed method ------------ 
beta.init.list <- gen_int_beta(n, p, q, whole.data)
beta.init <- beta.init.list$beta.init
lambda <- genelambda.obo()
result <- tuning.lambda(lambda, whole.data, n, q, p, beta.init)
index.list <- evaluation.tuning(n,q,p, result$admmres, result$abic.n, result$admmres2, result$abic.n2, Beta0, gr.GE=gr.GE, result$bic.var)
index.list$err.s

