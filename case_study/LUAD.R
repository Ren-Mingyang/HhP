################################################################################
# Codes for hierarchical heterogeneity analysis of LUAD image data
################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(Matrix)
library(MASS)
library(fmrs)
source("function.R")

# ------------- load data ----------------
data <- read.csv("LUAD-data.csv")
data_y <- data$y
conyi_names <- setdiff(names(data),c("SEX","TumorSN"))
data[,conyi_names] <- scale(data[,conyi_names])
data_x <- data[,-1]
data_y <- data[,1]
whole.data <- list()
whole.data$data_y <- as.matrix(data_y)
whole.data$data_x <- as.matrix(data_x)


# ------------- Set the necessary parameters ----------------
n <- dim(data_x)[1]
q <- 6
p <- dim(data_x)[2] - q
pp <- p+q
sample.index<-combn(n,2)
a.matrix<-kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(pp)))
a.matrix.main<-kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(q)))
a.matrix.GE<-kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(p)))
one.matrix<-bdiag(rep(list(rep(1,pp)),(n*(n-1)/2)))
one.matrix.main<-bdiag(rep(list(rep(1,q)),(n*(n-1)/2)))
one.matrix.GE<-bdiag(rep(list(rep(1,p)),(n*(n-1)/2)))

# ------------- Apply the proposed method ----------------
beta.init.list <- gen_int_beta(n, p, q, whole.data, ridge=T, gr.init=5, lambda.min=0.001)
beta.init.list$nC
beta.init.mat <- t(matrix(beta.init.list$beta.init,nrow = p+q))

ptm<-proc.time()
lambda <- genelambda.obo(nlambda1=6,lambda1_max=1.2,lambda1_min=0.8,nlambda2=4,lambda2_max=1.5,lambda2_min=1.2)
result <- tuning.lambda(lambda, whole.data, n, q, p, beta.init.list$beta.init, selection.sub=T)
t1 <- proc.time() - ptm

admm.sum <- result$admmres
num.sum <- matrix(0,length(admm.sum),2)
BIC <- c()
for (j in 1:length(admm.sum)){
  admmj <- admm.sum[[j]]
  BIC[j] <- admmj$BIC.var
  num.sum[j,1] <- admmj$num.main$gr.num
  num.sum[j,2] <- admmj$num.GE$gr.num
}
select.num <- which(num.sum[,1] > 1 & num.sum[,2] > num.sum[,1])
opt_n <- which(BIC == min(BIC[select.num]))
opt_admm <- admm.sum[[opt_n]]
cap_main <- opt_admm$num.main$capgfl.matrix2
cap_GE <- opt_admm$num.GE$capgfl.matrix2
sum(cap_main[cap_GE == 1] == 0)
b.main <- opt_admm$b.main
b.GE <- opt_admm$b.GE
b.main.mat <- t(matrix(b.main,ncol = n))
b.GE.mat <- t(matrix(b.GE,ncol = n))

# The type 1 features
num_subgroup <- unique(apply(cap_main, 1, function(a){which(a == 1)}))
K1 <- length(num_subgroup)
memb <- rep(0,n)
coe.main <- matrix(0,q,K1)
for (k in 1:K1) {
  memb[num_subgroup[[k]]] <- k
  coe.main[,k] <- b.main.mat[num_subgroup[[k]][1],]
}
coe.main <- as.data.frame(coe.main)
row.names(coe.main) <- names(data_x)[1:q]
names(coe.main) <- paste("Subgroup",c(1:K1))
memb_main <- memb
num_subgroup_main <- num_subgroup

# The type 2 features
num_subgroup <- unique(apply(cap_GE, 1, function(a){which(a == 1)}))
K2 <- length(num_subgroup)
memb <- rep(0,n)
coe.GE <- matrix(0,p,K2)
for (k in 1:K2) {
  memb[num_subgroup[[k]]] <- k
  coe.GE[,k] <- b.GE.mat[num_subgroup[[k]][1],]
}
coe.GE <- as.data.frame(coe.GE)
row.names(coe.GE) <- names(data_x)[(q+1):(q+p)]
names(coe.GE) <- paste("Subgroup",c(1:K2))
memb_GE <- memb
num_subgroup_GE <- num_subgroup

# Subgroup size
s1 <- c()
for (k in 1:length(num_subgroup_main)) {
  s1[k] <- length(num_subgroup_main[[k]])
}
s1
s2 <- c()
for (k in 1:length(num_subgroup_GE)) {
  s2[k] <- length(num_subgroup_GE[[k]])
}
s2

# Output the results
GE.coe <- rbind(coe.GE,s2)
GE.coe <- GE.coe[,]
write.csv(rbind(coe.main,s1),"LUAD.coe.main.csv")
write.csv(GE.coe,"LUAD.coe.GE.csv")

