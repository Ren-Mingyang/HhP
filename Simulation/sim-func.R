################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances,
# which are used to support numerical simulation studies.
################################################################################


################################################################################
# Functions for generation of simulated data
# This section includes four functions: 
# CorMat()    group_GE()    Generatebeta()    GenerateData()
################################################################################
CorMat <- function(m, sig.index = "En"){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: CorMat
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the correlation structure.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ m: Dimensions of the correlation matrix.
  ## @ sig.index: The selection of the correlation structure of the design matrix
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ sig.G: The correlation matrix.
  ## ---------------------------------------------------------------------------------------------------------------
  
  sig.G= matrix(0,m,m)
  if (sig.index != "En"){
    for (i in 1:m){
      for(j in i:m){
        if(sig.index=="AR1"){
          sig.G[i,j]=0.25^abs(i-j)
        }else if(sig.index=="AR2"){
          sig.G[i,j]=0.75^abs(i-j)
        }else if(sig.index=="B2"){
          if(abs(i-j) ==1){
            sig.G[i,j]= 0.6
          }else if (abs(i-j)==2){
            sig.G[i,j] = 0.3
          }
        }else if(sig.index=="B1"){
          if(abs(i-j) ==1){
            sig.G[i,j]= 0.3
          }
        }else if(sig.index=="AR_E"){
          sig.G[i,j]=0.2
        }
        sig.G[j,i] = sig.G[i,j]
      }
    }
  }
  diag(sig.G)= 1
  return(sig.G)
}

group_GE <- function(n, pr_GE=rep(1/4,4), num_main=list(c(1,2),c(3,4)))
{
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: group_GE
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the hierarchical subgroup structure.
  ## ---------------------------------------------------------------------------------------------------------------
  
  group.num_GE=list()
  n.num = sample(n,n,replace = F)
  group.num_GE[[1]] = sort(n.num[1:ceiling(n*pr_GE[1])])
  for (g in 2:(length(pr_GE))) {
    group.num_GE[[g]] = sort(n.num[(ceiling(n*sum(pr_GE[1:(g-1)]))+1):ceiling(n*sum(pr_GE[1:g]))])
  }
  
  group.num_main=list()
  for (g in 1:(length(num_main))) {
    num_maing = unlist(num_main[g])
    num_maing.vec = NULL
    for (gg in 1:length(num_maing)) {
      num_maing.vec = c(num_maing.vec,group.num_GE[[num_maing[gg]]])
    }
    group.num_main[[g]] = sort(num_maing.vec)
  }
  
  return(list(group.num_main=group.num_main,group.num_GE=group.num_GE))
  
}

Generatebeta <-function(n, q, p, group.num_main, group.num_GE, beta_main.vec, 
                        beta_GE.vec, sparse=F, S_E=1, S_GE=1)
{
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: Generatebeta
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the hierarchical regression coefficients.
  ## ---------------------------------------------------------------------------------------------------------------
  
  if(!sparse){S_E=q;S_GE=p}
  
  beta_main = matrix(0,length(group.num_main),S_E)
  for (i in 1:length(group.num_main)) {
    beta_main[i,] = rep(beta_main.vec[i],S_E)
  }
  beta_GE = matrix(0,length(group.num_GE),S_GE)
  for (i in 1:length(group.num_GE)) {
    beta_GE[i,] = rep(beta_GE.vec[i],S_GE)
  }
  
  pp = q+p
  setE <- sort(sample(1:q, size = S_E, replace = FALSE))
  setGE <- sort(sample(1:p, size = S_GE, replace = FALSE))
  
  # main: true beta
  E_n.sig = setE
  beta.true.main.list<-vector(mode = "list",length=n)
  beta.true.main.list0<-vector(mode = "list",length=n)
  for(i in 1:length(group.num_main)){
    for(j in 1:length(group.num_main[[i]])){
      betai = rep(0,pp)
      betai[c(E_n.sig)] = beta_main[i,]
      ab <- as.matrix(betai)
      beta.true.main.list[[group.num_main[[i]][j]]]<-as.matrix(ab)
      beta.true.main.list0[[group.num_main[[i]][j]]]<-as.matrix(ab[c(1:q)])
    }
  }
  
  # interaction: true beta
  GE_n.sig = setGE
  beta.true.GE.list<-vector(mode = "list",length=n)
  beta.true.GE.list0<-vector(mode = "list",length=n)
  for(i in 1:length(group.num_GE)){
    for(j in 1:length(group.num_GE[[i]])){
      betai = rep(0,pp)
      betai[q+GE_n.sig] = beta_GE[i,]
      ab <- as.matrix(betai)
      beta.true.GE.list[[group.num_GE[[i]][j]]]<-as.matrix(ab)
      beta.true.GE.list0[[group.num_GE[[i]][j]]]<-as.matrix(ab[setdiff(c(1:pp),c(1:q))])
    }
  }
  
  beta.true.list<-vector(mode = "list",length=n)
  for (g in 1:n) {
    beta.true.list[[g]] = beta.true.main.list[[g]] + beta.true.GE.list[[g]]
  }
  beta.true.mat<-t(do.call(cbind,beta.true.list))
  
  return(list(beta.true.list=beta.true.list,
              beta.true.main.list0=beta.true.main.list0,
              beta.true.GE.list0=beta.true.GE.list0,
              beta.true.mat=beta.true.mat,
              setGE=setGE,setE=setE))
}

GenerateData <- function(n, p, q, beta.true.list, maincorr="En", intercorr="En", 
                         corZX = F, sd.epsi=0.5, p_y=1, alpha0=0, E.cat = 0)
{
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: GenerateData
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the simulated data.
  ## ---------------------------------------------------------------------------------------------------------------
  
  Zsigmat <- CorMat(p, intercorr)
  Xsigmat <- CorMat(q, maincorr)
  sigmat <- as.matrix(Matrix::bdiag(Xsigmat,Zsigmat))
  
  if(corZX){
    sigmat[sigmat==0] <- min(abs(sigmat[sigmat>0]))/20
  }
  
  pp <- p+q
  XZ <- MASS::mvrnorm(n, rep(0,pp), sigmat)
  X <- XZ[,1:q]
  Z <- XZ[,(q+1):pp]
  
  if (E.cat > 0){
    X[, 1:E.cat] <- as.numeric(X[,1:E.cat] > 0 )
    Z[, 1:(E.cat*p/q)] <- as.numeric(Z[,1:(E.cat*p/q)] > 0 )
  }
  
  cova <- as.matrix(cbind(X,Z))
  
  # error
  err.sigma<-diag(p_y)*sd.epsi
  err<-MASS::mvrnorm(n, rep(0,p_y), err.sigma)
  
  # generate response
  respon<-matrix(0,nrow=n,ncol=p_y)
  for(i4 in 1:n){
    respon[i4,]<-cova[i4,]%*%beta.true.list[[i4]]+err[i4,]+alpha0
  }
  
  
  return(list(data_y=respon,data_x=cova,beta.true.list=beta.true.list,data_X=X,data_Z=Z,data_e=err))
  
}


################################################################################
# Functions for evaluating performances of proposed methods
################################################################################
evaluation <- function(n, q, p, beta.gfl,eta.gfl,Beta0)
{
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: evaluation
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Evaluating performances of proposed methods.
  ## -----------------------------------------------------------------------------------------------------------------
  
  pp=q+p
  main.n=c(1:q)
  GE.n=setdiff(1:pp,main.n)
  a.main.n = sort(rep((0:(n-1))*pp,q))+main.n
  a.GE.n = sort(rep((0:(n-1))*pp,p))+GE.n
  d.main.n =sort(rep((0:(n*(n-1)/2-1))*pp,q)) + main.n
  d.GE.n = sort(rep((0:(n*(n-1)/2-1))*pp,p)) + GE.n
  
  b.main = as.matrix(beta.gfl[a.main.n,])
  b.GE = as.matrix(beta.gfl[a.GE.n,])
  eta.main = as.matrix(eta.gfl[d.main.n,],sparse=T)
  eta.GE = as.matrix(eta.gfl[d.GE.n,],sparse=T)
  
  beta.true.list = Beta0$beta.true.list
  beta.true.main.list = Beta0$beta.true.main.list0
  beta.true.GE.list = Beta0$beta.true.GE.list0
  
  beta.main.esi<-t(matrix(b.main,nrow = q))
  beta.main.true<-t(do.call(cbind,beta.true.main.list))
  error.main <- sqrt(sum((beta.main.esi - beta.main.true)^2)) / sqrt(n*q)
  beta.GE.esi<-t(matrix(b.GE,nrow = p))
  beta.GE.true<-t(do.call(cbind,beta.true.GE.list))
  error.GE <- sqrt(sum((beta.GE.esi - beta.GE.true)^2)) / sqrt(n*p)
  
  RI.main<-getErrorRate(n, b.main,eta.main, beta.true.main.list, q, a.matrix.main)
  gr.main.est<-RI.main$gr.num
  RI.GE<-getErrorRate(n, b.GE,eta.GE, beta.true.GE.list, p, a.matrix.GE)
  gr.GE.est<-RI.GE$gr.num
  
  return(list(RI.main=RI.main, RI.GE=RI.GE, error.main=error.main, error.GE=error.GE))
}

evaluation.tuning <- function(n,q,p,admmres, abic.n, admmres2, abic.n2, 
                              Beta0, gr.GE=4, bic.var)
{
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: evaluation.tuning
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##    Evaluating performances of proposed methods under all tuning parameters,
  ##    and selecting the optional results.
  ## -----------------------------------------------------------------------------------------------------------------
  
  L1 <- length(admmres)
  L2 <- length(admmres2)
  err.all <- matrix(0,L1,6)
  for (l1 in 1:L1) {
    ad <- admmres[[l1]]
    reseva <- evaluation(n,q, p, ad$beta.gfl,ad$eta.gfl,Beta0)
    RI.main <- reseva$RI.main
    RI.GE <- reseva$RI.GE
    err.all[l1,] <- c(RI.main$gr.num,1-RI.main$group.error[1],reseva$error.main,RI.GE$gr.num,1-RI.GE$group.error[1],reseva$error.GE)
  }
  
  err.all2 <- matrix(0,L2,6)
  for (l2 in 1:L2) {
    ad <- admmres2[[l2]]
    reseva <- evaluation(n,q, p, ad$beta.gfl,ad$eta.gfl,Beta0)
    RI.main <- reseva$RI.main
    RI.GE <- reseva$RI.GE
    err.all2[l2,] <- c(RI.main$gr.num,1-RI.main$group.error[1],reseva$error.main,RI.GE$gr.num,1-RI.GE$group.error[1],reseva$error.GE)
  }
  
  n4 <- which(err.all[,4] == gr.GE)
  if(length(n4) > 0){
    bic.part <- bic.var[n4]
    err.part <- err.all[n4,]
    if(length(n4) == 1){err.part <- t(as.matrix(err.all[n4,]))}
    abic.part <- which(bic.part == min(bic.part[!is.na(bic.part)]))[1]
    err.s.or <- err.part[abic.part,]
  }
  if(length(n4) == 0){
    err.s.or <- err.all[abic.n,]
    err.s.or[4:6] <- 0
  }
  
  gr.main.n <- err.all[abic.n,1]
  n.big <- which(err.all[,4] > gr.main.n)
  if(length(n.big) > 0){
    bic.part <- bic.var[n.big]
    err.part <- err.all[n.big,]
    if(length(n.big) == 1){err.part <- t(as.matrix(err.all[n.big,]))}
    abic.part <- which(bic.part == min(bic.part[!is.na(bic.part)]))[1]
    err.s.big <- err.part[abic.part,]
  }
  if(length(n.big) == 0){
    err.s.big <- err.all[abic.n,]
    err.s.big[4:6] <- 0
  }
  
  result = list(err.s.or=err.s.or, err.s=err.s.big)
  
  return(result)
  
}

getErrorRate <- function(size, beta.gfl, eta.gfl, beta.true.list.aa, q, a.matrix.aa)
{
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: getErrorRate
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the indicators of grouping accuracy.
  ## -----------------------------------------------------------------------------------------------------------------
  
  # true group matrix
  beta.true.mat<-do.call(rbind,beta.true.list.aa)
  diffbeta<-a.matrix.aa%*%as.matrix(beta.true.mat)
  dim(diffbeta)<-c(q,size*(size-1)/2)
  epsi.betatrue<-apply(diffbeta,2,sum)
  
  captrue.matrix<-matrix(0,nrow=size,ncol=size)
  for(i in 1:length(epsi.betatrue)){
    if(epsi.betatrue[i]==0){
      captrue.matrix[sample.index[1,i],sample.index[2,i]]<-1
    }
  }
  
  diff.gfl<-apply(matrix(eta.gfl,nrow=q,ncol=size*(size-1)/2),2,sum)
  
  capgfl.matrix<-matrix(0,nrow=size,ncol=size)
  if(length(which(diff.gfl==0))==0){
    captrue.matrix2 = captrue.matrix+t(captrue.matrix); diag(captrue.matrix2)=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group_error_rate<-sum(abs(captrue.matrix-capgfl.matrix))/(size*(size-1)/2)
    group_tp<-length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==1)))/sum(captrue.matrix)
    group_fp<-length(intersect(which(as.vector(captrue.matrix)==0),which(as.vector(capgfl.matrix)==1)))/(size*(size-1)/2-sum(captrue.matrix))
    return(list(group.error=c(group_error_rate,group_tp,group_fp), gr.num=size, matrix=list(captrue.matrix2=captrue.matrix2,capgfl.matrix2=capgfl.matrix2)))
  }
  if(length(which(diff.gfl==0))==size*(size-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)]=1
    captrue.matrix2 = captrue.matrix+t(captrue.matrix); diag(captrue.matrix2)=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group_error_rate<-sum(abs(captrue.matrix-capgfl.matrix))/(size*(size-1)/2)
    group_tp<-length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==1)))/sum(captrue.matrix)
    group_fp<-length(intersect(which(as.vector(captrue.matrix)==0),which(as.vector(capgfl.matrix)==1)))/(size*(size-1)/2-sum(captrue.matrix))
    return(list(group.error=c(group_error_rate,group_tp,group_fp), gr.num=1, matrix=list(captrue.matrix2=captrue.matrix2,capgfl.matrix2=capgfl.matrix2)))
  }else{
    sample.index.gfl<-sample.index[,which(diff.gfl==0)]
    if(length(which(diff.gfl==0))==1){
      capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]]<-1
    }else{
      for(i in 1:length(which(diff.gfl==0))){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]]<-1
      }
    }
    
    captrue.matrix2 = captrue.matrix+t(captrue.matrix); diag(captrue.matrix2)=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(capgfl.matrix2))
    
    group_error_rate<-sum(abs(captrue.matrix-capgfl.matrix))/(size*(size-1)/2)
    group_tp<-length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==1)))/sum(captrue.matrix)
    group_fp<-length(intersect(which(as.vector(captrue.matrix)==0),which(as.vector(capgfl.matrix)==1)))/(size*(size-1)/2-sum(captrue.matrix))
    group_fn<-length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==0)))/sum(captrue.matrix)
    return(list(group.error=c(group_error_rate,group_tp,group_fp), gr.num=group.num.gf, matrix=list(captrue.matrix2=captrue.matrix2,capgfl.matrix2=capgfl.matrix2)))
  }
}

