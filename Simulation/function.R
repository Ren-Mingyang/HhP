#####################################################################################
# This document includes main functions of proposed methods, 
# which are used to support numerical simulation studies and real data analysis
# in the paper
# "Hierarchical Cancer Heterogeneity Analysis based on Histopathological Imaging Data"
#####################################################################################


############################# Functions for main algorithms ############################
ADMM<-function(n, q, p, data.total, beta.init, lambda1, lambda2, p_y=1,
               iter.max=50, epsi=0.1, gam.mcp=3, penal.para=1, merge.all=F)
{
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            The key function of hierarchical heterogeneity analysis:
  ##            The implementation of ADMM.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ n: The sample size.
  ## @ q: The dimension of type 1 features.
  ## @ p: The dimension of type 2 features.
  ## @ data.total: The input data analyzed (a list including the response and design matrix).
  ## @ beta.init: The Initial values of regression coefficients.
  ## @ lambda1: The tuning parameter controlling the number of refined subgroup.
  ## @ lambda2: the tuning parameter controlling the number of rough subgroup.
  ## @ p_y: The dimension of y.
  ## @ iter.max: int, Maximum number of cycles of the ADMM algorithm, the default setting is 5.
  ## @ epsi: a float value, algorithm termination threshold.
  ## @ gam.mcp: The regularization parameter in MCP, the default setting is 3.
  ## @ penal.para: The penalty parameter in ADMM algorithm, the default setting is 1.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including estimated regression coefficients, subgroups, BIC value, and so on 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  pp=q+p
  main.n=c(1:q)
  GE.n=setdiff(1:pp,main.n)
  a.main.n = sort(rep((0:(n-1))*pp,q))+main.n
  a.GE.n = sort(rep((0:(n-1))*pp,p))+GE.n
  d.main.n =sort(rep((0:(n*(n-1)/2-1))*pp,q)) + main.n
  d.GE.n = sort(rep((0:(n*(n-1)/2-1))*pp,p)) + GE.n
  
  cova<-data.total$data_x
  data.response<-data.total$data_y
  
  row_vec<-vector(mode="list",length=n)
  for(i5 in 1:n){
    row_vec[[i5]]<-matrix(cova[i5,],1,pp)
  }
  cova_diag<-bdiag(row_vec)
  
  #################################### Initialize beta #########################################
  # group
  eta.init<-a.matrix%*%beta.init
  gamma.init<-Matrix(0,nrow=p_y,ncol=pp*n*(n-1)/2,sparse = T)
  
  # iteration
  iter<-1
  # coe
  beta.est.list<-vector(mode="list",length=iter.max);beta.est.list[[iter]]<-beta.init
  # differences of coe
  eta.est.list<-vector(mode="list",length=iter.max);eta.est.list[[iter]]<-eta.init
  # the dual variables
  gamma.est.list<-vector(mode="list",length=iter.max);gamma.est.list[[iter]]<-gamma.init
  
  while(iter<=iter.max){
    iter<-iter+1
    if(iter>iter.max){
      break
    }
    
    beta.est.list[[iter]]<-solve(t(cova_diag)%*%cova_diag+penal.para*t(a.matrix)%*%a.matrix)%*%
      (t(cova_diag)%*%data.response-t(a.matrix)%*%t(gamma.est.list[[iter-1]])+penal.para*(t(a.matrix)%*%eta.est.list[[iter-1]]))
    
    eta.i = a.matrix%*%beta.est.list[[iter]]
    gamma.i = gamma.est.list[[iter-1]]
    eta.tem1<-eta.i+t(gamma.i)/penal.para
    eta.tem1.main<-as.matrix(eta.tem1[d.main.n,],sparse=T)
    eta.tem1.GE<-as.matrix(eta.tem1[d.GE.n,],sparse=T)
    eta.i1<-eta.i
    
    # 2-norm
    eta.tem2norm.main<-sqrt(t(one.matrix.main)%*%(eta.tem1.main^2))
    eta.tem2norm.GE<-sqrt(t(one.matrix.GE)%*%(eta.tem1.GE^2))
    eta.tem2norm<-sqrt(eta.tem2norm.GE^2 + eta.tem2norm.main^2)
    
    # all:non-shrinkage & main:non-shrinkage
    num_n_n<-which(eta.tem2norm > gam.mcp*lambda1 & eta.tem2norm.main > gam.mcp*lambda2)
    
    # all:shrinkage & main:non-shrinkage
    vu.coe<-apply(1-lambda1/penal.para/eta.tem2norm,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
    vv<- vu.coe * eta.tem2norm.main
    num_s_n<-which(eta.tem2norm <= gam.mcp*lambda1 & vv > gam.mcp*lambda2)
    
    # all:non-shrinkage & main:shrinkage
    vv.coe<-apply(1-lambda2/penal.para/eta.tem2norm.main,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
    ww<-sqrt(eta.tem2norm.GE^2 + ( vv.coe * eta.tem2norm.main)^2)
    num_n_s<-which(ww > gam.mcp*lambda1 & eta.tem2norm.main <= gam.mcp*lambda2)
    
    # all:shrinkage & main:shrinkage
    num_s_s<-setdiff(1:(n*(n-1)/2),Reduce(union,list(num_n_n,num_s_n,num_n_s)))
    
    if(length(num_n_n)<2){
      num_n_n2<-which(rowSums(as.matrix(one.matrix[,num_n_n]))!=0)
    }else{num_n_n2<-which(rowSums(one.matrix[,num_n_n])!=0)}
    if(length(num_s_n)<2){
      num_s_n2<-which(rowSums(as.matrix(one.matrix[,num_s_n]))!=0)
    }else{num_s_n2<-which(rowSums(one.matrix[,num_s_n])!=0)}
    if(length(num_n_s)<2){
      num_n_s2<-which(rowSums(as.matrix(one.matrix[,num_n_s]))!=0)
    }else{num_n_s2<-which(rowSums(one.matrix[,num_n_s])!=0)}
    if(length(num_s_s)<2){
      num_s_s2<-which(rowSums(as.matrix(one.matrix[,num_s_s]))!=0)
    }else{num_s_s2<-which(rowSums(one.matrix[,num_s_s])!=0)}
    
    # all:non-shrinkage & main:non-shrinkage
    if(length(num_n_n2) > 0){eta.i1[num_n_n2,]<-eta.i[num_n_n2,]}
    length(num_n_n2)
    # all:shrinkage & main:non-shrinkage
    num = num_s_n; num2 = num_s_n2;
    if(length(num) > 0){
      eta.tem3<-as.matrix(vu.coe[num])
      eta.tem4<-as.vector(apply(eta.tem3, 1, function(x) rep(x,pp)))
      eta.i1[num2,]<-eta.tem4*eta.i[num2,]
    }
    
    # all:non-shrinkage & main:shrinkage
    num = num_n_s; num2 = num_n_s2;
    if(length(num) > 0){
      num.main = num2[sort(rep((1:length(num)-1)*pp,q)) + main.n]
      num.GE = num2[sort(rep((1:length(num)-1)*pp,p)) + GE.n]
      eta.tem3<-as.matrix(vv.coe[num])
      eta.tem4<-as.vector(apply(eta.tem3, 1, function(x) rep(x,q)))
      eta.i1[num.main,]<-eta.tem4*eta.i[num.main,]
      eta.i1[num.GE,]<-eta.i[num.GE,]
    }
    
    # all:shrinkage & main:shrinkage
    num = num_s_s; num2 = num_s_s2;
    if(length(num) > 0){
      num.main = num2[sort(rep((1:length(num)-1)*pp,q)) + main.n]
      num.GE = num2[sort(rep((1:length(num)-1)*pp,p)) + GE.n]
      eta.i0norm<-sqrt(t(one.matrix)%*%(eta.est.list[[iter-1]]^2))
      eta.i0norm.main<-sqrt(t(one.matrix.main)%*%(eta.est.list[[iter-1]][d.main.n,]^2))
      mcp.u<-mcp_d(eta.i0norm[num],lambda1,gam.mcp)/penal.para/(eta.i0norm[num]+10^(-7))
      mcp.v<-mcp_d(eta.i0norm.main[num],lambda2,gam.mcp)/penal.para/(eta.i0norm.main[num]+10^(-7))
      eta.tem3.main<-as.matrix(mcp.u+mcp.v)
      eta.tem3.GE<-as.matrix(mcp.u)
      eta.tem4.main<-as.vector(apply(eta.tem3.main, 1, function(x) rep(x,q)))
      eta.tem4.GE<-as.vector(apply(eta.tem3.GE, 1, function(x) rep(x,p)))
      eta.i1[num.main,]<-eta.i[num.main,]/(1+eta.tem4.main)
      eta.i1[num.GE,]<-eta.i[num.GE,]/(1+eta.tem4.GE)
      eta.i1[abs(eta.i1) < 10^(-9)]=0
    }
    
    eta.est.list[[iter]] = eta.i1
    gamma.est.list[[iter]]<-gamma.est.list[[iter-1]]+penal.para*t(a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])
    eps.group = sqrt(sum((a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])^2))
    
    if(eps.group<epsi){
      break
    }
    
    beta.est.list[iter-1]<-list(NULL)
    eta.est.list[iter-1]<-list(NULL)
    gamma.est.list[iter-1]<-list(NULL)
  }
  
  if(iter>iter.max){
    beta.gfl<-beta.est.list[[iter-1]]
    eta.gfl<-eta.est.list[[iter-1]]
  }else{
    beta.gfl<-beta.est.list[[iter]]
    eta.gfl<-eta.est.list[[iter]]
  }
  
  b.main = as.matrix(beta.gfl[a.main.n,])
  b.GE = as.matrix(beta.gfl[a.GE.n,])
  eta.main = as.matrix(eta.gfl[d.main.n,],sparse=T)
  eta.GE = as.matrix(eta.gfl[d.GE.n,],sparse=T)
  
  rm(beta.est.list)
  rm(eta.est.list)
  rm(gamma.est.list)
  gc()
  ###################
  num.main <- get.gr.num(n, b.main,eta.main, q, merge.all=merge.all)
  num.GE <- get.gr.num(n, b.GE,eta.GE, p, merge.all=merge.all)
  gr.main.est<-num.main$gr.num
  gr.GE.est<-num.GE$gr.num
  residual <- log(sum((data.response-cova_diag%*%beta.gfl)^2/n))
  BIC.var<-residual+log(n*pp)*log(n)*(gr.main.est*(q)+gr.GE.est*(p))/n
  BIC.o<-residual+log(n)*(gr.main.est*q+gr.GE.est*p)/n
  
  return(list(beta.gfl=beta.gfl, b.main=b.main, b.GE=b.GE, eta.gfl=eta.gfl, eta.main=eta.main, eta.GE=eta.GE,
              residual=residual,BIC.var=BIC.var,BIC.o=BIC.o,num.main=num.main,num.GE=num.GE))
}

############################## Functions for tuning lambdas ############################
genelambda.obo = function(nlambda1=20,lambda1_max=0.5,lambda1_min=0.1,
                          nlambda2=5,lambda2_max=1.5,lambda2_min=0.1)
{
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: genelambda.obo
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating a sequence of the tuning parameters (lambda1 and lambda2).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ nlambda1, nlambda2: The numbers of lambda 1 2.
  ## @ lambda1_min, lambda2_min: The minimum values of lambda 1 2.
  ## @ lambda1_max, lambda2_max: The maximum values of lambda 1 2.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ lambda: a sequence of the tuning parameters (lambda1 and lambda2).   
  ## -----------------------------------------------------------------------------------------------------------------
  
  lambda1 = exp(seq(log(lambda1_max),log(lambda1_min),len= nlambda1))
  lambda2 =exp(seq(log(lambda2_max),log(lambda2_min),len= nlambda2))
  lambda = list(lambda1=lambda1,lambda2=lambda2)
  return(lambda)
}

tuning.lambda <- function(lambda, whole.data, n, q, p, beta.init, 
                          merge.all=F, trace=F, selection.sub=F)
{
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Searching and selecting the optional tuning parameters under
  ##            the adaptive BIC-type criterion using the proposed method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ lambda: The sequences of the tuning parameters (lambda1 and lambda2).
  ## @ whole.data: The input data analyzed (a list including the response and design matrix).
  ## @ n: The sample size.
  ## @ q: The dimension of type 1 features.
  ## @ p: The dimension of type 2 features.
  ## @ beta.init: The Initial values of regression coefficients.
  ## @ trace: the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including results corresponding all choices of given tuning parameters. 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  L1 = length(lambda1)
  L2 = length(lambda2)
  L = L1+L2
  
  lam2 = lambda2
  lam1 = lam1=rep(0,length(lam2))
  admmres2=vector(mode="list",length=length(lam2))
  bic.var2<-c()
  gr.nn2 <- c()
  for(l2 in 1:L2){
    if(trace){if(l2%%5==1){cat('-----------',l2,'-th lambda--------------\n')}}
    ad<-ADMM(n, q, p, whole.data, beta.init, lam1[l2], lam2[l2], merge.all=merge.all)
    bic.var2[l2]<-ad$BIC.var
    RI.main=ad$num.main
    RI.GE=ad$num.GE
    admmres2[[l2]] <- ad
    gr.nn2[l2] <- RI.main$gr.num
    if(trace){
      print(c(RI.main$gr.num,RI.GE$gr.num))
    }
  }
  
  if(sum(gr.nn2 > 1) > 0){
    abic.n2 = which(bic.var2 == min(bic.var2[(!is.na(bic.var2)) & (gr.nn2 > 1)]))[1]
  } else {
    abic.n2 = which(bic.var2 == min(bic.var2[(!is.na(bic.var2)) & (gr.nn2 > 0)]))[1]
  }
  
  lam20=lam2[abic.n2][1]
  lam1=lambda1
  lam2=rep(lam20,length(lam1))
  admmres=vector(mode="list",length=length(lam1))
  bic.var<-c()
  for(l1 in 1:L1){
    if(trace){if((L2+l1)%%5==1){cat('-----------',L2+l1,'-th lambda--------------\n')}}
    ad<-ADMM(n, q, p, whole.data, beta.init, lam1[l1], lam2[l1], merge.all=merge.all)
    bic.var[l1]<-ad$BIC.var
    RI.main=ad$num.main
    RI.GE=ad$num.GE
    admmres[[l1]] <- ad
    if(trace){
      print(c(RI.main$gr.num,RI.GE$gr.num))
    }
  }
  abic.n = which(bic.var == min(bic.var[!is.na(bic.var)]))[1]
  beta.over<-admmres[[abic.n]]$beta.gfl
  RI.main <- admmres[[abic.n]]$RI.main
  RI.GE <- admmres[[abic.n]]$RI.GE
  result = list(beta.over=beta.over,bic.var=bic.var,abic.n=abic.n, admmres=admmres,
                bic.var2=bic.var2, abic.n2=abic.n2, admmres2=admmres2)
  
  if(selection.sub){
    admm.sum <- admmres
    num.sum <- matrix(0,length(admm.sum),2)
    BIC <- c()
    for (j in 1:length(admm.sum)) {
      admmj <- admm.sum[[j]]
      BIC[j] <- admmj$BIC.var
      num.sum[j,1] <- admmj$num.main$gr.num
      num.sum[j,2] <- admmj$num.GE$gr.num
    }
    select.num <- which(num.sum[,1] > 1 & num.sum[,2] > num.sum[,1])
    if(length(select.num) == 0){
      result = list(beta.over=beta.over,bic.var=bic.var,abic.n=abic.n, admmres=admmres,
                    bic.var2=bic.var2, abic.n2=abic.n2, admmres2=admmres2)
    }
    if(length(select.num) > 0){
      abic.n <- which(BIC == min(BIC[select.num]))
      beta.over<-admmres[[abic.n]]$beta.gfl
      RI.main <- admmres[[abic.n]]$RI.main
      RI.GE <- admmres[[abic.n]]$RI.GE
      result = list(beta.over=beta.over,bic.var=bic.var,abic.n=abic.n, admmres=admmres,
                    bic.var2=bic.var2, abic.n2=abic.n2, admmres2=admmres2)
    }
  }
  return(result)
}

############################# Some fundamental supporting functions ############################
dMatrixFun <- function(indx){
  e.vec<-matrix(0,n,1)
  e.vec[indx[1],]<-1
  e.vec[indx[2],]<-(-1)
  return(e.vec)
}

dMatrixFunj <- function(indx){
  e.vec<-matrix(0,nj,1)
  e.vec[indx[1],]<-1
  e.vec[indx[2],]<-(-1)
  return(e.vec)
}

mcp_d <- function(x,lambda,a=3){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: mcp_d
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the derivative of the MCP
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ x: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ rho: the derivative of the MCP.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  if(lambda!=0){
    rho <- lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}

get.gr.num<-function(size, beta.gfl,eta.gfl, q, merge.all=F){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: get.gr.num
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##          compute the number of estimated subgroups.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input: It can refer to the output of the algorithm
  ## -----------------------------------------------------------------------------------------------------------------
  
  diff.gfl<-apply(matrix(eta.gfl,nrow=q,ncol=size*(size-1)/2),2,sum)
  
  capgfl.matrix<-matrix(0,nrow=size,ncol=size)
  if(length(which(diff.gfl==0))==0){
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    return(list(gr.num=size, capgfl.matrix2=capgfl.matrix2))
  }
  if(length(which(diff.gfl==0))==size*(size-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)]=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    return(list(gr.num=1, capgfl.matrix2=capgfl.matrix2))
  }else{
    sample.index.gfl<-sample.index[,which(diff.gfl==0)]
    if(length(which(diff.gfl==0))==1){
      capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]]<-1
    }else{
      for(i in 1:length(which(diff.gfl==0))){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]]<-1
      }
    }
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(capgfl.matrix2))
    
    if(merge.all){
      cap <- capgfl.matrix2
      num_subgroup <- unique(apply(cap, 1, function(a){which(a == 1)}))
      non_inter_list <- list()
      vv <- 1
      non_inter <- c(1:length(num_subgroup))
      repeat{
        a <- num_subgroup[[non_inter[1]]]
        KK_k <- setdiff(non_inter,non_inter[1])
        non_inter <- c()
        i=1
        for (k2 in KK_k) {
          if(length(intersect(a,num_subgroup[[k2]])) > 0){
            a <- union(a,num_subgroup[[k2]])
          } else {
            non_inter[i] <- k2
            i=i+1
          }
        }
        non_inter_list[[vv]] <- a
        vv <- vv+1
        if(length(non_inter) == 0){break}
      }
      
      for (i in 1:dim(cap)[1]) {
        for (k in 1:length(non_inter_list)) {
          if(length(match(cap[i,],non_inter_list[[k]])) > 0){
            cap[i,non_inter_list[[k]]] <- 1
          }
        }
      }
      capgfl.matrix2 <- cap
      group.num.gf = nrow(unique(capgfl.matrix2))
    }
    
    return(list(gr.num=group.num.gf,capgfl.matrix2=capgfl.matrix2))
  }
  
}

############################# Functions for generating the initial values ############################
gen_int_beta <- function(n, p, q, whole.data, subgroup=c(2,4), 
                         ridge = F, gr.init=5, lambda.min=0.0001){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: gen_int_beta
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the initial values using FMR or ridge fusion method.
  ## -----------------------------------------------------------------------------------------------------------------
  
  pp<-p+q
  n<-length(whole.data$data_y)
  if(!ridge){
    fmr.bic<-c()
    fmr.res0<-c()
    res.mle0.list<-list()
    nC.vec<-subgroup
    for (ii in 1:length(nC.vec)) {
      nC<-nC.vec[ii]
      b00 = c(rep(-1,(nC*pp+nC)/2),rep(1,nC*pp+nC-floor((nC*pp+nC)/2)))
      res.mle <- fmrs.mle(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
                          nComp =nC, disFamily = "norm",
                          # initCoeff = rnorm(nC*pp+nC),
                          initCoeff = b00,
                          initDispersion = rep(1,nC),
                          initmixProp = rep(1/nC, nC),nIterNR = 200)
      res.lam <- fmrs.tunsel(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
                             nComp =nC, disFamily = "norm",
                             initCoeff = c(coefficients(res.mle)),
                             initDispersion = dispersion(res.mle),
                             initmixProp = mixProp(res.mle),
                             penFamily = "adplasso",nIterNR = 200)
      res.var <- fmrs.varsel(y = whole.data$data_y, x =whole.data$data_x, delta = rep(0,length(whole.data$data_y)),
                             nComp = ncomp(res.mle), disFamily = "norm",
                             initCoeff=c(coefficients(res.mle)),
                             initDispersion = dispersion(res.mle),
                             initmixProp = mixProp(res.mle),
                             penFamily = "adplasso",
                             lambPen = slot(res.lam, "lambPen"),nIterNR = 200)
      res.mle0<-res.var
      res.mle0.list[[ii]]<-res.mle0
      
      fmr.res=apply(abs(res.mle0@residuals),1,which.min)
      fmr.resj=c()
      for (jj in 1:n) {
        fmr.resj[jj]=res.mle0@residuals[jj,fmr.res[jj]]
      }
      fmr.res0[ii]<-sqrt(sum(fmr.resj^2))
      fmr.bic[ii]<-res.mle0@BIC
    }
    fmr.coe<-res.mle0@coefficients[-1,]
    fmr.w<-apply(res.mle0@weights,1,which.max)
    fmr.r<-apply(abs(res.mle0@residuals),1,which.min)
    
    # weight
    res.mle0 <- res.mle0.list[[which.max(fmr.bic)]]
    nC<-nC.vec[which.max(fmr.bic)]
    beta.int.fmr<-list()
    for (i in 1:n) {
      beta.int.fmr[[i]]<-as.matrix(as.numeric(fmr.coe[,fmr.w[i]]))
    }
    beta.init=do.call(rbind,beta.int.fmr)
  }
  
  if(ridge){
    ridge.list <- InitialBeta(whole.data, pp, n, lambda.min=lambda.min,gr.init=gr.init)
    beta.init <- ridge.list$beta.init
    beta.init.mat <- t(matrix(beta.init,nrow = pp))
    nC <- length(unique(apply(beta.init.mat, 1, function(a){floor(sum(a))})))
  }
  return(list(beta.init=beta.init,nC=nC))
}

InitialBeta <- function(data.total, q, size, lambda.min=0.001, gr.init=10){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: InitialBeta
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the initial values using the ridge fusion method.
  ## -----------------------------------------------------------------------------------------------------------------
  
  cova<-data.total$data_x
  data.response<-data.total$data_y
  
  # equal beta
  row_vec<-vector(mode="list",length=size)
  for(i5 in 1:size){
    row_vec[[i5]]<-matrix(cova[i5,],1,q)
  }
  cova_diag<-bdiag(row_vec)
  beta.ridge<-solve(t(cova_diag)%*%cova_diag+lambda.min*(t(a.matrix)%*%a.matrix))%*%t(cova_diag)%*%data.response
  
  beta.ridge.mat<-matrix(beta.ridge,nrow=q,ncol=size)
  beta.ridge.med<-apply(beta.ridge.mat,2,median)
  
  group.size<-size/gr.init
  
  # lasso
  beta.order<-order(beta.ridge.med)
  cova_group<-vector(mode="list",length=gr.init)
  respon_group<-vector(mode="list",length=gr.init)
  beta.init.list<-vector(mode="list",length=size)
  for(i in 1:gr.init){
    ind<-beta.order[c(((i-1)*group.size+1):(i*group.size))]
    cova_group[[i]]<-cova[ind,]
    respon_group[[i]]<-data.response[ind,]
    beta.ols<-unname(solve(t(cova_group[[i]])%*%cova_group[[i]])%*%t(cova_group[[i]])%*%respon_group[[i]])
    for(j in 1:length(ind)){
      beta.init.list[[ind[j]]]<-beta.ols
    }
  }
  beta.init <- NULL
  for (i in 1:size) {
    bi <- beta.init.list[[i]]
    if(length(bi) > 0){
      beta.init <- c(beta.init,bi)
    } else {
      bi <- rep(0,q)
      beta.init <- c(beta.init,bi)
    }
  }
  beta.init <- as.matrix(beta.init)
  return(list(beta.init=beta.init))
}






