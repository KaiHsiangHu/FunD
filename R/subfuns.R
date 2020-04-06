#======FD Empirical======
datainf <- function(data, datatype){
  if(datatype == "abundance"){
    f1_10 = sapply(1:10, function(i) sum(data==i))
    a1 <- c(sum(data), sum(data != 0), f1_10)
    names(a1) <- c("n", "S.obs", paste0("f", 1:10))
    return(a1)
  }else if(datatype == 'incidence_freq'){
    f1_10 = sapply(1:10, function(i) sum(data[-1]==i))
    a1 <- matrix(c(data[1], sum(data[-1] != 0), f1_10), ncol = 1)
    names(a1) <- c("T", "S.obs", paste0("Q", 1:10))
    return(a1)
  }
}
FD_mle <- function(ai_vi, q){
  V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])
  n <- round(V_bar)
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    a <- ai_vi$ai[,i];v = ai_vi$vi[,i]
    sapply(q,function(qq){
      if(qq==1){
        exp(sum(-v*a/n*log(a/n)))
      }else{
        (sum(v*(a/n)^qq))^(1 / (1-qq))
      }
    })
  })
 
  matrix(out,nrow = length(q),ncol = ncol(ai_vi$ai))
}
EstiBootComm.Func = function(data, distance, datatype){
  #ifelse(datatype=="incidence_freq", data <- data[-1], data <- data)
  if (datatype=="incidence_freq") {
    n <- data[1]
    data <- data[-1]
    u=sum(data)
  } else if (datatype=="abundance") {
    n = sum(data)
    data=data
  }
  distance = as.matrix(distance)
  dij = distance[data!=0, data!=0]
  
  X = data[data>0]
  f1 <- sum(X == 1) ; f2 <- sum(X == 2)
  f0.hat <- ceiling(ifelse(f2>0, ((n-1)/n)*f1^2/2/f2, ((n-1)/n)*f1*(f1-1)/2))
  if (datatype=="abundance") {
    C1 = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
    W <- (1 - C1)/sum(X/n*(1-X/n)^n)
    Prob.hat.Unse <- rep((1-C1)/f0.hat, f0.hat)
  } else if (datatype=="incidence_freq") {
    C1 = ifelse(f2>0, 1-f1/u*(n-1)*f1/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/u/((n-1)*(f1-1)+2))
    W <- (1 - C1)/sum(X/u*(1-X/n)^n)
    Prob.hat.Unse <- rep(u/n*(1-C1)/f0.hat, f0.hat)
  }
  
  Prob.hat <- X/n*(1-W*(1-X/n)^n)
  Prob <- c(Prob.hat, Prob.hat.Unse)
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  #
  #F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(2*F.1)/2))
  #F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11)/(2 * n * (n-1))) )
  #
  if (datatype=="abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype=="incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }
  
  if (f0.hat==0) {
    d=dij
  } else if (f0.hat==1) {
    random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
    d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
    d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
    
    fo.num = (f0.hat * (f0.hat-1) )/2
    random_d00 = as.vector(rmultinom(1, 1000, rep(1/fo.num, fo.num) ) )/1000
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)*random_d00
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  
  return(list("pi" = Prob,"dij" = d))
}
FDtable_mle <- function(datalist, dij, tau, q, datatype, nboot = 30, conf = 0.95){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(datatype=='abundance'){
    out <- lapply(datalist, function(x){
      data_aivi <- data_transform(data = x,dij = dij,tau = tau)
      emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        n=sum(x)
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau)
          FD_mle(ai_vi = Boot_aivi,q = q) %>% as.numeric()
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype == 'incidence_freq'){
    output
    ### to be added
  }
  sites_tmp <- rep(sites,each = length(q)*length(tau))
  tau_tmp <- rep(rep(tau,each = length(q)),length(sites))
  Output <- tibble(order = rep(q,length(tau)*length(sites)), Empirical = out[,1],
                       LCL = out[,2], UCL = out[,3],
                       tau = tau_tmp,Community = sites_tmp)
  Output
}
AUCtable_mle <- function(datalist, dij, q = c(0,1,2), knots = 100, tau=NULL, datatype,
                         nboot=0, conf=0.95) {
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  dmin <- min(dij[dij>0])
  dmax <- max(dij)
  if(is.null(tau)){
    tau <- seq(dmin,dmax,length.out = knots)
  }
  #q_int <- c(0, 1, 2)
  AUC <- FDtable_mle(datalist,dij,tau,q,datatype,nboot = nboot,conf = conf)
  AUC <- AUC %>% group_by(Community,order) %>% 
    summarise(AUC_L = sum(Empirical[seq_along(Empirical[-1])]*diff(tau)),
              AUC_R = sum(Empirical[-1]*diff(tau)),
              AUC_LCL_L = sum(LCL[seq_along(LCL[-1])]*diff(tau)),
              AUC_LCL_R = sum(LCL[-1]*diff(tau)),
              AUC_UCL_L = sum(UCL[seq_along(UCL[-1])]*diff(tau)),
              AUC_UCL_R = sum(UCL[-1]*diff(tau))) %>% 
    mutate(AUC = (AUC_L+AUC_R)/2, AUC_LCL = (AUC_LCL_L+AUC_LCL_R)/2,
           AUC_UCL = (AUC_UCL_L+AUC_UCL_R)/2) %>% select(Community,order,AUC,AUC_LCL,AUC_UCL)
  AUC
}

#======FD Estimated======
data_transform <- function(data,dij,tau){
  dij <- dij[data>0,data>0]
  data <- data[data>0]
  out <- lapply(tau,function(tau_){
    dij_ <- dij
    dij_[which(dij_>tau_,arr.ind = T)] <- tau_
    a <- as.vector((1 - dij_/tau_) %*% data )
    data <- data[a!=0]
    a <- a[a!=0]
    v <- data/a
    cbind(a,v)
  }) 
  out_a <- matrix(sapply(out, function(x) x[,1]),ncol = length(tau))
  out_v <- matrix(sapply(out, function(x) x[,2]),ncol = length(tau))
  colnames(out_a) <- colnames(out_v) <- paste0('tau_',round(tau,3))
  # output <- array(data = 0,dim = c(nrow(out_a),ncol(out_a),2),dimnames = list(
  #   NULL,
  #   paste0('tau_',round(tau,3)),
  #   c('ai','vi')
  # ))
  # output[,,1] <- out_a
  # output[,,2] <- out_v
  output = list(ai = out_a, vi = out_v)
}
FD_est = function(ai_vi, q){ # ai_vi is array containing two elements: ai and vi
  V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])
  n <- round(V_bar)
  Sub <- function(q,FD_obs,n,f1,f2,h1,h2,A,av,avtab,deltas){
    if(q==0){
      ans <- FD_obs+FDq0(n,f1,f2,h1,h2,A)
    }else if(q==1){
      h_est_2 <- FDq1_1(n,h1,A)
      h_est_1 <- av %>% filter(ai<=(n-1)) %>% mutate(diga = digamma(n)-digamma(ai)) %>%
        apply(., 1, prod) %>% sum(.)/n
      ans <- exp(h_est_1+h_est_2)
    }else if(q==2){
      ans <- FDq2(as.matrix(avtab),n)
    }else{
      k <- 0:(n-1)
      a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
      b <- ifelse(h1==0|A==1,0,(h1*((1-A)^(1-n))/n)*(A^(q-1)-sum(choose(q-1,k)*(A-1)^k)))
      ans <- (a+b)^(1/(1-q))
    }
    return(ans)
  }
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    av = tibble(ai = ceiling(ai_vi$ai[,i]), vi = ai_vi$vi[,i])
    FD_obs <- sum(av[,2])
    f1 <- sum(av[,1]==1); h1 <- ifelse(f1>0,sum(av[av[,1]==1,2]),0)
    f2 <- sum(av[,1]==2); h2 <- ifelse(f2>0,sum(av[av[,1]==2,2]),0)
    if(f2 > 0){
      A = 2*f2/((n-1)*f1+2*f2)
    }else if(f2 == 0 & f1 > 0){
      A = 2/((n-1)*(f1-1)+2)
    }else{
      A = 1
    }
    if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
      avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
      deltas <- sapply(0:(n-1), function(k){
        del_tmp <- avtab[avtab[,1]<=(n-k),,drop=FALSE]
        delta(del_avtab = del_tmp,k,n)
      })
    }else{
     deltas <- 0
     avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
    }
    c(n,nrow(ai_vi$ai),FD_obs,f1,f2,h1,h2,
      sapply(q, function(qq) Sub(qq,FD_obs,n,f1,f2,h1,h2,A,av,avtab,deltas)))
  }) 
  out = matrix(out,ncol = ncol(ai_vi$ai))
  info = t(out[(1:7),])
  colnames(info) = c('n','S.obs','FD.obs','f1','f2','h1','h2')
  list(est = out[-(1:7),,drop=F], info = info)
}
FDtable_est <- function(datalist, dij, tau, q, datatype, nboot = 30, conf = 0.95){#change final list name
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(datatype=="abundance"){
    out <- lapply(datalist, function(x){
      data_aivi <- data_transform(data = x,dij = dij,tau = tau)
      est_info <- FD_est(ai_vi = data_aivi,q = q)
      est <- est_info$est %>% as.numeric()
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        n=sum(x)
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau)
          FD_est(ai_vi = Boot_aivi,q = q)$est %>% as.numeric()
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(est))
      }
      output <- cbind(est,est-qtile*ses,est+qtile*ses)
      output[output[,2]<0,2] <- 0
      list(estimates = output,info = est_info$info)
    }) 
  }else if(datatype=="incidence_freq"){
   # to be added
  }
  info <- lapply(out, function(x) x[[2]]) %>% 
    do.call(rbind,.) %>% as_tibble %>% 
    mutate(Community = rep(names(datalist),each = length(tau)),
           tau = rep(tau,length(datalist))) %>% 
    select(Community,tau,n,S.obs,f1,f2,h1,h2)
  
  
  Estoutput <- lapply(out, function(x) x[[1]]) %>%
    do.call(rbind,.) #%>% mutate(Community = rep(names(datalist),each = length(q)))
  sites_tmp <- rep(sites,each = length(q)*length(tau))
  tau_tmp <- rep(rep(tau,each = length(q)),length(sites))
  Estoutput <- tibble(order = rep(q,length(tau)*length(sites)), Estimated = Estoutput[,1],
                   LCL = Estoutput[,2], UCL = Estoutput[,3],
                   tau = tau_tmp,Community = sites_tmp)
  Estoutput$LCL[Estoutput$LCL<0] = 0
  return(list(Estoutput = Estoutput, info = info))
}
AUCtable_est <- function(datalist, dij, q = c(0,1,2), knots = 100, tau=NULL, datatype,
                         nboot=0, conf=0.95) {
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  dmin <- min(dij[dij>0])
  dmax <- max(dij)
  if(is.null(tau)){
    tau <- seq(dmin,dmax,length.out = knots)
  }
  #q_int <- c(0, 1, 2)
  AUC <- FDtable_est(datalist,dij,tau,q,datatype,nboot = nboot,conf = conf)$Estoutput
  AUC <- AUC %>% group_by(Community,order) %>% 
    summarise(AUC_L = sum(Estimated[seq_along(Estimated[-1])]*diff(tau)),
              AUC_R = sum(Estimated[-1]*diff(tau)),
              AUC_LCL_L = sum(LCL[seq_along(LCL[-1])]*diff(tau)),
              AUC_LCL_R = sum(LCL[-1]*diff(tau)),
              AUC_UCL_L = sum(UCL[seq_along(UCL[-1])]*diff(tau)),
              AUC_UCL_R = sum(UCL[-1]*diff(tau))) %>% 
    mutate(AUC = (AUC_L+AUC_R)/2, AUC_LCL = (AUC_LCL_L+AUC_LCL_R)/2,
           AUC_UCL = (AUC_UCL_L+AUC_UCL_R)/2) %>% select(Community,order,AUC,AUC_LCL,AUC_UCL)
  AUC
}

#===============iNextFD==================
FD.m.est = function(ai_vi, m, q){
  EFD = function(m,qs,obs,asy,beta,av){
    m = m-n
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2 ) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 ){
        1/sum( (av[,2])*((1/(n+m))*(av[,1]/n)+((n+m-1)/(n+m))*(av[,1]*(av[,1]-1)/(n*(n-1)))) )
      } 
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])
  n <- round(V_bar)
  asy <- FD_est(ai_vi,q)$est
  obs <- FD_mle(ai_vi,q)
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    av = cbind(ai = ceiling(ai_vi$ai[,i]), vi =  ai_vi$vi[,i])
    RFD_m = RFD(av, n, n-1, q)
    beta <- rep(0,length(q))
    #asymptotic value; observed value
    asy_i <- asy[,i];obs_i <- obs[,i]
    beta0plus <- which( asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus]-RFD_m[beta0plus])/(asy_i[beta0plus]-RFD_m[beta0plus])
    sapply(m, function(mm){
      if(mm<n){
        RFD(av,n,mm,q) 
      }else if(mm==n){
        obs_i
      }else{
        EFD(m = mm,qs = q,obs = obs_i,asy = asy_i,beta = beta,av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out,ncol = ncol(ai_vi$ai))
}
Coverage = function(data, datatype, m){
  n <- ifelse(datatype=='incidence_freq', data[1], sum(data) )
  if(datatype == "incidence_freq"){
    x <- data[-1]
    u<-sum(x)
  }else if(datatype == "abundance"){
    x <- data
  }
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  Sub2 <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/u*A
    if(m > n) out <- 1-f1/u*A^(m-n+1)
    out
  }
  sapply(m, FUN = function(i){
    ifelse(datatype=='incidence', Sub2(i), Sub(i) )
  })
}

iNextFD = function(datalist, dij, q = c(0,1,2), datatype, tau, nboot, conf, m){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  #if(datatype=="abundance") ns <- sapply(datalist, sum)
  #if(datatype=='incidence_freq') ns <- sapply(datalist, function(y) y[1])
  if(datatype=="abundance"){
    length_tau <- length(tau)
    out <- lapply(1:length(datalist), function(i){
      x <- datalist[[i]]
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau)
      qFDm <- FD.m.est(ai_vi = data_aivi,m = m[[i]],q = q) %>%
         as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau)
          qFDm_b <- FD.m.est(ai_vi = Boot_aivi,m = m[[i]],q = q) %>%
            t() %>% as.numeric()
          covm_b = Coverage(Boot.X[,B], datatype, m[[i]])
          return(c(qFDm_b,covm_b))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(c(qFDm,covm)))
      }
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      method <- rep(method,length(q)*length(tau))
      orderq <- rep(rep(q,each = length(m[[i]])),length(tau))
      threshold <- rep(tau,each = length(q)*length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_cov <- rep(ses_cov,each = length(q)*length(tau))
      ses_fd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      covm <- rep(covm,length(q)*length(tau))
      tibble(m=rep(m[[i]],length(q)*length(tau)),method=method,order=orderq,
             qFD=qFDm,qFD.LCL=qFDm-qtile*ses_fd,qFD.UCL=qFDm+qtile*ses_fd,
             SC=covm,SC.LCL=covm-qtile*ses_cov,SC.UCL=covm+qtile*ses_cov,
             site = sites[i],threshold = threshold) %>% 
      arrange(threshold,order)
    }) %>% do.call(rbind, .)
   
  }else if(datatype=="incidence_freq"){
    # to be added
  }
  return(out)
}

