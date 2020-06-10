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
#===== old version ======
# FD_mle <- function(ai_vi, q){
#   V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])
#   n <- round(V_bar)
#   out <- sapply(1:ncol(ai_vi$ai), function(i){
#     a <- ceiling(ai_vi$ai[,i]);v = ai_vi$vi[,i]
#     sapply(q,function(qq){
#       if(qq==1){
#         exp(sum(-v*a/n*log(a/n)))
#       }else{
#         (sum(v*(a/n)^qq))^(1 / (1-qq))
#       }
#     })
#   })
#  
#   matrix(out,nrow = length(q),ncol = ncol(ai_vi$ai))
# }
#===== new version ======
FD_mle <- function(ai_vi, q){
  v_bar <- round(sum(ai_vi$ai[,1]*ai_vi$vi[,1]))
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    a <- ai_vi$ai[,i]
    a[a<1] <- 1
    a <- round(a)
    v = ai_vi$vi[,i]
    sapply(q,function(qq){
      if(qq==1){
        exp(sum(-v*a/v_bar*log(a/v_bar)))
      }else{
        (sum(v*(a/v_bar)^qq))^(1 / (1-qq))
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
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
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
    out <- lapply(datalist, function(x){
      nT = x[1]
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          FD_mle(ai_vi = Boot_aivi,q = q) %>% as.numeric()
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
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
  # dmin <- min(dij[dij>0])
  # dmax <- max(dij)
  # if(is.null(tau)){
  #   tau <- seq(dmin,dmax,length.out = knots)
  # }
  if(is.null(tau)){
    tau <- seq(0,1,length.out = knots)
  }
  #q_int <- c(0, 1, 2)
  
  AUC <- FDtable_mle(datalist,dij,tau,q,datatype,nboot = 0) %>%
    group_by(Community,order) %>% 
    summarise(AUC_L = sum(Empirical[seq_along(Empirical[-1])]*diff(tau)),
              AUC_R = sum(Empirical[-1]*diff(tau))) %>% ungroup %>% 
    mutate(AUC = (AUC_L+AUC_R)/2) %>% select(Community,order,AUC)
  if(datatype=='abundance'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        ses <- FDtable_mle(Boot.X,dij_boot,tau,q,datatype,nboot = 0) %>% 
          group_by(Community,order) %>% 
          summarise(AUC_L = sum(Empirical[seq_along(Empirical[-1])]*diff(tau)),
                    AUC_R = sum(Empirical[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Community = rep(sites,each = length(q)))
    }else{
      ses <- tibble(order = rep(q,length(datalist)), se = NA, Community = rep(sites,each = length(q)))
    }
  }else if(datatype=='incidence_freq'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>%
          split(., col(.))
        ses <- FDtable_mle(Boot.X,dij_boot,tau,q,datatype,nboot = 0) %>% 
          group_by(Community,order) %>% 
          summarise(AUC_L = sum(Empirical[seq_along(Empirical[-1])]*diff(tau)),
                    AUC_R = sum(Empirical[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Community = rep(sites,each = length(q)))
    }else{
      ses <- tibble(order = rep(q,length(datalist)), se = NA, Community = rep(sites,each = length(q)))
    }
  }
 
  AUC <- left_join(x = AUC, y = ses, by = c('Community','order')) %>% mutate(LCL = AUC - se * qtile,
                                                                    UCL = AUC + se * qtile) %>% 
    select(-se)
  AUC$LCL[AUC$LCL<0] <- 0
  AUC
}

#======FD Estimated======
#===== old version ======
# data_transform <- function(data,dij,tau){
#   dij <- dij[data>0,data>0]
#   data <- data[data>0]
#   out <- lapply(tau,function(tau_){
#     dij_ <- dij
#     if(tau_==0){
#       dij_[dij_>0] <- 1
#       a <- as.vector((1 - dij_/1) %*% data )
#     }else{
#       dij_[which(dij_>tau_,arr.ind = T)] <- tau_
#       a <- as.vector((1 - dij_/tau_) %*% data )
#     }
#     data <- data[a!=0]
#     a <- a[a!=0]
#     v <- data/a
#     cbind(a,v)
#   })
#   out_a <- matrix(sapply(out, function(x) x[,1]),ncol = length(tau))
#   out_v <- matrix(sapply(out, function(x) x[,2]),ncol = length(tau))
#   colnames(out_a) <- colnames(out_v) <- paste0('tau_',round(tau,3))
# 
#   
#   
#   # output <- array(data = 0,dim = c(nrow(out_a),ncol(out_a),2),dimnames = list(
#   #   NULL,
#   #   paste0('tau_',round(tau,3)),
#   #   c('ai','vi')
#   # ))
#   # output[,,1] <- out_a
#   # output[,,2] <- out_v
#   output = list(ai = out_a, vi = out_v)
#   output
# }
#===== new version ======
data_transform <- function(data,dij,tau,datatype,truncate = TRUE){
  if(datatype == 'abundance'){
    dij <- dij[data>0,data>0]
    data <- data[data>0]
    out <- lapply(tau,function(tau_){
      dij_ <- dij
      if(tau_==0){
        dij_[dij_>0] <- 1
        a <- as.vector((1 - dij_/1) %*% data )
      }else{
        dij_[which(dij_>tau_,arr.ind = T)] <- tau_
        a <- as.vector((1 - dij_/tau_) %*% data )
      }
      data <- data[a!=0]
      a <- a[a!=0]
      v <- data/a
      cbind(a,v)
    })
  }else if (datatype == 'incidence_freq'){
    nT = data[1]
    data <- data[-1]
    dij <- dij[data>0,data>0]
    data <- data[data>0]
    if(truncate == TRUE){
      out <- lapply(tau,function(tau_){
        dij_ <- dij
        if(tau_==0){
          dij_[dij_>0] <- 1
          a <- as.vector((1 - dij_/1) %*% data )
        }else{
          dij_[which(dij_>tau_,arr.ind = T)] <- tau_
          a <- as.vector((1 - dij_/tau_) %*% data )
        }
        data <- data[a!=0]
        a <- a[a!=0]
        a[a>nT] <- nT
        v <- data/a
        cbind(a,v)
      })
    }else{
      out <- lapply(tau,function(tau_){
        dij_ <- dij
        if(tau_==0){
          dij_[dij_>0] <- 1
          a <- as.vector((1 - dij_/1) %*% data )
        }else{
          dij_[which(dij_>tau_,arr.ind = T)] <- tau_
          a <- as.vector((1 - dij_/tau_) %*% data )
        }
        data <- data[a!=0]
        a <- a[a!=0]
        v <- data/a
        cbind(a,v)
      })
    }
  }
  out_a <- matrix(sapply(out, function(x) x[,1]),ncol = length(tau))
  out_v <- matrix(sapply(out, function(x) x[,2]),ncol = length(tau))
  colnames(out_a) <- colnames(out_v) <- paste0('tau_',round(tau,3))
  
  output = list(ai = out_a, vi = out_v)
  output
}

#=====old version=====
# FD_est = function(ai_vi, q){ # ai_vi is array containing two elements: ai and vi
# 
#   V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])
#   n <- round(V_bar)
#   Sub <- function(q,FD_obs,n,f1,f2,h1,h2,A,av,avtab,deltas){
#     if(q==0){
#       ans <- FD_obs+FDq0(n,f1,f2,h1,h2,A)
#     }else if(q==1){
#       h_est_2 <- FDq1_1(n,h1,A)
#       h_est_1 <- av %>% filter(ai<=(n-1)) %>% mutate(diga = digamma(n)-digamma(ai)) %>%
#         apply(., 1, prod) %>% sum(.)/n
#       ans <- exp(h_est_1+h_est_2)
#     }else if(q==2){
#       ans <- FDq2(as.matrix(avtab),n)
#     }else{
#       k <- 0:(n-1)
#       a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
#       b <- ifelse(h1==0|A==1,0,(h1*((1-A)^(1-n))/n)*(A^(q-1)-sum(choose(q-1,k)*(A-1)^k)))
#       ans <- (a+b)^(1/(1-q))
#     }
#     return(ans)
#   }
#   out <- sapply(1:ncol(ai_vi$ai), function(i){
#     av = tibble(ai = ceiling(ai_vi$ai[,i]), vi = ai_vi$vi[,i])
#     FD_obs <- sum(av[,2])
#     f1 <- sum(av[,1]==1); h1 <- ifelse(f1>0,sum(av[av[,1]==1,2]),0)
#     f2 <- sum(av[,1]==2); h2 <- ifelse(f2>0,sum(av[av[,1]==2,2]),0)
#     if(f2 > 0){
#       A = 2*f2/((n-1)*f1+2*f2)
#     }else if(f2 == 0 & f1 > 0){
#       A = 2/((n-1)*(f1-1)+2)
#     }else{
#       A = 1
#     }
#     if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
#       avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
#       deltas <- sapply(0:(n-1), function(k){
#         del_tmp <- avtab[avtab[,1]<=(n-k),,drop=FALSE]
#         delta(del_avtab = del_tmp,k,n)
#       })
#     }else{
#      deltas <- 0
#      avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
#     }
#     c(n,nrow(ai_vi$ai),FD_obs,f1,f2,h1,h2,
#       sapply(q, function(qq) Sub(qq,FD_obs,n,f1,f2,h1,h2,A,av,avtab,deltas)))
#   }) 
#   out = matrix(out,ncol = ncol(ai_vi$ai))
#   info = t(out[(1:7),])
#   colnames(info) = c('n','S.obs','FD.obs','f1','f2','h1','h2')
#   list(est = out[-(1:7),,drop=F], info = info)
# }
#=====new version=====
FD_est = function(ai_vi, q, nT){ # ai_vi is array containing two elements: ai and vi
  V_bar <- round(sum(ai_vi$ai[,1]*ai_vi$vi[,1]))/nT
  Sub <- function(q,FD_obs,nT,f1,f2,h1,h2,A,av,avtab,deltas){
    if(q==0){
      ans <- FD_obs+FDq0(nT,f1,f2,h1,h2,A)
    }else if(q==1){
      h_est_2 <- FDq1_1(nT,h1,A)
      h_est_1 <- av %>% filter(ai<=(nT-1)) %>% mutate(diga = digamma(nT)-digamma(ai)) %>%
        apply(., 1, prod) %>% sum(.)/nT
      ans <- V_bar*exp((h_est_1+h_est_2)/V_bar)
    }else if(q==2){
      ans <- FDq2(as.matrix(avtab),nT)*V_bar^2
    }else{
      k <- 0:(nT-1)
      a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
      b <- ifelse(h1==0|A==1,0,(h1*((1-A)^(1-nT))/nT)*(A^(q-1)-sum(choose(q-1,k)*(A-1)^k)))
      ans <- ((a+b)/(V_bar^q))^(1/(1-q))
    }
    return(ans)
  }
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    ai <- ai_vi$ai[,i]
    ai[ai<1] <- 1
    av = tibble(ai = round(ai), vi = ai_vi$vi[,i])
    FD_obs <- sum(av[,2])
    f1 <- sum(av[,1]==1); h1 <- ifelse(f1>0,sum(av[av[,1]==1,2]),0)
    f2 <- sum(av[,1]==2); h2 <- ifelse(f2>0,sum(av[av[,1]==2,2]),0)
    if(f2 > 0){
      A = 2*f2/((nT-1)*f1+2*f2)
    }else if(f2 == 0 & f1 > 0){
      A = 2/((nT-1)*(f1-1)+2)
    }else{
      A = 1
    }
    if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
      avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
      deltas <- sapply(0:(nT-1), function(k){
        del_tmp <- avtab[avtab[,1]<=(nT-k),,drop=FALSE]
        delta(del_avtab = del_tmp,k,nT)
      })
    }else{
      deltas <- 0
      avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
    }
    c(nT,nrow(ai_vi$ai),FD_obs,f1,f2,h1,h2,
      sapply(q, function(qq) Sub(qq,FD_obs,nT,f1,f2,h1,h2,A,av,avtab,deltas)))
  }) 
  out = matrix(out,ncol = ncol(ai_vi$ai))
  info = t(out[(1:7),])
  colnames(info) = c('nT','S.obs','FD.obs','f1','f2','h1','h2')
  list(est = out[-(1:7),,drop=F], info = info)
}
FDtable_est <- function(datalist, dij, tau, q, datatype, nboot = 30, conf = 0.95){#change final list name
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(datatype=="abundance"){
    out <- lapply(datalist, function(x){
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      est_info <- FD_est(ai_vi = data_aivi,q = q,nT = n)
      est <- est_info$est %>% as.numeric()
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
          FD_est(ai_vi = Boot_aivi,q = q,nT = n)$est %>% as.numeric()
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(est))
      }
      output <- cbind(est,est-qtile*ses,est+qtile*ses)
      output[output[,2]<0,2] <- 0
      list(estimates = output,info = est_info$info)
    }) 
  }else if(datatype=="incidence_freq"){
    out <- lapply(datalist, function(x){
      nT=x[1]
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      est_info <- FD_est(ai_vi = data_aivi,q = q,nT = nT)
      est <- est_info$est %>% as.numeric()
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          FD_est(ai_vi = Boot_aivi,q = q,nT = nT)$est %>% as.numeric()
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(est))
      }
      output <- cbind(est,est-qtile*ses,est+qtile*ses)
      output[output[,2]<0,2] <- 0
      list(estimates = output,info = est_info$info)
    })
  }
  info <- lapply(out, function(x) x[[2]]) %>% 
    do.call(rbind,.) %>% as_tibble %>% 
    mutate(Community = rep(names(datalist),each = length(tau)),
           tau = rep(tau,length(datalist))) %>% 
    select(Community,tau,nT,S.obs,f1,f2,h1,h2)
  
  
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
  # dmin <- min(dij[dij>0])
  # dmax <- max(dij)
  # if(is.null(tau)){
  #   tau <- seq(dmin,dmax,length.out = knots)
  # }
  if(is.null(tau)){
    tau <- seq(0,1,length.out = knots)
  }
  #q_int <- c(0, 1, 2)
  AUC <- FDtable_est(datalist,dij,tau,q,datatype,nboot = 0)$Estoutput %>%
    group_by(Community,order) %>% 
    summarise(AUC_L = sum(Estimated[seq_along(Estimated[-1])]*diff(tau)),
              AUC_R = sum(Estimated[-1]*diff(tau))) %>% ungroup %>% 
    mutate(AUC = (AUC_L+AUC_R)/2) %>% select(Community,order,AUC)
  if(datatype=='abundance'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        ses <- FDtable_est(Boot.X,dij_boot,tau,q,datatype,nboot = 0)$Estoutput %>% 
          group_by(Community,order) %>% 
          summarise(AUC_L = sum(Estimated[seq_along(Estimated[-1])]*diff(tau)),
                    AUC_R = sum(Estimated[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Community = rep(sites,each = length(q)))
    }else{
      ses <- tibble(order = rep(q,length(datalist)), se = NA, Community = rep(sites,each = length(q)))
    }
  }else if(datatype=='incidence_freq'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>%
          split(., col(.))
        ses <- FDtable_est(Boot.X,dij_boot,tau,q,datatype,nboot = 0)$Estoutput %>% 
          group_by(Community,order) %>% 
          summarise(AUC_L = sum(Estimated[seq_along(Estimated[-1])]*diff(tau)),
                    AUC_R = sum(Estimated[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Community = rep(sites,each = length(q)))
    }else{
      ses <- tibble(order = rep(q,length(datalist)), se = NA, Community = rep(sites,each = length(q)))
    }
  }
  AUC <- left_join(x = AUC, y = ses, by = c('Community','order')) %>% mutate(LCL = AUC - se * qtile,
                                                                             UCL = AUC + se * qtile) %>% 
    select(-se)
  AUC$LCL[AUC$LCL<0] <- 0
  AUC
}

#===============iNextFD==================
FD.m.est = function(ai_vi, m, q, nT){
  EFD = function(m,qs,obs,asy,beta,av){
    m = m-nT
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2 ) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 ){
        V_bar^2/sum( (av[,2])*((1/(nT+m))*(av[,1]/nT)+((nT+m-1)/(nT+m))*(av[,1]*(av[,1]-1)/(nT*(nT-1)))) )
      } 
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])/nT
  asy <- FD_est(ai_vi,q,nT)$est
  obs <- FD_mle(ai_vi,q)
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    ai <- ai_vi$ai[,i]
    ai[ai<1] <- 1
    av = cbind(ai = round(ai), vi = ai_vi$vi[,i])
    # av = cbind(ai = ceiling(ai_vi$ai[,i]), vi =  ai_vi$vi[,i])
    RFD_m = RFD(av, nT, nT-1, q,V_bar)
    beta <- rep(0,length(q))
    #asymptotic value; observed value
    asy_i <- asy[,i];obs_i <- obs[,i]
    asy_i <- sapply(1:length(q), function(j){
      max(asy_i[j],obs_i[j])
    })
    beta0plus <- which( asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus]-RFD_m[beta0plus])/(asy_i[beta0plus]-RFD_m[beta0plus])
    sapply(m, function(mm){
      if(mm<nT){
        RFD(av,nT,mm,q,V_bar) 
      }else if(mm==nT){
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
  x <- x[x>0]
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
    ifelse(datatype!='abundance', Sub2(i), Sub(i) )
  })
}

iNextFD = function(datalist, dij, q = c(0,1,2), datatype, tau, nboot, conf = 0.95, m){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  length_tau <- length(tau)
  #if(datatype=="abundance") ns <- sapply(datalist, sum)
  #if(datatype=='incidence_freq') ns <- sapply(datalist, function(y) y[1])
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      x <- datalist[[i]]
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      qFDm <- FD.m.est(ai_vi = data_aivi,m = m[[i]],q = q,nT = n) %>%
         as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
          qFDm_b <- FD.m.est(ai_vi = Boot_aivi,m = m[[i]],q = q,nT = n) %>%
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
    out <- lapply(1:length(datalist), function(i){
      x <- datalist[[i]]
      nT=x[1]
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      qFDm <- FD.m.est(ai_vi = data_aivi,m = m[[i]],q = q,nT = nT) %>%
        as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          qFDm_b <- FD.m.est(ai_vi = Boot_aivi,m = m[[i]],q = q,nT = nT) %>%
            t() %>% as.numeric()
          covm_b = Coverage(Boot.X, datatype, m[[i]])
          return(c(qFDm_b,covm_b))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(c(qFDm,covm)))
      }
      method <- ifelse(m[[i]]>nT,'Extrapolation',ifelse(m[[i]]<nT,'Rarefaction','Observed'))
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
  }
  return(out)
}
AUCtable_iNextFD <- function(datalist, dij, q = c(0,1,2), knots = 100, datatype, tau=NULL,
                         nboot=0, conf=0.95, m) {
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  # dmin <- min(dij[dij>0])
  # dmax <- max(dij)
  # if(is.null(tau)){
  #   tau <- seq(dmin,dmax,length.out = knots)
  # }
  if(is.null(tau)){
    tau <- seq(0,1,length.out = knots)
  }
  AUC <- iNextFD(datalist,dij,q,datatype,tau,nboot = 0,m = m) %>%
    group_by(site,order,m) %>% 
    summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
              AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC),method = unique(method)) %>% ungroup %>% 
    mutate(AUC = (AUC_L+AUC_R)/2) %>% select(site,order,m,AUC,SC,method)
  if(datatype == 'abundance'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Community_ <- rep(sites[[i]],length(q)*length(m[[i]]))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        m_boot <- lapply(1:nboot, function(b) m[[i]])
        ses <- iNextFD(Boot.X,dij_boot,q,datatype,tau,nboot = 0,m = m_boot) %>% 
          group_by(site,order,m) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order,m) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(site = Community_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(site,order,m) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }else if (datatype == 'incidence_freq'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Community_ <- rep(sites[[i]],length(q)*length(m[[i]]))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>% 
          split(., col(.))
        m_boot <- lapply(1:nboot, function(b) m[[i]])
        ses <- iNextFD(Boot.X,dij_boot,q,datatype,tau,nboot = 0,m = m_boot) %>% 
          group_by(site,order,m) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order,m) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(site = Community_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(site,order,m) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }
  
  AUC <- left_join(x = AUC, y = ses, by = c('site','order','m')) %>% mutate(
    AUC.LCL = AUC - AUC_se * qtile, AUC.UCL = AUC + AUC_se * qtile,
    SC.LCL = SC - SC_se * qtile, SC.UCL = SC + SC_se * qtile) %>% 
    select(site,order,m,method,AUC,AUC.LCL,AUC.UCL,SC.LCL,SC.UCL)
  AUC$AUC.LCL[AUC$AUC.LCL<0] <- 0
  AUC$SC.LCL[AUC$SC.LCL<0] <- 0
  AUC
}
#======EstimateFD========
invChatFD_abu <- function(ai_vi, data_, q, Cs, tau){
  n <- sum(data_)
  refC = Coverage(data_, 'abundance', n)
  f <- function(m, cvrg) abs(Coverage(data_, 'abundance', m) - cvrg)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, cvrg = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
      mm <- round(mm)
    }else if (refC <= cvrg) {
      f1 <- sum(data_ == 1)
      f2 <- sum(data_ == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - cvrg))/log(A) - 1
      mm <- n + mm
      mm <- round(mm)
    }
    mm
  })
  mm[mm==0] <- 1
  SC <- Coverage(data_, 'abundance', mm)
  out <- FD.m.est(ai_vi = ai_vi,m = mm,q = q,nT = n)
  out <- as.vector(out)
  method <- ifelse(mm>n,'Extrapolated',ifelse(mm<n,'Interpolated','Observed'))
  method <- rep(method,length(q)*length(tau))
  m <- rep(mm,length(q)*length(tau))
  order <- rep(rep(q,each = length(mm)),length(tau))
  SC <- rep(SC,length(q)*length(tau))
  goalSC <- rep(Cs,length(q)*length(tau))
  threshold <- rep(tau,each = length(q)*length(mm))
  tibble(m = m,method = method,order = order,
         qFD = out,SC=SC,goalSC = goalSC, threshold = threshold)
}

invChatFD_inc <- function(ai_vi, data_, q, Cs, tau){
  n <- data_[1]
  refC = Coverage(data_, 'incidence_freq', n)
  f <- function(m, cvrg) abs(Coverage(data_, 'incidence_freq', m) - cvrg)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, cvrg = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
      mm <- round(mm)
    }else if (refC <= cvrg) {
      f1 <- sum(data_ == 1)
      f2 <- sum(data_ == 2)
      U <- sum(data_)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - cvrg))/log(A) - 1
      mm <- n + mm
      mm <- round(mm)
    }
    mm
  })
  mm[mm==0] <- 1
  SC <- Coverage(data_, 'incidence_freq', mm)
  out <- FD.m.est(ai_vi = ai_vi,m = mm,q = q,nT = n)
  out <- as.vector(out)
  method <- ifelse(mm>n,'Extrapolated',ifelse(mm<n,'Interpolated','Observed'))
  method <- rep(method,length(q)*length(tau))
  m <- rep(mm,length(q)*length(tau))
  order <- rep(rep(q,each = length(mm)),length(tau))
  SC <- rep(SC,length(q)*length(tau))
  goalSC <- rep(Cs,length(q)*length(tau))
  threshold <- rep(tau,each = length(q)*length(mm))
  tibble(m = m,method = method,order = order,
         qFD = out,SC=SC,goalSC = goalSC, threshold = threshold)
}
invChatFD <- function(datalist, dij, q, datatype, level, nboot, conf = 0.95, tau){
  qtile <- qnorm(1-(1-conf)/2)

  if(datatype=='abundance'){
    out <- lapply(datalist,function(x_){
      data_aivi <- data_transform(data = x_,dij = dij,tau = tau,datatype = datatype)
      #n_sp_samp <- sum(aL_table$tgroup=='Tip')
      est <- invChatFD_abu(ai_vi = data_aivi,data_ = x_,q = q,Cs = level,tau = tau)
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x_,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        n=sum(x_)
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
          invChatFD_abu(ai_vi = Boot_aivi,data_ = Boot.X[,B],q = q,Cs = level,tau = tau)$qFD
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,nrow(est))
      }
      est <- est %>% mutate(qFD.LCL=qFD-qtile*ses,qFD.UCL=qFD+qtile*ses) 
    }) %>% do.call(rbind,.)
  }else if(datatype=='incidence_freq'){
    out <- lapply(datalist,function(x_){
      nT=x_[1]
      data_aivi <- data_transform(data = x_,dij = dij,tau = tau,datatype = datatype)
      est <- invChatFD_inc(ai_vi = data_aivi,data_ = x_,q = q,Cs = level,tau = tau)
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x_,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          invChatFD_inc(ai_vi = Boot_aivi,data_ = Boot.X,q = q,Cs = level,tau = tau)$qFD
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,nrow(est))
      }
      est <- est %>% mutate(qFD.LCL=qFD-qtile*ses,qFD.UCL=qFD+qtile*ses) 
    }) %>% do.call(rbind,.)
  }
  Community = rep(names(datalist), each = length(q)*length(level)*length(tau))
  out <- out %>% mutate(site = Community) %>% select(
    site,m,method,order,SC,qFD,qFD.LCL,qFD.UCL,goalSC,threshold
  )
  rownames(out) <- NULL
  out
}
AUCtable_invFD <- function(datalist, dij, q = c(0,1,2), knots = 100, datatype, level, nboot = 0, conf = 0.95, tau=NULL){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(is.null(tau)){
    tau <- seq(0,1,length.out = knots)
  }
  AUC <- invChatFD(datalist,dij,q,datatype,level,nboot = 0,tau = tau) %>%
    group_by(site,order,goalSC) %>% 
    summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
              AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC),m = mean(m),method = unique(method)) %>% 
    ungroup %>% mutate(AUC = (AUC_L+AUC_R)/2) %>% select(site,order,goalSC,m,method,AUC,SC)
  if(datatype == 'abundance'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Community_ <- rep(sites[[i]],length(q)*length(level))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        ses <- invChatFD(Boot.X,dij_boot,q,datatype,level,nboot = 0,tau = tau) %>% 
          group_by(site,order,goalSC) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order,goalSC) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(site = Community_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(site,order,goalSC) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }else if(datatype == 'incidence_freq'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Community_ <- rep(sites[[i]],length(q)*length(level))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>% 
          split(., col(.))
        ses <- invChatFD(Boot.X,dij_boot,q,datatype,level,nboot = 0,tau = tau) %>% 
          group_by(site,order,goalSC) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(order,goalSC) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(site = Community_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(site,order,goalSC) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }
  
  AUC <- left_join(x = AUC, y = ses, by = c('site','order','goalSC')) %>% mutate(
    AUC.LCL = AUC - AUC_se * qtile, AUC.UCL = AUC + AUC_se * qtile,
    SC.LCL = SC - SC_se * qtile, SC.UCL = SC + SC_se * qtile) %>% 
    select(site,order,goalSC,m,method,AUC,AUC.LCL,AUC.UCL,SC.LCL,SC.UCL)
  AUC$AUC.LCL[AUC$AUC.LCL<0] <- 0
  AUC$SC.LCL[AUC$SC.LCL<0] <- 0
  AUC
}
