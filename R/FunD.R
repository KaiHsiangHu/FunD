#' Interpolation and extrapolation of functional diversity
#'
#' \code{iNEXTFD}: the seamless rarefaction and extrapolation sampling curves of functional diversity(FD) for q = 0, 1 and 2.
#' See Chao et al. (2019) for pertinent background and methods.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr Type (1) abundance data: a S by N matrix/data.frame
#' where N is the number of assemblages. The element in i-th row and k-th is the abundance of species i in assemblage k. Please note
#' that the rownames of data must be the species names matching the species names in distance matrix and thus can't be empty.\cr
#' Type (2) incidence frequency data: the sampling unit is quadrat or transect, the observed species was only recorded as presence(detection)/absence(non-detection)
#' data in each sampling unit. Likewise, the rownames of data must be the species names matching the species names in phylogeny tree and thus can't be empty. \cr
#' @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
#' @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
#' @param endpoint an positive interger specifying the endpoint for rarefaction and
#' extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum reference sample size. It will be ignored if \code{size} is given. \cr
#' @param knots a positive integer specifying the number of knots between 1 and the \code{endpoint}. Default is 40.\cr
#' @param size a sequence of positive integers specifying the sample sizes for which PD estimates will be calculated. If \code{NULL}, then PD estimates will be
#' calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}. \cr
#' @param plot.type a positive integer vector specifying types of curves. Three types of plots: sample-size-based rarefaction and extrapolation curve (\code{plot.type = 1});
#' coverage-based rarefaction and extrapolation curve (\code{plot.type = 2}); sample completeness curve (\code{plot.type = 3}). Default is \code{c(1,2,3)}. \cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param threshold a positive value between 0 and 1 specifying tau. If \code{NULL}, \code{threshold} = (dmean+dmin)/2. Default is \code{NULL}.
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rmultinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @return a list of one objects: \cr\cr
#' \code{$inext} a table of FD estimates and sample completeness for interpolated or extrapolated sample sizes along with their confidence intervals (if \code{nboot > 0}). \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data[-1,]
#' dij <-  FunDdata.inc$dij
#' out <- iNEXTFD(data = data, distM = dij,datatype = "abundance")
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). functional diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of functional diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive functional diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
iNEXTFD <- function(data, distM, datatype = "abundance", q = c(0,1,2), endpoint = NULL, 
                    knots = 40, size = NULL, plot.type = 1:3, conf = 0.95, nboot = 50, threshold = NULL) {

  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,]
  }
  data <- data[rowSums(data)>0,,drop=FALSE]
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]

  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }

  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp <- rowMeans(matrix(sapply(dat, function(x) x/sum(x)),ncol = length(dat)))  
    }else if(datatype=='incidence_freq'){
      tmp <- rowMeans(matrix(sapply(dat, function(x) x[-1]/sum(x[-1])), ncol = length(dat)))
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  }else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to (dmean+dmin)/2.",call. = FALSE)
  }
  

  # if(class(mydata) == "list"){
  #   infos <- sapply(mydata, function(x){
  #     datainf(data = x, datatype, phylotr = mytree,reft = reft)})
  # }else{
  #   return(NULL)
  # }

  ############output2
  if(length(knots)!=length(dat)) knots <- rep(knots,length(dat))
  if(is.null(size)){
    if(is.null(endpoint)){
      if(datatype == "abundance") {
        endpoint <- sapply(dat, function(x) 2*sum(x))
      }else if(datatype == "incidence_freq"){
        endpoint <- sapply(dat, function(x) 2*x[1])
      }
    }else{
      if(length(endpoint)!=length(dat)){
        endpoint <- rep(endpoint,length(dat))
      }
    }
    size <- lapply(1:length(dat),function(i){
      ni <- sum(dat[[i]])
      if(endpoint[i] <= ni){
        mi <- floor(seq(1,endpoint[i],length.out = knots[i]))
      }else{
        mi <- floor(c(seq(1,ni,length.out = floor(knots[i]/2)),
                      seq(ni+1,endpoint[i],length.out = knots[i]-floor(knots[i]/2))))
      }
      unique(mi)
    })
  }else{
    if(class(size)=="numeric"|class(size)=="integer"|class(size)=="double"){
      size <- list(size = size)
    } 
    if(length(size)!=length(dat)) lapply(1:length(dat), function(x) size[[1]])
    size <- lapply(1:length(dat),function(i){
      ni <- sum(dat[[i]])
      if( sum(size[[i]] == ni) == 0 ) mi <- sort(c(ni,size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }

  FUN <- function(e){
    if(class(dat)=="list"){
      temp = iNextFD(datalist = dat,dij = distM,q = q,datatype = datatype,tau = threshold,
                     nboot = nboot,conf = conf,m = size)
      temp$qFD.LCL[temp$qPD.LCL<0] <- 0;temp$SC.LCL[temp$SC.LCL<0] <- 0
      temp$SC.UCL[temp$SC.UCL>1] <- 1
      return(temp)
    }else{
      return(NULL)
    }
  }
  RE.table <- tryCatch(FUN(e), error = function(e){return()})
  ###############output3
  # TYPE <- 1:3
  # FUN2 <- function(e){
  #   if(is.na(sum(pmatch(plot.type, TYPE))) == F){
  #     temp2 <- lapply(plot.type, function(j) RE_plot(RE.table, datatype, j))
  #     allname <- c("RE.plot.size", "RE.plot.C", "RE.plot.sizeC")
  #     names(temp2) <- allname[plot.type]
  #     temp2
  #   }else{
  #     return("invalid plot type", call. = FALSE)
  #   }
  # }
  # RE.plot <- tryCatch(FUN2(e), error = function(e){return()})

  # ans <- list(summary = infos,reference_time = reft, inext = RE.table, figure = RE.plot)
  ans <- list(inext = RE.table)
  ans
  
}

#' Compute functional diversity with particular sample coverages
#'
#'\code{EstimateFD}: computes functional diversity(FD) with particular user-specified levels of sample coverages.
#' See Chao et al. (2019) for pertinent background and methods.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTFD}} for data details.
#' @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
#' @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
#' @param level a positive sequence < 1 specifying a particular values of sample coverages.
#' If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each site to its double sample size. Default is \code{NULL}.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param threshold a positive sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}.
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rmultinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats optimize
#' @return a table including the sample size, sample coverage,
#' method (Interpolated or Extrapolated), and diversity estimates with each \code{q} for the user-specified sample coverages. \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data[-1,]
#' dij <-  FunDdata.inc$dij
#' out <- EstimateFD(data = data, distM = dij,datatype = "abundance")
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). functional diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of functional diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive functional diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
EstimateFD <- function(data, distM, datatype = "abundance", q = c(0,1,2), level = NULL, nboot = 50,conf = 0.95,threshold = NULL){
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,]
  }
  data <- data[rowSums(data)>0,,drop=FALSE]
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp <- rowMeans(sapply(dat, function(x) x/sum(x)))
    }else if(datatype=='incidence_freq'){
      tmp <- rowMeans(sapply(dat, function(x) x[-1]/sum(x[-1])))
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  }else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to (dmean+dmin)/2.",call. = FALSE)
  }
  
  if(is.null(level)){
    if(datatype=='abundance'){
      level <- sapply(dat,function(x){
        ni <- sum(x)
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
      
    }else if(datatype=='incidence_freq'){
      level <- sapply(dat,function(x){
        ni <- x[1]
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
    }
    level <- min(level)
  }
  
  out <- invChatFD(datalist = dat, dij = distM, q = q, datatype = datatype,
                   level=level, nboot = nboot, conf = conf, tau = threshold)
  out$qFD.LCL[out$qFD.LCL<0] <- 0
  if((nboot>1)==FALSE){
    out$qFD.LCL <- NA
    out$qFD.LCL <- NA
  }
  out
  
}


#' @useDynLib FunD
#' @importFrom Rcpp sourceCpp
NULL

