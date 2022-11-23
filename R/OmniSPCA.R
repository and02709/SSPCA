#' This function provides all options for analysis using SSPCA
#' @param xtrain training predictor dataset
#' @param ytrain training response dataset
#' @param xtest testing predictor dataset
#' @param ytest testing response dataset
#' @param resp.names parameter containing response names
#' @param ycont flag indicates response data is continuous
#' @param ybinary flags for whether response is binary
#' @param sepAnalysis flag for whether analysis should be performed separately for
#'   each response
#' @param K number of eigenevectors allowed to be calculated
#' @param sumabsv causes shrinkage to occur to the loadings
#' @param strictEV flag whether to impose strict requiresments to accept number
#'   of eigenevectors or to select a more appropriate number
#' @param Sens.analysis not currently supported
#' @param not currently supported
#' @param binClassCombine not currently support
#' @param NoExtraEV unknown function
#' @param ydensity calculates density plot
#' @param dpyhat calculates y hat plot
#' @param dsq calculates squared error plot
#' @param qqplot generates qqplot
#' @param opplot calculates observed vs predicted plot
#' @param Evec.plot generates eigenevector plot
#' @param manhattan.plot generates basic manhattan plot
#' @param cifti.man.plot generates cifti specific manhattan plot
#' @param vol.man.plot generates brain volume specific manhattan plot
#' @param man.thresh visual horizontal bar for manhattan plot
#' @param nresp number of response vectors present in datasets
#' @param niter number of iterations
#' @param orth determines whether orthogonal vectors are required
#' @param trace this allows debugging
#' @param v focuses on the loadings vector
#' @param center this centers the predictor columns
#' @param test_metrics flag for calculating test metrics
#' @param cnames this gives column names
#' @param vpos this determines whether the matrix is to be positive only
#' @param vneg this determines whether the matrix is to be negative only
#' @param compute.pve computes cv using missing values
#' @keywords Omni
#' @export
#' @examples OmniSPCA <- function(xtrain, ytrain, xtest, ytest, resp.names = NULL, 
#'   ycont=T, ybinary=F, sepAnalysis = F, 
#'   K=1, sumabsv=0, strictEV=T, Sens.analysis=F,
#'   sepclass=T, binClassCombine=F, NoExtraEV = F, ydensity=F,
#'   dpyhat=F, dsq=F, qqplot=F, opplot=F, 
#'   Evec.plot = F, manhattan.plot=F, cifti.man.plot=F, 
#'   man.thresh = 0, niter=20, orth=TRUE, 
#'   trace=TRUE, v=NULL, center=TRUE, test_metrics=T,
#'   cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE)


OmniSPCA <- function(xtrain, ytrain, xtest, ytest, resp.names = NULL, ycont=T, ybinary=F, sepAnalysis = F, 
                     K=1, sumabsv=0, strictEV=T, Sens.analysis=F,
                     sepclass=T, binClassCombine=F, NoExtraEV = F, ydensity=F,
                     dpyhat=F, dsq=F, qqplot=F, opplot=F, 
                     Evec.plot = F, manhattan.plot=F, cifti.man.plot=F, 
                     vol.man.plot=F,
                     man.thresh = 0, specific.groups=NULL, exclude.groups=NULL, 
                     special.colors=NULL, 
                     niter=20, orth=TRUE, 
                     trace=TRUE, v=NULL, center=TRUE, test_metrics=T,
                     cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE){
  
  
  ytrain <- as.matrix(ytrain)
  ytest <- as.matrix(ytest)
  
  if(dim(as.matrix(xtrain))[[1]]==1 && dim(as.matrix(xtest))[[1]]==1){
    xtrain <- t(xtrain)
    xtest <- t(xtest)
  }
  
  
  
  
  nobs <- nrow(ytrain)
  npred <- ncol(xtrain)
  nresp <- ncol(ytrain)
  
  if(!is.null(resp.names) && length(resp.names)!=nresp) stop("Incorrect length of response names to response matrix provided")
  
  if(!is.null(resp.names)){
    colnames(ytrain) <- resp.names
    colnames(ytest) <- resp.names
  }
  else{
    resp.names <- paste0("y", 1:nresp)
  }
  
  if(center){
    xmeans <- colMeans(xtrain)
    xtrain <- t(apply(xtrain, 1, function(x) x-xmeans))
    xtest <- t(apply(xtest, 1, function(x) x-xmeans))
  }
  
  if(dim(as.matrix(xtrain))[[1]]==1 && dim(as.matrix(xtest))[[1]]==1){
    xtrain <- t(xtrain)
    xtest <- t(xtest)
  }
  
  Dual <- ifelse(npred>nobs, TRUE, FALSE)
  
  if(sepAnalysis){
    if(sumabsv==0){
      if(Dual){
        Z <- LoopDualSPC(ytrain=ytrain, xtrain=xtrain, ytest=ytest, xtest=xtest, 
                         ycont=ycont, nresp=nresp, strictEV=strictEV, nTopEvecs=K)
      }
      else{
        Z <- LoopSPC(ytrain=ytrain, xtrain=xtrain, ytest=ytest, xtest=xtest, 
                     ycont=ycont, nresp=nresp, strictEV=strictEV, nTopEvecs=K)
      }
    }
    else{
      Z <- LoopSSPC(xtrain=xtrain, ytrain=ytrain, xtest=xtest, ycont=ycont, 
                    nresp=nresp, sumabsv=sumabsv, niter=niter, K=K, orth=orth, 
                    trace=trace, v=v, center=center, cnames=cnames, vpos=vpos, 
                    vneg=vneg, compute.pve=compute.pve, strictEV=strictEV)
    }
  }
  else{
    if(sumabsv==0){
      if(Dual){
        Z <- DualSPC(ytrain=ytrain, xtrain=xtrain, ytest=ytest, xtest=xtest, 
                     ycont=ycont, ybinary=ybinary, nresp=nresp, strictEV=strictEV, 
                     nTopEvecs=K)
      }
      else{
        Z <- SPC(ytrain=ytrain, xtrain=xtrain, ytest=ytest, xtest=xtest, 
                 ycont=ycont, ybinary=ybinary, nresp=nresp, strictEV=strictEV, nTopEvecs=K)
      }
    }
    else{
      Z <- SSPC(xtrain=xtrain, ytrain=ytrain, xtest=xtest, ycont=ycont, ybinary=ybinary, 
                nresp=nresp, sumabsv=sumabsv, 
                niter=niter, K=K, orth=orth, trace=trace, v=v, center=center, 
                cnames=cnames, vpos=vpos, vneg=vneg, compute.pve=compute.pve, 
                strictEV=strictEV)
    }
  }
  
  if(sepAnalysis){
    if(ycont){
      y_hat <- LMcontSep(ytrain=ytrain, ytest=ytest, Z=Z, nresp=nresp, resp.names=resp.names)
      
    }
    if(!ycont && ybinary){
      y_hat <- GLMbinSep(ytrain=ytrain, ytest=ytest, Z=Z, nresp=nresp, resp.names=resp.names)
    }
    if(!ycont && !ybinary){
      y_hat <- GLMmultiSep(ytrain=ytrain, ytest=ytest, Z=Z, nresp=nresp, resp.names=resp.names)
    }
    
    
  }
  else{
    if(ycont){
      y_hat <- LMcont(ytrain=ytrain, ytest=ytest, ztrain=Z$Z, ztest=Z$z, nresp=nresp, resp.names=resp.names)
      
    }
    if(!ycont && ybinary){
      y_hat <- GLMbin(ytrain=ytrain, ztrain=Z$Z, ytest=ytest, ztest=Z$z, nresp=nresp, resp.names=resp.names)
    }
    if(!ycont && !ybinary){
      y_hat <- GLMmulti(ytrain=ytrain, ztrain=Z$Z, ytest=ytest, ztest=Z$z, nresp=nresp, resp.names=resp.names)
    }
  }
  
  if(test_metrics){
    if(ycont){
      TrTs <- continuous_metrics(ytrain=ytrain, ytest=ytest, y_hat=y_hat, nresp=nresp, resp.names=resp.names)
    }
    else{
      TrTs <- categorical_metrics(ytrain=ytrain, ytest=ytest, y_hat=y_hat, nresp=nresp, resp.names=resp.names)
    }
  }
  else{
    TrTs <- NULL
  }
  
  if(Evec.plot){
    plot.list <- EvecPlot(Zlist=Z, sepAnalysis=sepAnalysis, nresp=nresp, resp.names=resp.names)
  }
  else{
    plot.list <- NULL
  }
  
  if(manhattan.plot){
    man.plot.list <- ManPlot(Zlist=Z, sepAnalysis=sepAnalysis, nresp=nresp, man.thresh=man.thresh, resp.names=resp.names)
  }
  else{
    man.plot.list <- NULL
  }
  if(cifti.man.plot){
    if(npred != 61776) stop("Wrong number of predictors for a CIFTI format")
    Cifti.manhattan.plot <- CifManPlot(Zlist=Z, sepAnalysis=sepAnalysis, nresp=nresp, man.thresh=man.thresh, 
                                       resp.names=resp.names, specific.groups=specific.groups, exclude.groups=exclude.groups, 
                                       special.colors=special.colors)
  }
  else{
    Cifti.manhattan.plot <- NULL
  }
  if(vol.man.plot){
    if(npred != 379) stop("Wrong number of predictors for a volume format")
    Vol.manhattan.plot <- VManPlot(Zlist=Z, sepAnalysis=sepAnalysis, nresp=nresp, man.thresh=man.thresh, 
                                       resp.names=resp.names, specific.groups=specific.groups, exclude.groups=exclude.groups, 
                                       special.colors=special.colors)
  }
  else{
    Vol.manhattan.plot <- NULL
  }
  if(ydensity && ycont){
    y.density.list <- ydens(ytrain=ytrain, ytest=ytest, nresp=nresp, resp.names=resp.names)
  }
  else{
    y.density.list <- NULL
  }
  if(dpyhat && ycont){
    yhat.density.list <- densplotyhat(y_hat=y_hat, nresp=nresp, resp.names=resp.names)
  }
  else{
    yhat.density.list <- NULL
  }
  if(dsq && ycont){
    ysq.density.list <- densplotSq(TrTs = TrTs, nresp=nresp, resp.names=resp.names)
  }
  else{
    ysq.density.list <- NULL
  }
  if(qqplot && ycont){
    QQ.list <- QQplot(ytrain=ytrain, ytest=ytest, nresp=nresp, resp.names=resp.names)
  }
  else{
    QQ.list <- NULL
  }
  if(opplot && ycont){
    OP.list <- OPplot(ytrain=ytrain, ytest=ytest, yhat=y_hat, nresp=nresp, resp.names=resp.names)
  }
  else{
    OP.list <- NULL
  }
  
  if(sepAnalysis){
    Z_restructured <- restructuring.func(Z=Z, resp.names=resp.names, nresp=nresp)
    Z <- Z_restructured
  }
  
  
  
  return(list(Z.U.lambda=Z, y.hat=y_hat, metrics=TrTs, plot.list=plot.list, man.plot.list=man.plot.list, 
              Cifti.manhattan.plot=Cifti.manhattan.plot, Vol.manhattan.plot=Vol.manhattan.plot,
              y.density.list=y.density.list, 
              yhat.density.list=yhat.density.list, ysq.density.list=ysq.density.list, 
              QQ.list=QQ.list))  
  
}