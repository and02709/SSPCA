# This function encompasses all possible options used for conducting SPCA and
#   sparse SPCA.  This includes generation of the eigenvectors, encoding of
#   of the data, performance metrics, and several visualization options



OmniSPCA <- function(xtrain, ytrain, xtest, ytest, resp.names = NULL, ycont=T, ybinary=F, sepAnalysis = F, 
                     K=1, sumabsv=0, strictEV=T, Sens.analysis=F,
                     sepclass=T, binClassCombine=F, NoExtraEV = F, ydensity=F,
                     dpyhat=F, dsq=F, qqplot=F, opplot=F, 
                     Evec.plot = F, manhattan.plot=F, cifti.man.plot=F, 
                     man.thresh = 0, niter=20, orth=TRUE, 
                     trace=TRUE, v=NULL, center=TRUE, test_metrics=T,
                     cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE){
  
  if(is_tibble(ytrain) || is_tibble(ytest)){
    ytrain <- as.matrix(ytrain)
    ytest <- as.matrix(ytest)
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
    xtrain <- xtrain-xmeans
    xtest <- xtest-xmeans
  }
  
  nobs <- nrow(xtrain)
  npred <- ncol(xtrain)
  nresp <- ncol(ytrain)
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
    Cifti.manhattan.plot <- CifManPlot(Zlist=Z, sepAnalysis=sepAnalysis, nresp=nresp, man.thresh=man.thresh, resp.names=resp.names)
  }
  else{
    Cifti.manhattan.plot <- NULL
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
  
  
  
  return(list(Z, y_hat, TrTs, plot.list, man.plot.list, Cifti.manhattan.plot, y.density.list, yhat.density.list, ysq.density.list, QQ.list))  
  
}