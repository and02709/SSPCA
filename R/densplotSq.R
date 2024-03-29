#' This function provides visualization for the responses
#' @param TrTs this contains the fitted/predicted responses
#' @param sepAnalysis tells whether the eigenvectors should be displayed separately
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords fitted square distance density plot
#' @export
#' @examples densplotSq(TrTs, nresp, resp.names)

densplotSq <- function(TrTs, nresp, resp.names){
  y.hat.train <- TrTs[[1]]
  y.hat.test <- TrTs[[2]]
  
  for(i in 1:nresp){
    if(nresp==1){
      y.tr <- y.hat.train
      y.ts <- y.hat.test
    }
    else{
      y.tr <- y.hat.train[,i]
      y.ts <- y.hat.test[,i]
    }
    
    dftemp <- data.frame(yhat=y.tr)
    
    fname <- paste0("/outputfolder/Sqare Error ", resp.names[i], " train density.png")
    if(i==1){
      plot.list <- list(fname)
    }
    else{
      plot.list <- append(plot.list, list(fname))
    }
    gg_plot_title <- paste0("Square Error ", resp.names[i], " train density")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=yhat, colour="black", fill="red")) +
      geom_density(alpha=.4, show.legend = F) + 
      ggtitle(gg_plot_title) +
      xlab(paste0("Square Error ", resp.names[i])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
    dftemp <- data.frame(yhat=y.ts)
    
    fname <- paste0("/outputfolder/Sqare Error ", resp.names[i], " test density.png")
    
    plot.list <- append(plot.list, list(fname))
    
    gg_plot_title <- paste0("Square Error ", resp.names[i], " test density")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=yhat, colour="black", fill="red")) +
      geom_density(alpha=.4, show.legend = F) + 
      ggtitle(gg_plot_title) +
      xlab(paste0("Square Error ", resp.names[i])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
  }
  
  return(plot.list)
  
}