#' This function provides a comparison for the training and testing data
#'   with the fitted values
#' @param ytrain training response dataset
#' @param ytest testing response dataset
#' @param yhat this is the fitted data to be compared to the training and testing
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords observed predicted plot
#' @export
#' @examples QQplot(ytrain, ytest, nresp, resp.names)

OPplot <- function(ytrain, ytest, yhat, nresp, resp.names){
  
  y.hat.train <- yhat[[1]]
  y.hat.test <- yhat[[2]]
  
  for(i in 1:nresp){
    if(nresp==1){
      y.tr <- ytrain
      y.ts <- ytest
      y.hat.tr <- y.hat.train
      y.hat.ts <- y.hat.test
    }
    else{
      y.tr <- ytrain[,i]
      y.ts <- ytest[,i]
      y.hat.tr <- y.hat.train[,i]
      y.hat.ts <- y.hat.test[,i]
    }
    
    dftemp <- data.frame(y=y.tr, yhat=y.hat.tr)
    
    fname <- paste0("/outputfolder/Observed vs Predicted ", resp.names[i], " train.png")
    if(i==1){
      plot.list <- list(fname)
    }
    else{
      plot.list <- append(plot.list, list(fname))
    }
    gg_plot_title <- paste0("Observed vs Predicted ", resp.names[i], " train")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=yhat, y=y)) +
      geom_point(colour = "black", size = 1.2, show.legend = F) +
      geom_smooth(method="loess", span=.9, se=TRUE) +
      #geom_hline(yintercept = 0, color = "red", linetype = "solid") +
      ggtitle(gg_plot_title) +
      xlab("Fitted Value") + 
      ylab("Observed Value") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
    dftemp <- data.frame(y=y.ts, yhat=y.hat.ts)
    
    fname <- paste0("/outputfolder/Observed vs Predicted ", resp.names[i], " test.png")
    
    plot.list <- append(plot.list, list(fname))
    
    gg_plot_title <- paste0("Observed vs Predicted ", resp.names[i], " test")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=yhat, y=y)) +
      geom_point(colour = "black", size = 1.2, show.legend = F) +
      geom_smooth(method="loess", span=.9, se=TRUE) +
      #geom_hline(yintercept = 0, color = "red", linetype = "solid") +
      ggtitle(gg_plot_title) +
      xlab("Fitted Value") + 
      ylab("Observed Value") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
  }
  
  return(plot.list)
  
}