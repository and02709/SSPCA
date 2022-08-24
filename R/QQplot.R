#' This function provides a QQ plot for the training and testing data
#' @param ytrain training response dataset
#' @param ytest testing response dataset
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords 
#' @export
#' @examples QQplot(ytrain, ytest, nresp, resp.names)

QQplot <- function(ytrain, ytest, nresp, resp.names){
  
  
  for(i in 1:nresp){
    if(nresp==1){
      y.tr <- ytrain
      y.ts <- ytest
    }
    else{
      y.tr <- ytrain[,i]
      y.ts <- ytest[,i]
    }
    
    dftemp <- data.frame(y=y.tr)
    
    fname <- paste0("QQ ", resp.names[i], " train.png")
    if(i==1){
      plot.list <- list(fname)
    }
    else{
      plot.list <- append(plot.list, list(fname))
    }
    gg_plot_title <- paste0("QQ ", resp.names[i], " train")
    
    densplot <- ggplot2::ggplot(data=dftemp, aes(sample = y)) +
      stat_qq(distribution = qnorm) +
      stat_qq_line(distribution = qnorm) +
      ggtitle(gg_plot_title) +
      xlab(paste0("QQ ", resp.names[i])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
    dftemp <- data.frame(y=y.ts)
    
    fname <- paste0("QQ ", resp.names[i], " test.png")
    
    plot.list <- append(plot.list, list(fname))
    
    gg_plot_title <- paste0("QQ ", resp.names[i], " test")
    
    densplot <- ggplot2::ggplot(data=dftemp, aes(sample = y)) +
      stat_qq(distribution = qnorm) +
      stat_qq_line(distribution = qnorm) +
      ggtitle(gg_plot_title) +
      xlab(paste0("QQ ", resp.names[i])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
  }
  
  return(plot.list)
  
}