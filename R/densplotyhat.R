#' This function provides visualization for the responses
#' @param y_hat this contains the fitted/predicted responses
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords 
#' @export
#' @examples densplotyhat(y_hat, nresp, resp.names)

densplotyhat<- function(y_hat, nresp, resp.names){
  y.hat.train <- y_hat[[1]]
  y.hat.test <- y_hat[[2]]
  
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
    
    fname <- paste0("Fitted ", resp.names[i], " train density.png")
    if(i==1){
      plot.list <- list(fname)
    }
    else{
      plot.list <- append(plot.list, list(fname))
    }
    gg_plot_title <- paste0("Fitted ", resp.names[i], " train density")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=yhat, colour="black", fill="red")) +
      geom_density(alpha=.4, show.legend = F) + 
      ggtitle(gg_plot_title) +
      xlab(paste0("Fitted ", resp.names[i])) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
    dftemp <- data.frame(yhat=y.ts)
    
    fname <- paste0("Fitted ", resp.names[i], " test density.png")
    
    plot.list <- append(plot.list, list(fname))
    
    gg_plot_title <- paste0("Fitted ", resp.names[i], " test density")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=yhat, colour="black", fill="red")) +
      geom_density(alpha=.4, show.legend = F) + 
      ggtitle(gg_plot_title) +
      xlab(paste0("Fitted ", resp.names[i])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
  }
  
  return(plot.list)
  
}