#' This function provides visualization for the responses.
#' @param ytrain training response dataset
#' @param ytest testing response dataset
#' @param sepAnalysis tells whether the eigenvectors should be displayed separately
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords response density plot
#' @export
#' @examples ydens(ytrain, ytest, nresp, resp.names)

ydens<- function(ytrain, ytest, nresp, resp.names){
  
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
    
    fname <- paste0(resp.names[i], " train density.png")
    if(i==1){
      plot.list <- list(fname)
    }
    else{
      plot.list <- append(plot.list, list(fname))
    }
    gg_plot_title <- paste0(resp.names[i], " train density")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=y, colour="black", fill="red")) +
      geom_density(alpha=.4, show.legend = F) + 
      ggtitle(gg_plot_title) +
      xlab(resp.names[i]) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
    dftemp <- data.frame(y=y.ts)
    
    fname <- paste0(resp.names[i], " test density.png")
    
    plot.list <- append(plot.list, list(fname))
    
    gg_plot_title <- paste0(resp.names[i], " test density")
    
    densplot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=y, colour="black", fill="red")) +
      geom_density(alpha=.4, show.legend = F) + 
      ggtitle(gg_plot_title) +
      xlab(resp.names[i]) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    png(filename=fname)
    print(densplot)
    dev.off()
    
  }
  return(plot.list)
  
}