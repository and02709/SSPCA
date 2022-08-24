#' This function provides visualization for the loadings of each eigenvector
#'   provided to the function.
#' @param Zlist this contains the eigenvectors
#' @param sepAnalysis tells whether the eigenvectors should be displayed separately
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords eigenvector plot
#' @export
#' @examples EvecPlot(Zlist, sepAnalysis, nresp, resp.names)

EvecPlot <- function(Zlist, sepAnalysis, nresp, resp.names){
  if(sepAnalysis){
    hold.pattern <- Zlist
    n.resp <- length(hold.pattern)
    if(n.resp!=nresp) warning("Discrepancy between number of responses and decomposition process")
    if(n.resp==1){
      hold.U <- hold.pattern[[1]]$U
      Umat <- as.matrix(hold.U)
      num.vec <- ncol(Umat)
      num.load <- nrow(Umat)
      
      for(j in 1:num.vec){
        if(num.vec==1){
          dftemp <- data.frame(xvar=c(1:num.load),yvar=Umat)
        }
        else{
          dftemp <- data.frame(xvar=c(1:num.load),yvar=Umat[,j])
        }
        
        fname <- paste0(resp.names, " Eigenvector ", j, ".png")
        if(j==1){
          plot.list <- list(fname)
        }
        else{
          plot.list <- append(plot.list, list(fname))
        }
        gg_plot_title <- paste0(resp.names, " Eigenvector ", j)
        
        reg_e_vec_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
          geom_point(colour = "black", show.legend = F) +
          geom_smooth(method="loess", span=.9, se=FALSE) +
          geom_hline(yintercept = 0, color = "red", linetype = "solid") +
          ggtitle(gg_plot_title) +
          xlab("Loading Index") + 
          ylab("Loading Value") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
        png(filename=fname)
        print(reg_e_vec_plot)
        dev.off()
        
      }
      
      
    }
    else{
      
      for(i in 1:n.resp){
        hold.U <- hold.pattern[[i]]$U
        Umat <- as.matrix(hold.U)
        num.vec <- ncol(Umat)
        num.load <- nrow(Umat)
        
        for(j in 1:num.vec){
          if(num.vec==1){
            dftemp <- data.frame(xvar=c(1:num.load),yvar=Umat)
          }
          else{
            dftemp <- data.frame(xvar=c(1:num.load),yvar=Umat[,j])
          }
          
          fname <- paste0(resp.names[i], " Eigenvector ", j, ".png")
          if(i==1 && j==1){
            plot.list <- list(fname)
          }
          else{
            plot.list <- append(plot.list, list(fname))
          }
          gg_plot_title <- paste0(resp.names[i], " Eigenvector ", j)
          
          reg_e_vec_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
            geom_point(colour = "black", show.legend = F) +
            geom_smooth(method="loess", span=.9, se=FALSE) +
            geom_hline(yintercept = 0, color = "red", linetype = "solid") +
            ggtitle(gg_plot_title) +
            xlab("Loading Index") + 
            ylab("Loading Value") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
          png(filename=fname)
          print(reg_e_vec_plot)
          dev.off()
          
        }
      }
      
    }
  }
  else{
    Umat <- as.matrix(Zlist$U)
    num.vec <- ncol(Umat)
    num.load <- nrow(Umat)
    for(j in 1:num.vec){
      if(num.vec==1){
        dftemp <- data.frame(xvar=c(1:num.load),yvar=Umat)
      }
      else{
        dftemp <- data.frame(xvar=c(1:num.load),yvar=Umat[,j])
      }
      
      fname <- paste0("Eigenvector_", j, ".png")
      if(j==1){
        plot.list <- list(fname)
      }
      else{
        plot.list <- append(plot.list, list(fname))
      }
      gg_plot_title <- paste0("Eigenvector_", j)
      
      reg_e_vec_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
        geom_point(colour = "black", show.legend = F) +
        geom_smooth(method="loess", span=.9, se=FALSE) +
        geom_hline(yintercept = 0, color = "red", linetype = "solid") +
        ggtitle(gg_plot_title) +
        xlab("Loading Index") + 
        ylab("Loading Value") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      png(filename=fname)
      print(reg_e_vec_plot)
      dev.off()
      
    }
  }
  
  return(plot.list)
  
  
}