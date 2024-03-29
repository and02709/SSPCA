#' This function provides basic manhattan plots for the loadings of each eigenvector
#'   provided to the function.
#' @param Zlist this contains the eigenvectors
#' @param sepAnalysis tells whether the eigenvectors should be displayed separately
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords manhattan plot
#' @export
#' @examples ManPlot(Zlist, sepAnalysis, nresp, man.thresh, resp.names)

ManPlot <- function(Zlist, sepAnalysis, nresp, man.thresh, resp.names){
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
          dftemp <- data.frame(xvar=c(1:num.load),yvar=abs(Umat))
        }
        else{
          dftemp <- data.frame(xvar=c(1:num.load),yvar=abs(Umat[,j]))
        }
        
        fname <- paste0("/outputfolder/",resp.names, " Manhattan_Eigenvector ", j, ".png")
        if(j==1){
          plot.list <- list(fname)
        }
        else{
          plot.list <- append(plot.list, list(fname))
        }
        gg_plot_title <- paste0(resp.names, " Manhattan_Eigenvector ", j)
        
        if(man.thresh > 0){
          man_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
            geom_point(colour = "black") +
            geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
            #labs(size = "Loading \nAbsolute \nValue") + 
            ggtitle(gg_plot_title) +
            xlab("Loading Index") + 
            ylab("Loading Absolute Value") + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
        }
        else{
          man_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
            geom_point(colour = "black") +
            #geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
            #labs(size = "Loading \nAbsolute \nValue") + 
            ggtitle(gg_plot_title) +
            xlab("Loading Index") + 
            ylab("Loading Absolute Value") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
        }
        
        
        
        png(filename=fname)
        print(man_plot)
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
            dftemp <- data.frame(xvar=c(1:num.load),yvar=abs(Umat))
          }
          else{
            dftemp <- data.frame(xvar=c(1:num.load),yvar=abs(Umat[,j]))
          }
          
          fname <- paste0("/outputfolder/",resp.names[i], " Manhattan Eigenvector ", j, ".png")
          if(i==1 && j==1){
            plot.list <- list(fname)
          }
          else{
            plot.list <- append(plot.list, list(fname))
          }
          gg_plot_title <- paste0(resp.names[i], " Manhattan Eigenvector ", j)
          
          if(man.thresh > 0){
            man_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
              geom_point(colour = "black") +
              geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
              #labs(size = "Loading \nAbsolute \nValue") + 
              ggtitle(gg_plot_title) +
              xlab("Loading Index") + 
              ylab("Loading Absolute Value") + 
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
          }
          else{
            man_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
              geom_point(colour = "black") +
              #geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
              #labs(size = "Loading \nAbsolute \nValue") + 
              ggtitle(gg_plot_title) +
              xlab("Loading Index") + 
              ylab("Loading Absolute Value") + 
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
          }
          
          
          png(filename=fname)
          print(man_plot)
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
        dftemp <- data.frame(xvar=c(1:num.load),yvar=abs(Umat))
      }
      else{
        dftemp <- data.frame(xvar=c(1:num.load),yvar=abs(Umat[,j]))
      }
      
      fname <- paste0("/outputfolder/Manhattan_Eigenvector_", j, ".png")
      if(j==1){
        plot.list <- list(fname)
      }
      else{
        plot.list <- append(plot.list, list(fname))
      }
      gg_plot_title <- paste0("Manhattan_Eigenvector_", j)
      
      if(man.thresh > 0){
        man_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
          geom_point(colour = "black") +
          geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
          #labs(size = "Loading \nAbsolute \nValue") + 
          ggtitle(gg_plot_title) +
          xlab("Loading Index") + 
          ylab("Loading Absolute Value") + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
      }
      else{
        man_plot <- ggplot2::ggplot(data=dftemp, mapping=aes(x=xvar, y=yvar)) +
          geom_point(colour = "black") +
          #geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
          #labs(size = "Loading \nAbsolute \nValue") + 
          ggtitle(gg_plot_title) +
          xlab("Loading Index") + 
          ylab("Loading Absolute Value") + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
      }
      
      
      png(filename=fname)
      print(man_plot)
      dev.off()
      
    }
  }
  
  return(plot.list)
  
  
}