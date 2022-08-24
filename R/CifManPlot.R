#' This function provides a manhattan plot unique to CIFTI data
#' @param Zlist contains the eigenevectors corresponding to the loadings
#' @param sepAnalysis tells whether the eigenvectors should be displayed separately
#' @param nresp number of response vectors present in datasets
#' @param man.thresh flag for whether this is a manhattan plot
#' @param resp.names gives list of names for each response
#' @param specific.groups allows certain groups to be selected
#' @param exclude.groups allows certain groups to be excluded
#' @param color scheme for plot
#' @keywords 
#' @export
#' @examples CifManPlot(Zlist, sepAnalysis, nresp, man.thresh, resp.names, 
#'   specific.groups=NULL, exclude.groups=NULL, special.colors=NULL)


CifManPlot <- function(Zlist, sepAnalysis, nresp, man.thresh, resp.names, 
                       specific.groups=NULL, exclude.groups=NULL, special.colors=NULL){
  if(sepAnalysis){
    hold.pattern <- Zlist
    n.resp <- length(hold.pattern)
    if(n.resp!=nresp) warning("Discrepancy between number of responses and decomposition process")
    if(n.resp==1){
      hold.U <- hold.pattern[[1]]$U
      Umat <- as.matrix(hold.U)
      num.vec <- ncol(Umat)
      num.load <- nrow(Umat)
      dftemp <- read.table("ROI.txt", header=T)
      dftemp$group <- as.factor(dftemp$group)
      
      if(!is.null(special.colors)){
        dfcolor <- data.frame(special.colors)
        colnames(dfcolor) <- c("group", "color")
        dfcolor$group <- as.factor(dfcolor$group)
        #dftemp <- inner_join(dftemp, dfcolor, by="group")
      }
      
      for(j in 1:num.vec){
        if(num.vec==1){
          dftemp$U <- abs(Umat)
        }
        else{
          dftemp$U <- abs(Umat[,j])
        }
        
        # not convinced these variables are needed for anything
        #dgroup <- levels(dftemp$group)
        #dftemp$plot.index <- seq(1:length(dftemp$index))
        
        if(!is.null(specific.groups) || !is.null(exclude.groups)){
          if(!is.null(specific.groups)){
            dftemp <- dftemp %>% dplyr::filter(group %in% specific.groups)
          }
          if(!is.null(exclude.groups)){
            dftemp <- dftemp %>% dplyr::filter(!group %in% exclude.groups)
          }
        }
        
        if(!is.null(special.colors)){
          if(!is.null(specific.groups)){
            dfcolor <- dfcolor %>% dplyr::filter(group %in% specific.groups)
          }
          if(!is.null(exclude.groups)){
            dfcolor <- dfcolor %>% dplyr::filter(!group %in% exclude.groups)
          }
        }
        
        # next code :
        # make groups by variable "group"
        # within a group, order by index value.
        # make a new variable sequence 1 to n within group called "newindx"
        # order dataframe by group then newindx within group
        # save in dftemp
        dftemp<- dftemp %>% dplyr::group_by(group) %>%
          dplyr::arrange(index) %>%
          dplyr::mutate(newindx = 1:dplyr::n()) %>%
          dplyr::arrange(group,newindx) 
        
        # now need to have x axis sequencing change smoothly across group levels
        # compute max value of newindx within each group
        # cumsum() makes a running total of the max_idx values that came prior
        # make new variable idx_add which contains the value of running index that
        # was at the end of the previous group - - this makes a value of start index (idx_add)
        # for this group's row in this dframe.
        dftemp2 <-  dftemp %>% dplyr::group_by(group) %>%
          dplyr::summarise(max_idx = max(newindx)) %>%
          dplyr::mutate(idx_add = dplyr::lag(cumsum(max_idx), default = 0)) %>% 
          dplyr::select(group, idx_add) %>%
          dplyr::ungroup() 
        # join the start index dframe (dftemp2) with the plot dataframe dftemp
        # make a new variable idx_cum that creates the running index position for each
        # observation. The plot_data dataframe will be used in ggplot.
        plot_data <- dftemp %>% dplyr::inner_join(dftemp2, by="group") %>%
          dplyr::mutate(idx_cum = newindx + idx_add)
        # find the middle index position to put the group index label on x axis:
        
        
        
        axis_set <- plot_data %>% 
          dplyr::group_by(group) %>% 
          dplyr::summarize(center = mean(idx_cum))
        # next code makes a column with upper limit for y plotting 
        ylim <- plot_data %>% 
          dplyr::mutate(ylim = (max(U) + 0.2*max(U))) %>% 
          dplyr::pull(ylim)
        
        fname <- paste0(resp.names, " Manhattan Cifti Eigenvector ", j, ".png")
        if(j==1){
          plot.list <- list(fname)
        }
        else{
          plot.list <- append(plot.list, list(fname))
        }
        gg_plot_title <- paste0(resp.names, " Manhattan Cifti Eigenvector ", j)
        
        if(!is.null(special.colors)){
          dfcolor2 <- dfcolor[order(as.numeric(dfcolor$group)),]
        }
        else{
          dfcolor2 <- data.frame(color=rep(c("#276FBF", "#183059"), unique(length(axis_set$group))))
          
        }
        
        
        if(man.thresh > 0){
          manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = U, 
                                                                  color = as_factor(group))) +
            #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
            geom_point(alpha = 0.75, size=.85) +
            geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
            scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
            scale_y_continuous(expand = c(0,0), limits = c(0, ylim[1])) +
            #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
            scale_color_manual(values = dfcolor2$color) +
            scale_size_continuous(range = c(0.25,2.5)) +
            labs(x = NULL, 
                 y = "Absolute Value \nEigenvector \nLoadings") + 
            ggtitle(gg_plot_title) +
            theme_minimal() + 
            theme( 
              legend.position = "none",
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title.y = ggtext::element_markdown(),
              axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
            )
        }
        else{
          manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = U, 
                                                                  color = as_factor(group))) +
            #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
            geom_point(alpha = 0.75, size=.85) +
            scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
            scale_y_continuous(expand = c(0,0), limits = c(0, ylim[1])) +
            #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
            scale_color_manual(values = dfcolor2$color) +
            scale_size_continuous(range = c(0.25,2.5)) +
            labs(x = NULL, 
                 y = "Absolute Value \nEigenvector \nLoadings") + 
            ggtitle(gg_plot_title) +
            theme_minimal() + 
            theme( 
              legend.position = "none",
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title.y = ggtext::element_markdown(),
              axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
            )
        }
        
        
        png(filename=fname)
        print(manhplot)
        dev.off()
        
        
      }
    }
    
    else{
      
      for(i in 1:n.resp){
        hold.U <- hold.pattern[[i]]$U
        Umat <- as.matrix(hold.U)
        num.vec <- ncol(Umat)
        num.load <- nrow(Umat)
        dftemp <- read.table("ROI.txt", header=T)
        dftemp$group <- as.factor(dftemp$group)
        
        if(!is.null(special.colors)){
          dfcolor <- data.frame(special.colors)
          colnames(dfcolor) <- c("group", "color")
          dfcolor$group <- as.factor(dfcolor$group)
          #dftemp <- inner_join(dftemp, dfcolor, by="group")
        }
        
        for(j in 1:num.vec){
          if(num.vec==1){
            dftemp$U <- abs(Umat)
          }
          else{
            dftemp$U <- abs(Umat[,j])
          }
          
          # not convinced these variables are needed for anything
          #dgroup <- levels(dftemp$group)
          #dftemp$plot.index <- seq(1:length(dftemp$index))
          
          if(!is.null(specific.groups) || !is.null(exclude.groups)){
            if(!is.null(specific.groups)){
              dftemp <- dftemp %>% dplyr::filter(group %in% specific.groups)
            }
            if(!is.null(exclude.groups)){
              dftemp <- dftemp %>% dplyr::filter(!group %in% exclude.groups)
            }
          }
          
          if(!is.null(special.colors)){
            if(!is.null(specific.groups)){
              dfcolor <- dfcolor %>% dplyr::filter(group %in% specific.groups)
            }
            if(!is.null(exclude.groups)){
              dfcolor <- dfcolor %>% dplyr::filter(!group %in% exclude.groups)
            }
          }
          
          # next code :
          # make groups by variable "group"
          # within a group, order by index value.
          # make a new variable sequence 1 to n within group called "newindx"
          # order dataframe by group then newindx within group
          # save in dftemp
          dftemp<- dftemp %>% dplyr::group_by(group) %>%
            dplyr::arrange(index) %>%
            dplyr::mutate(newindx = 1:dplyr::n()) %>%
            dplyr::arrange(group,newindx) 
          
          # now need to have x axis sequencing change smoothly across group levels
          # compute max value of newindx within each group
          # cumsum() makes a running total of the max_idx values that came prior
          # make new variable idx_add which contains the value of running index that
          # was at the end of the previous group - - this makes a value of start index (idx_add)
          # for this group's row in this dframe.
          dftemp2 <-  dftemp %>% dplyr::group_by(group) %>%
            dplyr::summarise(max_idx = max(newindx)) %>%
            dplyr::mutate(idx_add = dplyr::lag(cumsum(max_idx), default = 0)) %>% 
            dplyr::select(group, idx_add) %>%
            dplyr::ungroup() 
          # join the start index dframe (dftemp2) with the plot dataframe dftemp
          # make a new variable idx_cum that creates the running index position for each
          # observation. The plot_data dataframe will be used in ggplot.
          plot_data <- dftemp %>% dplyr::inner_join(dftemp2, by="group") %>%
            dplyr::mutate(idx_cum = newindx + idx_add)
          # find the middle index position to put the group index label on x axis:
          
          
          
          
          axis_set <- plot_data %>% 
            dplyr::group_by(group) %>% 
            dplyr::summarize(center = mean(idx_cum))
          # next code makes a column with upper limit for y plotting 
          ylim <- plot_data %>% 
            dplyr::mutate(ylim = (max(U) + 0.2*max(U))) %>% 
            dplyr::pull(ylim)
          
          fname <- paste0(resp.names[i], " Manhattan Cifti Eigenvector ", j, ".png")
          if(i==1 && j==1){
            plot.list <- list(fname)
          }
          else{
            plot.list <- append(plot.list, list(fname))
          }
          gg_plot_title <- paste0(resp.names[i], " Manhattan Cifti Eigenvector ", j)
          
          if(!is.null(special.colors)){
            dfcolor2 <- dfcolor[order(as.numeric(dfcolor$group)),]
          }
          else{
            dfcolor2 <- data.frame(color=rep(c("#276FBF", "#183059"), unique(length(axis_set$group))))
            
          }
          
          if(man.thresh > 0){
            manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = U, 
                                                                    color = as_factor(group))) +
              #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
              geom_point(alpha = 0.75, size=.85) +
              geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
              scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
              scale_y_continuous(expand = c(0,0), limits = c(0, ylim[1])) +
              #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
              scale_color_manual(values = dfcolor2$color) +
              scale_size_continuous(range = c(0.25,2.5)) +
              labs(x = NULL, 
                   y = "Absolute Value \nEigenvector \nLoadings") + 
              ggtitle(gg_plot_title) +
              theme_minimal() + 
              theme( 
                legend.position = "none",
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                axis.title.y = ggtext::element_markdown(),
                axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
              )
          }
          else{
            manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = U, 
                                                                    color = as_factor(group))) +
              #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
              geom_point(alpha = 0.75, size=.85) +
              scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
              scale_y_continuous(expand = c(0,0), limits = c(0, ylim[1])) +
              #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
              scale_color_manual(values = dfcolor2$color) +
              scale_size_continuous(range = c(0.25,2.5)) +
              labs(x = NULL, 
                   y = "Absolute Value \nEigenvector \nLoadings") + 
              ggtitle(gg_plot_title) +
              theme_minimal() + 
              theme( 
                legend.position = "none",
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                axis.title.y = ggtext::element_markdown(),
                axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
              )
          }
          
          
          png(filename=fname)
          print(manhplot)
          dev.off()
          
          
        }
      }
    }
  }
  
  
  else{
    Umat <- as.matrix(Zlist$U)
    num.vec <- ncol(Umat)
    num.load <- nrow(Umat)
    
    #dftemp <- read.table("ROI.txt", header=T)
    #dftemp$group <- as.factor(dftemp$group)
    
    
    
    # use existing mapping file.
    dftemp <- read.table("ROI.txt", header=TRUE)
    dftemp$group <- as.factor(dftemp$group)
    
    if(!is.null(special.colors)){
      dfcolor <- data.frame(special.colors)
      colnames(dfcolor) <- c("group", "color")
      dfcolor$group <- as.factor(dfcolor$group)
      #dftemp <- inner_join(dftemp, dfcolor, by="group")
    }
    
    # add y axis values to dataframe
    
    for(j in 1:num.vec){
      if(num.vec==1){
        dftemp$U <- abs(Umat)
      }
      else{
        dftemp$U <- abs(Umat[,j])
      }
      
      # not convinced these variables are needed for anything
      #dgroup <- levels(dftemp$group)
      #dftemp$plot.index <- seq(1:length(dftemp$index))
      
      # next code :
      # make groups by variable "group"
      # within a group, order by index value.
      # make a new variable sequence 1 to n within group called "newindx"
      # order dataframe by group then newindx within group
      # save in dftemp
      
      if(!is.null(specific.groups) || !is.null(exclude.groups)){
        if(!is.null(specific.groups)){
          dftemp <- dftemp %>% dplyr::filter(group %in% specific.groups)
        }
        if(!is.null(exclude.groups)){
          dftemp <- dftemp %>% dplyr::filter(!group %in% exclude.groups)
        }
      }
      
      if(!is.null(special.colors)){
        if(!is.null(specific.groups)){
          dfcolor <- dfcolor %>% dplyr::filter(group %in% specific.groups)
        }
        if(!is.null(exclude.groups)){
          dfcolor <- dfcolor %>% dplyr::filter(!group %in% exclude.groups)
        }
      }
      
      dftemp<- dftemp %>% dplyr::group_by(group) %>%
        dplyr::arrange(index) %>%
        dplyr::mutate(newindx = 1:dplyr::n()) %>%
        dplyr::arrange(group,newindx) 
      
      
      
      # now need to have x axis sequencing change smoothly across group levels
      # compute max value of newindx within each group
      # cumsum() makes a running total of the max_idx values that came prior
      # make new variable idx_add which contains the value of running index that
      # was at the end of the previous group - - this makes a value of start index (idx_add)
      # for this group's row in this dframe.
      dftemp2 <-  dftemp %>% dplyr::group_by(group) %>%
        dplyr::summarise(max_idx = max(newindx)) %>%
        dplyr::mutate(idx_add = dplyr::lag(cumsum(max_idx), default = 0)) %>% 
        dplyr::select(group, idx_add) %>%
        dplyr::ungroup() 
      # join the start index dframe (dftemp2) with the plot dataframe dftemp
      # make a new variable idx_cum that creates the running index position for each
      # observation. The plot_data dataframe will be used in ggplot.
      plot_data <- dftemp %>% dplyr::inner_join(dftemp2, by="group") %>%
        dplyr::mutate(idx_cum = newindx + idx_add)
      # find the middle index position to put the group index label on x axis:
      
      # Filter for connection groups provided by user
      
      
      
      axis_set <- plot_data %>% 
        dplyr::group_by(group) %>% 
        dplyr::summarize(center = mean(idx_cum))
      # next code makes a column with upper limit for y plotting 
      ylim <- plot_data %>% 
        dplyr::mutate(ylim = (max(U) + 0.2*max(U))) %>% 
        dplyr::pull(ylim)
      
      fname <- paste0("Manhattan Cifti Eigenvector ", j, ".png")
      if(j==1){
        plot.list <- list(fname)
      }
      else{
        plot.list <- append(plot.list, list(fname))
      }
      gg_plot_title <- paste0("Manhattan Cifti Eigenvector ", j)
      
      if(!is.null(special.colors)){
        dfcolor2 <- dfcolor[order(as.numeric(dfcolor$group)),]
      }
      else{
        dfcolor2 <- data.frame(color=rep(c("#276FBF", "#183059"), unique(length(axis_set$group))))
        
      }
      
      if(man.thresh > 0){
        manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = U, 
                                                                color = as_factor(group))) +
          #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
          geom_point(alpha = 0.75, size=.85) +
          geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
          scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
          scale_y_continuous(expand = c(0,0), limits = c(0, ylim[1])) +
          #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
          scale_color_manual(values = dfcolor2$color) +
          scale_size_continuous(range = c(0.25,2.5)) +
          labs(x = NULL, 
               y = "Absolute Value \nEigenvector \nLoadings") + 
          ggtitle(gg_plot_title) +
          theme_minimal() + 
          theme( 
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.y = ggtext::element_markdown(),
            axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
          )
      }
      else{
        manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = U, 
                                                                color = as_factor(group))) +
          #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
          geom_point(alpha = 0.75, size=.85) +
          scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
          scale_y_continuous(expand = c(0,0), limits = c(0, ylim[1])) +
          #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
          scale_color_manual(values = dfcolor2$color) +
          scale_size_continuous(range = c(0.25,2.5)) +
          labs(x = NULL, 
               y = "Absolute Value \nEigenvector \nLoadings") + 
          ggtitle(gg_plot_title) +
          theme_minimal() + 
          theme( 
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.y = ggtext::element_markdown(),
            axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
          )
      }
      
      
      png(filename=fname)
      print(manhplot)
      dev.off()
      
      
      
      
    }
  }
  
  return(plot.list)
  
  
}