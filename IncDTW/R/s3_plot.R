
plot.rundtw <- function(x, knn = TRUE, minima = TRUE, 
                        scale = c("none", "01", "z"),  
                        selDim = 1, lix = 1, Q = NULL, C = NULL, normalize = c("none", "01", "z"), ...){

   if (!missing("normalize")){
      warning("Argument 'normalize' is deprecated, use 'scale' instead. 
              The parameter 'scale' is set equal the parameter 'normalize'.")
      scale <- normalize
   }
   
   if(is.list(x$dist)){
      
      x$dist <- x$dist[[lix]]
      if(!is.null(x$knn_indices)){
         x$knn_values  <- x$knn_values [x$knn_list_indices == lix]
         x$knn_indices <- x$knn_indices[x$knn_list_indices == lix]
         if(length(x$knn_values) == 0){
            x$knn_values  <- NULL
            x$knn_indices <- NULL
         }
      }
      if(!is.null(x$C)) x$C <- x$C[[lix]]
      
      return(
         plot(x, knn = knn, minima = minima, scale = scale, 
              selDim = selDim, lix = 1, Q = Q, C = C)
      )
      
   }
   if(!is.null(Q) & !is.null(x$Q)){
      if(!identical(Q, x$Q)){
         warning("Q is set and x$Q is not NULL, the set Q is plotted and x$Q ignored")
      }
   }
   if(!is.null(C) & !is.null(x$C)){
      if(!identical(C, x$C)){
         warning("C is set and x$C is not NULL, the set C is plotted and x$C ignored")
      }
   }
   if(is.null(Q) & !is.null(x$Q)) Q <- x$Q
   if(is.null(C) & !is.null(x$C)) C <- x$C
   group <- fct <- x_time <- y_value <- categ <- NULL
   
   scale <- match.arg(scale, c("none", "01", "z"))
   
   dfp1 <- data.frame(x_time = 1:length(x$dist),
                      y_value = x$dist,
                      categ = "dist",
                      fct = "DTW",
                      stringsAsFactors = FALSE)
   
   if(!is.null(C) & is.list(C)){
      warning("plot() for rundtw() with C as list is not supported yet")
      C <- NULL
   }
         
   mycols <- c("C, no fit" = "#000000", "kNN" =  "#FF3333", "minimum" = "#3399FF", 
               "dist" = "#808080")
   
   
   if(!is.null(Q)) nQ <- ifelse(is.matrix(Q), nrow(Q), length(Q))
   if(!is.null(C) & !is.list(C)){
      nC <- ifelse(is.matrix(C), nrow(C), length(C))
   }
   
   if(minima & !is.null(Q)){
      ix_minima <- find_peaks(x$dist, w = nQ, get_min = TRUE)   
      dfp1[ix_minima, "categ"] <- "minimum"
   }else{
      minima <- FALSE
   }
   
   if(!is.null(x$knn_indices) & knn){
      dfp1[x$knn_indices, "categ"] <- "kNN"
   }else{
      knn <- FALSE
   }
   
   
   
   if(is.null(C)){
      
      gg <- ggplot(dfp1, ...) + 
         geom_line(aes_string(x = 'x_time', y = 'y_value', col = "'dist'"), na.rm=TRUE, ...) +
         geom_point(aes_string(x = 'x_time', y = 'y_value', col = 'categ'), na.rm=TRUE, ...) +
         guides(col = guide_legend(title = NULL))+
         scale_color_manual(values = mycols) + 
         ylab("Value") + xlab("Time")
      
      
   }else{
      
      C <- as.matrix(C)
      if(is.null(selDim)) selDim <- 1:ncol(C)
      C <- C[, selDim]
      # C_norm <- C
      C <- as.data.frame(C)
      C$x_time <- 1:nrow(C)
      C$categ <- "C, no fit"
      
      # colour the peaks
      if(minima & !is.null(Q)){
         for(mm in ix_minima){
            C[C$x_time >= mm & C$x_time <= (mm + nQ -2), "categ"] <- "minimum"       
         }
      }
      
      # colour the kNN
      if(knn & !is.null(Q)){
         for(ix_knn in x$knn_indices){
            C[C$x_time >= ix_knn & C$x_time <= (ix_knn + nQ -2), "categ"] <- "kNN"
         }
      }
      C <- as.data.table(C)
      dfp2 <- data.table::melt(C, id.vars = c("x_time", "categ"), value.name = "y_value")
      data.table::setnames(dfp2, "variable", "group")
      dfp1 <- data.table::as.data.table(dfp1)
      dfp1[, group := "DTW"]
      dfp2[, fct := "Time series C"]
      
      dfp <- rbind(dfp1[, list(x_time, y_value, group, categ, fct)],
                   dfp2[, list(x_time, y_value, group, categ, fct)])
      
      if(scale != "none"){
         dfp3 <- dfp2
         dfp3[,fct:= "Scaled fits"]
         dfp3[categ == "C, no fit", y_value:= NA]
         if(minima & !is.null(Q)){
            for(mm in ix_minima){
               dfp3[x_time >= mm & x_time <= (mm + nQ -1) & categ == "minimum", 
                    y_value:= IncDTW::scale(.SD[, y_value], scale), by = "group"]       
            }
         }
         
         # colour the kNN
         if(knn & !is.null(Q)){
            for(ix_knn in x$knn_indices){
               dfp3[x_time >= ix_knn & x_time <= (ix_knn + nQ -1) & categ == "kNN", 
                    y_value:= IncDTW::scale(.SD[, y_value], scale), by = "group"]       
            }
         }
         
         dfp <- rbind(dfp[, list(x_time, y_value, group, categ, fct)],
                      dfp3[, list(x_time, y_value, group, categ, fct)])
         
      }
      
      
      dfp[, y_value := as.numeric(y_value)]
      # dfp_DTW <- dfp[fct == "DTW"]
      dfp <- as.data.frame(dfp)
      # dfp_DTW <- as.data.frame(dfp_DTW)
      
      gg <- ggplot(dfp, ...) + 
         geom_line(aes_string(x = 'x_time', y = 'y_value', 
                              group = 'group', col = 'categ'), na.rm=TRUE, ...) +
         geom_point(data = dfp,#[fct == "DTW"], 
                    aes_string(x = 'x_time', y = 'y_value', col = 'categ'), na.rm=TRUE, ...) +
         facet_grid(fct~., scales = "free_y")+ 
         guides(col = guide_legend(title = NULL))+
         scale_color_manual(values = mycols) + 
         ylab("Value") + xlab("Time")
      
      
      
   }
   
   return(gg)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plot.idtw <- function(x, type = c("QC", "warp"), partial = NULL, selDim = 1, ...) {
   
   type <- match.arg(type)
   pm <- pmatch(type, c("QC", "warp"))
   
   switch(pm,
          plotQC(x, partial = partial, selDim = selDim, ...),
          plotWarp(x, partial = partial, selDim = selDim, ...)
   )
}

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plot.planedtw <- function(x, type = c("QC", "warp"), partial = NULL, selDim = 1, ...) {
   
   y <- dtw(x$Q, x$C, dist_method = x$control$dist_method,
            ws = x$control$ws, step_pattern = x$control$step_pattern, 
            return_QC = TRUE, return_wp = TRUE)
   
   plot(y, type = type, parital = partial, selDim = selDim, ...)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plot.dba <- function(x, type = c("barycenter", "m2m", "m2lot"), ...) {
   
   type <- match.arg(type)
   pm <- pmatch(type, c("barycenter", "m2m", "m2lot"))
   
   switch(pm,
          plotBary(x, ...),
          plotM2m(x, ...),
          plotM2lot(x, ...)
   )
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


## all aliases
plot_rundtw <- plot.rundtw;
plot_idtw <- plot.idtw;
plot_planedtw <- plot.planedtw;
plot_dba <- plot.dba;


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

plotM2m <- function(x, ...){
   dfp <- data.frame(Iterations = 1:x$input$iterMax, y = x$iterDist_m2m)
   gg  <- ggplot(dfp, ...) + 
      geom_line(aes_string(x = 'Iterations', y = 'y')) + 
      geom_point(aes_string(x = 'Iterations', y = 'y')) +
      ylab(bquote('Distances of' ~m[n]~to~m[n-1]))
           
   return(gg)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plotM2lot <- function(x, ...){
   if(x$input$step_pattern == "symmetric1"){
      dfp <- data.frame(iter = 1:x$input$iterMax,
                        meanDist = sapply(x$iterDist_m2lot, mean),
                        sdDist   = sapply(x$iterDist_m2lot, sd))
   }else{
      dfp <- data.frame(iter = 1:x$input$iterMax,
                        meanDist = sapply(x$iterDist_m2lot_norm, mean),
                        sdDist   = sapply(x$iterDist_m2lot_norm, sd))
   }
   
   dfp$ymin <- dfp$meanDist - dfp$sdDist
   dfp$ymax <- dfp$meanDist + dfp$sdDist
   
   gg  <- ggplot(dfp, ...) + geom_line(aes_string(x = 'iter', y = 'meanDist')) + 
      geom_point(aes_string(x = 'iter', y = 'meanDist')) +
      geom_ribbon(aes_string(x = 'iter', ymin = 'ymin', ymax = 'ymax'), alpha = 0.2) +
      xlab("Iterations")
   
   if(x$input$step_pattern == "symmetric1"){
      gg <- gg + ylab(bquote('Avg. dist of '~m[n]~'to lot'))
      ylab(bquote('Distances of' ~m[n]~to~m[n-1]))
   }else{
      gg <- gg + ylab(bquote('Avg. normalized dist of '~m[n]~'to lot'))
   }
   return(gg)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plotBary <- function(x, ...){
   dfp <- x$iterations
   dfp <- lapply(seq_along(dfp), function(i){
      cbind(dfp[[i]], Iter = i, j = 1:nrow(dfp[[i]]))
   })
   dfp <- melt(as.data.table(do.call(rbind, dfp)), id.vars = c("Iter", "j"))
   
   gg <- ggplot(dfp, ...) + 
      geom_line(aes_string(x = 'j', y = 'value', group = 'variable', col = 'Iter')) +
      xlab("Time") + ylab("Value")+
      facet_grid( ~variable)
   return(gg)
   
   # to adjust existing ggplot object
   # q <- ggplot_build(gg)
   # q$data[[1]]$size <- 1
   # gg4 <- ggplot_gtable(q)
   # plot(gg4)
}





#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plotQC <- function(x, Q = NULL, C = NULL, partial = NULL, selDim = 1, ...){
   # Q <- sin(1:10)#+rnorm(20)
   # # C <- sin(1:30)+rnorm(30)
   # C <- sin(-2:10)#+rnorm(15)
   # # C <- c(Q, cumsum(rnorm(20))) + rnorm(120, 0, 0.5)
   # x <- dtw(Q = Q, C = C,  return_diffM = FALSE, return_QC = T)
   
   if(is.null(Q) & !is.null(x$Q)) Q <- x$Q
   if(is.null(C) & !is.null(x$C)) C <- x$C
   if(is.null(Q) | is.null(C)){
      stop("Q and C either need to be returned from dtw() or idtw(), or defined manually")
   }
   if(!is.null(partial)){
      # adjust to partial matching
      new_wp <- BACKTRACK_cpp(x$dm[  partial$rangeQ[1]:partial$rangeQ[2],
                                     partial$rangeC[1]:partial$rangeC[2]  ])
      x$ii <- rev(new_wp$ii)
      x$jj <- rev(new_wp$jj)
   } 
   
   if(is.matrix(Q)) Q <- Q[, selDim]
   if(is.matrix(C)) C <- C[, selDim]
   
   # this makes the two time series to be plotted one over the other
   # Q <- Q-max(Q)
   # C <- C-min(C)
   
   tmp1 <- data.frame(val = c(Q, C),
                      x = c(1:length(Q), 1:length(C)),
                      id = c(rep("Q", length(Q)),
                             rep("C", length(C))
                      )
   )
   
   tmp2 <- data.frame(x = c(x$ii, x$jj), 
                      y = c(Q[x$ii], C[x$jj]),
                      id = rep(1:length(x$ii), 2))
   
   gg1 <- ggplot(tmp1, ...) +
      geom_line(aes_string(x = 'x', y = 'val', group = 'id', col = 'id')) + 
      xlab("Index") + ylab("Value") +
      geom_point(aes_string(x = 'x', y = 'val', col = 'id')) +
      geom_line(data = tmp2, aes_string(x = 'x', y = 'y', group = 'id'), lty = 2) +
      guides(col = guide_legend(title = NULL))
   ret <- gg1
   
   return(ret)
}

if(FALSE){
   Q <- sin(1:20)
   C <- cos(1:50)+1
   tmp <- IncDTW::dtw(Q = Q, C = C, return_QC = TRUE)
   par <- IncDTW::dtw_partial(tmp, partial_Q = FALSE, partial_C = TRUE)
   plot(tmp, partial = par, type="QC")
   plot(tmp, partial = par, type="warp")
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


plotWarp <- function(x, Q = NULL, C = NULL, partial = NULL, selDim = 1, ...){
   if(is.null(Q) & !is.null(x$Q)) Q <- x$Q
   if(is.null(C) & !is.null(x$C)) C <- x$C
   if(is.null(Q) | is.null(C)){
      stop("Q and C either need to be returned from dtw() or idtw(), or defined manually")
   }
   if(!is.null(partial)){
      # adjust to partial matching
      new_wp <- BACKTRACK_cpp(x$dm[  partial$rangeQ[1]:partial$rangeQ[2],
                                     partial$rangeC[1]:partial$rangeC[2]  ])
      x$ii <- rev(new_wp$ii)
      x$jj <- rev(new_wp$jj)
   } 
   
   if(is.matrix(Q)) Q <- Q[, selDim]
   if(is.matrix(C)) C <- C[, selDim]
   
   
   tmp3 <- data.frame(ii=x$ii, jj = x$jj)
   tmp1 <- data.frame(val = c(Q, C),
                      x = c(1:length(Q), 1:length(C)),
                      id = c(rep("Q", length(Q)),
                             rep("C", length(C))
                      )
   )
   text_yQ <- mean(Q)
   text_yC <- mean(C)
   text_xwp <- length(Q)
   text_ywp <- length(C)
   
   gg_wp <- ggplot(tmp3, ...) +
      geom_line(aes_string(x='ii', y='jj')) +
      geom_point(aes_string(x = 'ii', y = 'jj')) + 
      scale_x_continuous(position = "top", name = "", 
                         breaks = scales::pretty_breaks(n=5),
                         minor_breaks = NULL) +
      scale_y_continuous(trans = "reverse", name = "", 
                         breaks = scales::pretty_breaks(n=5),
                         minor_breaks = NULL) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
      geom_text(label="Warping path: ii", aes(y=1, x = text_xwp), hjust = 1, vjust = 1) +
      geom_text(label="Warping path: jj", aes(x=1, y = text_ywp), hjust = 0, vjust = 1, angle = 90)
   
   gg_Q <- ggplot(tmp1[tmp1$id == "Q", ]) + 
      geom_line(aes_string(x = 'x', y = 'val')) + ylab("")+
      geom_point(aes_string(x = 'x', y = 'val')) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      theme(panel.grid.major.y =  element_line(colour = "#eaeaea"))+
      theme(panel.grid.minor.y =  element_line(colour = "#eaeaea"))+
      scale_x_continuous( breaks = scales::pretty_breaks(n=5),
                          minor_breaks = NULL) +
      geom_text(label="Q", aes(x=1, y=text_yQ))
   
   gg_C <- ggplot(tmp1[tmp1$id == "C", ]) + 
      geom_line(aes_string(x = 'x', y = 'val')) + 
      geom_point(aes_string(x = 'x', y = 'val')) 
   gg_C <- gg_C + coord_flip() + 
      scale_x_continuous(position = "top", trans = "reverse", name="", 
                         breaks = scales::pretty_breaks(n=5), 
                         minor_breaks = NULL) +
      theme(panel.grid.major.x =  element_line(colour = "#eaeaea"))+
      theme(panel.grid.minor.x =  element_line(colour = "#eaeaea"))+
      geom_text(label="C", aes(x=1, y=text_yC))
   gg_C <- gg_C + theme(axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank())
   
   gg0 <- ggplot()+ geom_blank()+theme_minimal()
   gg2 <- gridExtra::grid.arrange(gg_wp, gg_C, gg_Q, gg0,
                                  layout_matrix = rbind(c(1,1,1,2), 
                                                        c(1,1,1,2),
                                                        c(1,1,1,2),
                                                        c(3,3,3,4)))
   ret <- gg2
   return(ret)
}
