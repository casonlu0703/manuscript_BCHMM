#' This function generates day profile plots
#' originally from Huang et. al., 2018
#' @param one_day_prob: selected one day circadian state probabilities
#' @param id: subject id used in caption
#' @param index: x coordinates of points in the plot
#' @param index_position/index_lable (X-axis parameters): position of tic marks and labels(clock time used here) 
#' 
figure_day_profile  <-  function(
    id, one_day_prob, index, index_position, index_lable, title_text
    ){
  ########24h oscillated state probability plot ################ 
  SP <- function(one_day_prob, n_states=3, L=288){
    L <- length(index) 
    plot_p <- matrix(NA,nrow=L-1,ncol=2*n_states)
    a <- 0
    for(i in 2:L){
      for(j in 1:n_states){ 
        plot_p[i-1,(j*2-1)] <- one_day_prob[i-1,j] 
        plot_p[i-1,j*2] <- one_day_prob[i,j]
        if (j==1){col_states <- "blue"#"#240691"
        } 
        if (j==2){col_states <- "pink"#"#31B57B"
        } 
        if (j==3){col_states <- "red"#"#FDE725"
        }
        if (j==1){
          point_1 <- a 
          point_2 <- point_1+plot_p[i-1,(j*2-1)] 
          point_4 <- a 
          point_3 <- point_4+plot_p[i-1,(j*2)] }
        if (j==2){ 
          point_1 <- a+plot_p[i-1,(j-1)*2-1] 
          point_2 <- point_1+plot_p[i-1,(j*2-1)] 
          point_4 <- a+plot_p[i-1,(j-1)*2] 
          point_3 <- point_4+plot_p[i-1,(j*2)] }
        if (j==3){ 
          point_1 <- a+plot_p[i-1,(j-2)*2-1]+plot_p[i-1,(j-1)*2-1] 
          point_2 <- point_1+plot_p[i-1,(j*2-1)] 
          point_4 <- a+plot_p[i-1,(j-2)*2]+plot_p[i-1,(j-1)*2] 
          point_3 <- point_4+plot_p[i-1,(j*2)]}
        polygon(
          c(index[i-1],index[i- 1],index[i],index[i]),c(point_1,point_2,point_3,point_4)
          , col = col_states
          , border = NA
          )
        lines(c(index[i-1],index[i]),c(point_2,point_3), col=col_states) }
    }
    axis(side=2, at=seq(0,1,0.2), labels=T, cex.axis = 1.5, las=1) 
    axis(side=1, at=index_position, labels=index_lable, cex.axis=1.5, las=1)
  }
  plot(index, rep(-2,length(index)), ylim=c(0,1)
       , xlab='', xaxt='n'
       , ylab='Probability', cex.lab = 1.5 , yaxt='n'
       , cex=1.5, las=1)
  SP(one_day_prob)
  title(title_text, cex.main = 2.0, font = 2
        , adj = 0 # left side
        ) 
  mtext('Clock Time (One Day)', side = 1,cex = 1.5, line = 2)
  # mtext('Probability', side=2, line = 1, cex=1.5, las = 3, outer = T) 
}
