library (tidyverse)

do_plot_fig3 = function(nr,nc,cost1,cost2,cost3,name){
  
  data <- expand.grid(X=nr, Y=nc)
  data$C1 <- c(as.matrix(cost1))
  data$C2 <- c(as.matrix(cost2))
  data$C3 <- c(as.matrix(cost3))
  dat.temp = data[,c(2,1,3,4,5)]
  names (dat.temp) <- c('X','Y',"(a) uniform cost","(b) west-to-east cost","(c) altitude cost")
  kk2 = pivot_longer(dat.temp,cols=c(3,4,5),names_to = 'cost_type',values_to = "Cost") 
  kk2$temperature = factor (kk2$cost_type,levels = c("(a) uniform cost","(b) west-to-east cost","(c) altitude cost"))
  plt_cost = ggplot() + coord_fixed()+
    geom_raster(data = kk2 , aes(x = X, y=Y,fill = Cost)) + theme(panel.background = element_blank(),strip.background = element_blank() ,line = element_blank(),axis.title = element_blank(),axis.text = element_blank(), axis.ticks = element_blank())+
    facet_wrap(temperature~.)+
    scale_fill_gradientn(colours = viridis(10),limits=c(100,500))
  ggsave(name, width = 15, height = 5)
  return(plt_cost)
}
