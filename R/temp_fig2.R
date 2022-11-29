data <- expand.grid(X=nr, Y=nc)
data$Z <- c(height_map)
pal <- wes_palette("Zissou1", 10, type = "continuous")

# plot height map
data$Z1 <- c(height_map)
plt1 = ggplot() +
  geom_raster(data = data , aes(x = X, y=Y,fill = Z1)) + 
  scale_fill_stepsn(n.breaks = 12, colours = terrain.colors(12),na.value="lightgrey")+
  ggtitle("Altitude map of the simulated island")+
  coord_fixed()+
  theme_minimal()

# plot 1st temperature map
data$Z2 <- c(climate_maps[[1]])
plt2 = ggplot() +
  geom_raster(data = data , aes(x = X, y=Y,fill = Z2)) + 
  scale_fill_gradientn(limits= c(0,16),colours = pal)+
  ggtitle("Temperature map - first time step")+
  theme_minimal()
### scale_fill_stepsn(n.breaks = 12, colours = viridis(12),na.value="lightgrey")

# plot last temperature map
data$Z3 <- c(climate_maps[[N_cycles]])

plt3 = ggplot() +
  geom_raster(data = data , aes(x = X, y=Y,fill = Z3)) + 
  scale_fill_gradientn(limits= c(0,16),colours = pal)+
  ggtitle("Temperature map - last time step")+
  theme_minimal()
###scale_fill_stepsn(n.breaks = 12, colours = viridis(12),na.value="lightgrey")+

#plt4 = grid.arrange(plt1,plt2,plt3)
plt4 = plot_grid(plt1, plt2, plt3, labels=c("(A)", "(B)", "(C)"), ncol = 3, nrow = 1)


plot(plt4)


ggsave(name, width = 15, height = 5)


#Pep starts
data
head(data)
dat.temp = data[,c("X","Y","Z2","Z3")]
names (dat.temp) <- c('X','Y','Initital timestep','Final timestep')

dat.temp$`Initital timestep`
dat.temp$`Final timestep`

library (tidyverse)
kk = pivot_longer(dat.temp,cols=c(3,4),names_to = 'temperature',values_to = "C") 

t1temp = kk %>% filter (temperature == 'Initital timestep') 
summary(t1temp)
t1temp$value

kk$temperature = factor (kk$temperature,levels = c('Initital timestep','Final timestep'))

p = ggplot() +
  geom_raster(data = kk , aes(x = X, y=Y,fill = C)) + 
  facet_wrap(temperature~.)+
  scale_fill_gradientn(limits= c(0,16),colours = pal)+
  ggtitle("Temperatures change")+
  theme(panel.background = element_blank(),strip.background = element_blank() ,line = element_blank(),axis.title = element_blank(),axis.text = element_blank(), axis.ticks = element_blank())+
  coord_fixed()

plt4 = plot_grid(plt1, p, labels=c("(A)", "(B)"), ncol = 2, nrow = 1)

plot(plt4)
?pivot_longer

library (arrangeGrob)
plot5<-grid.arrange(plot4, arrangeGrob(plt4, p, ncol=1), 
                    ncol=2, widths=c(1,1.2))
