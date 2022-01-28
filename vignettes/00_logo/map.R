

library(reshape2)
library(purrr)
library(tibble)
library(ggplot2)
library(egg)
library(ggforce)
theme_set(theme_grey())


read_map <- function(f='map.dat'){
  d <- read.table(f)
  names(d) <- c("x","y","z","I")
  invisible(d)
}

if(!file.exist('map.rds')){
  d <- read_map()
  saveRDS(d, file='map.rds')
}
d <- readRDS('map.rds')
input <- readLines('input')
rg <- paste(input[(grep("Scatterer",input)+1):(length(input)-1)], collapse="\n")

add_geometry <- function(descriptor){
  geometry <- setNames(read.table(textConnection(descriptor), 
                                  colClasses = c("character","double","double","double","double")), 
                       c("tag", 'x','y','z','r'))
  geometry$label <- paste0(geometry$tag, seq(1,nrow(geometry)))
  
  geometry
}

ge <- add_geometry(rg)

geometry <- ge

p <- ggplot(d, aes(x,y)) +
  # facet_wrap(~case, nrow=1) +
  geom_raster(aes(fill=log10(I))) +
  coord_equal() +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='black', alpha=0.5, lwd=0.2,lty=2) +
  scale_fill_viridis_c(option = 'plasma', 
                       na.value = 'grey90') +
  theme_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_void()+
  theme(legend.position = 'none', legend.direction = 'horizontal')


# p
ggsave("nf_map.png", p, width=200, height=120, units = 'mm')
