library(here)
setwd(here('vignettes/15_quadrimer'))

## ----load---
library(terms)
library(purrr)
library(ggplot2)
library(egg)
theme_set(theme_grey())

## ---run---

system("../../build/terms input_dimer_ff > log_dimer_ff")
system("../../build/terms input_quadrimer_ff > log_quadrimer_ff")
system("../../build/terms input_twodimers_ff > log_twodimers_ff")

# get the get collective Tmat at 400nm  
system("../../build/terms input_dimer_ff_400 > log_dimer_ff_400")

# run the near-field at 400nm 
system("../../build/terms input_dimer_nf > log_dimer_nf")
system("../../build/terms input_quadrimer_nf > log_quadrimer_nf")
system("../../build/terms input_twodimers_nf > log_twodimers_nf")


## ----read----

read_map <- function(f){
  setNames(data.frame(h5read(f, 'Near-Field')$map_E), c('wavelength','x','y','z','scatID', 'volID','I'))
}
d1 <- read_map('map_dimer.h5')
d2 <- read_map('map_quadrimer.h5')
d3 <- read_map('map_twodimers.h5')

glimpse(d1)
## ----nf----

d <- rbind(cbind(d1, case = 'dimer'),
           cbind(d2, case = 'quadrimer'),
           cbind(d3, case = 'two dimers'))

rg1 <- "Ag 0.0 -30.0 0.0 25.0
Ag 0.0  30.0 0.0 25.0"

rg2 <- "Ag -30.0 -30.0 0.0 25.0
Ag -30.0  30.0 0.0 25.0
Ag 30.0 -30.0 0.0 25.0
Ag 30.0  30.0 0.0 25.0"

add_geometry <- function(descriptor){
  geometry <- setNames(read.table(textConnection(descriptor), 
                                  colClasses = c("character","double","double","double","double")), 
                       c("tag", 'x','y','z','r'))
  geometry$label <- paste0(geometry$tag, seq(1,nrow(geometry)))
  
  geometry
}

ge1 <- add_geometry(rg1)
ge2 <- add_geometry(rg2)

geometry <- rbind(cbind(ge1, case = 'dimer'),
                  cbind(ge2, case = 'quadrimer'),
                  cbind(ge2, case = 'two dimers'))

p1 <- ggplot(d, aes(y,z)) +
  facet_wrap(~case, nrow=1) +
  geom_raster(aes(fill=log10(I))) +
  coord_equal() +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='black', alpha=0.5, lwd=0.2,lty=2) +
  scale_fill_viridis_c(option = 'plasma', 
                       limits=c(-4,max(log10(d2$I))), 
                       na.value = 'grey90') +
  theme_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="y /nm", y="z /nm", fill="log(I)", 
       title="Near-field intensity map at 400 nm",
       subtitle = "incidence along x, circular polarisation") +
  theme(legend.position = 'bottom', legend.direction = 'horizontal')

p1

## ----ff----

d1 <- consolidate_xsec('cross_sections_dimer.h5')
d2 <- consolidate_xsec('cross_sections_quadrimer.h5')
d3 <- consolidate_xsec('cross_sections_twodimers.h5')


all <- rbind(cbind(d1$mCOA, simulation = 'dimer'),
             cbind(d2$mCOA, simulation = 'quadrimer'),
             cbind(d3$mCOA, simulation = 'twodimers'))
glimpse(all)

p2 <- ggplot(all %>% filter(variable=='total', simulation != 'dimer'), 
       aes(wavelength, average, linetype=simulation, colour=crosstype)) +
  geom_line()

p2
