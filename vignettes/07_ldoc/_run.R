library(here)
setwd(here('vignettes/07_ldoc'))

## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())


## ----tpl----

tpl <- "ModeAndScheme 1 0
MultipoleCutoff 8
Wavelength 540
Medium 1.7689 # water

OutputFormat HDF5 ldoc
Incidence file incidence 1  # linear polarisation, 45 & 135
SpacePoints -15 15 90 -25 25 150 0 0 0
MapQuantity C

Scatterers 2
Au 0 12.5 0 10
Au 0 -12.5 0 10
"

tpl_oa <- "ModeAndScheme 1 3
MultipoleCutoff 8
Wavelength 540
Medium 1.7689 # water

Incidence 0 0 1 1 # LP -> both circular polarisation in OA
OutputFormat HDF5 ldoc_oa
SpacePoints -15 15 90 -25 25 150 0 0 0
MapOaQuantity E C

Scatterers 2
Au 0 12.5 0 10
Au 0 -12.5 0 10
"



## ----run----
cat(tpl, file='input')
system('../../build/terms input > log')

cat(tpl_oa, file='input_oa')
system('../../build/terms input_oa > log_oa')

## ----read----

library(rhdf5)
h5ls('ldoc.h5')
ld <- h5read('ldoc.h5', "Near-Field")

glimpse(ld$normalised_ldoc)
# lambda, x1, y1, z1, avg, inc1, inc2 ...
d <- data.frame(ld$normalised_ldoc[,c(1,2,3,4,6,7)])
names(d) <- c('wavelength', 'x','y','z','L45','L135')

m <- d %>% pivot_longer(c('L45','L135'))

## ----comparison----

ge <- get_geometry('input')

both <- ggplot(m, aes(x, y, fill=value)) +
  facet_grid(~name)+
  geom_raster() +
  coord_equal() +
  geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
  scale_fill_distiller(palette = 'PiYG') +
  labs(x='x /nm', y='y /nm',
       fill='LDOC',linetype = 'region') +
  theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
        panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
        panel.spacing.x = unit(5,'mm'))

both

## ----oa----

h5ls('ldoc_oa.h5')
ld <- h5read('ldoc_oa.h5', "Near-Field")
# glimpse(ld)

d <- data.frame(ld$mapOaQuantity)
names(d) <- c('wavelength','x','y','z','E','R','L')
head(d)
m <- d %>% pivot_longer(c('R','L'))

p <- ggplot(m, aes(x,y,fill=value)) +
  facet_grid(~name)+
  geom_raster() +
  geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, inherit.aes = FALSE) +
  scale_fill_distiller(palette = 'PiYG') +
  coord_equal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x='x /nm', y='y /nm',
     fill='OA-LDOC',linetype = 'region') +
  theme(panel.spacing.x = unit(5,'mm'))
  
p

# ----other----

# d <- data.frame(ld$map_B)
# glimpse(d)
# names(d) <- c('wavelength','x','y','z','a','b','c')
# head(d)
# m <- d %>% pivot_longer(c('a','b','c'))
# ggplot(m, aes(x,y,fill=value)) +
#   facet_grid(~name)+
#   geom_raster() +
#   coord_equal() 
# 
# p <- ggplot(m, aes(x, y, fill=value)) +
#   facet_grid(~name)+
#   geom_raster() +
#   coord_equal() +
#   geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, inherit.aes = FALSE) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
#   scale_fill_distiller(palette = 'PiYG') +
#   labs(x='x /nm', y='y /nm',
#        fill='LDOC',linetype = 'region') +
#   theme_grey() +
#   theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
#         panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3), panel.spacing.x = unit(5,'mm'))
# p
