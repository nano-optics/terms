library(here)
setwd(here('vignettes/107_WGM'))


## ----load----
library(terms)
theme_set(theme_grey())

# ----tpl----

tpl <- "ModeAndScheme 2 3
MultipoleCutoff 15
Wavelength 700 900 50
Verbosity 1
Medium 1.0

Incidence 0 1.570796 0 2
OutputFormat HDF5 xsec

Scatterers 2
Si 0.0  -1050.0 0.0 1000.0
Si 0.0  1050 0.0 1000.0
"

cat(tpl)

## ----run----
cat(glue(tpl), file = 'input')
system("../../build/terms input > log")


## ----tplnf----

tpl <- "ModeAndScheme 1 0
MultipoleCutoff 20
Wavelength 600 900 2
Verbosity 1
Medium 1.0

Incidence 0 1.570796 0 2
OutputFormat HDF5 map
SpacePoints -1500 1500 200 -2500 2500 400 0 0 0 0 0 0 
MapQuantity 2 E 

Scatterers 2
Si 0.0  -1050.0 0.0 1000.0
Si 0.0  1050 0.0 1000.0
"

cat(tpl)

## ----runnf----
cat(glue(tpl), file = 'input_map')
system("../../build/terms input_map > log")

## ----read----

xs <- consolidate_xsec('xsec.h5')

## ----oa----
mCOA <- xs$mCOA

mCOAt <- subset(mCOA, variable == 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average, colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = '', linetype='') +
  ggtitle("Orientation-averaged cross-sections")

p1


## ----intensity----

ld <- h5read('map.h5', "Near-Field")

d1 <- data.frame(ld$map_E)
names(d1) <- c('wavelength', 'x','y','z','Eavg','Ex','Ey')

ge <- get_geometry('input_map')

p2 <- ggplot(d1, aes(x/1e3, y/1e3, fill=log10(Ey))) +
  facet_wrap(~wavelength,nrow=1) +
  geom_raster() +
  coord_equal() +
  geom_circle(data=ge, aes(x0=x/1e3,y0=y/1e3,r=r/1e3), fill=NA, lty=1, lwd=0.5, col='white', inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0,0),lim=c(-1.45,1.45)) +
  scale_y_continuous(expand = c(0,0)) +
  # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
  scale_fill_viridis_c(option = 'C') +
  # scale_fill_distiller(palette = 'PiYG') +
  labs(x='x /µm', y='y /µm',
       fill=expression(log[10](E^2)),linetype = 'region') +
  theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
        panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
        panel.spacing.x = unit(2,'mm'))

p2

# p2 <- ggplot(d1 %>% filter((x^2 + (y-1000)^2) > 1000^2, (x^2 + (y+1000)^2) > 1000^2), aes(x, y, fill=Ey)) +
#   facet_wrap(~wavelength,nrow=1) +
#   geom_raster() +
#   coord_equal() +
#   geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill=NA, lty=2, inherit.aes = FALSE) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
#   scale_fill_viridis_c(option = 'C') +
#   # scale_fill_distiller(palette = 'PiYG') +
#   labs(x='x /nm', y='y /nm',
#        fill=expression(E^2),linetype = 'region') +
#   theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
#         panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
#         panel.spacing.x = unit(0,'mm'))
# p2

