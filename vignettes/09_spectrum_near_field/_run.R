library(here)
setwd(here('vignettes/09_spectrum_near_field'))

## ----load----
library(terms)
library(ggforce)
theme_set(theme_grey())

## ----run----
system("../../build/terms input > log")

## ----read----

library(rhdf5)
lf <- h5ls('map.h5')
lf

ld <- h5read('map.h5', "Near-Field")

d1 <- data.frame(ld$map_E)
names(d1) <- c('wavelength', 'x','y','z','Eavg','Ex','Ey')

## ----intensity----

ge <- get_geometry('input')

p1 <- ggplot(d1 %>% filter(wavelength %in% seq(400, 600, by=20)), aes(x, y, fill=Ey)) +
  facet_wrap(~wavelength,nrow=1) +
  geom_raster() +
  coord_equal() +
  # geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0,0),breaks=seq(-10,10,by=10)) +
  scale_y_continuous(expand = c(0,0)) +
  # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
  scale_fill_viridis_c(option = 'C') +
  labs(x='x /nm', y='y /nm',
       fill=expression(E^2),linetype = 'region') +
  theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
        panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
        panel.spacing.x = unit(0,'mm'))

p1


## ----ldoc----

# lambda, x1, y1, z1, avg, inc1, inc2 ...
d2 <- data.frame(ld$normalised_ldoc)
names(d2) <- c('wavelength', 'x','y','z','C')

p2 <- ggplot(d2 %>% filter(wavelength %in% seq(400, 600, by=20)), aes(x, y, fill=C)) +
  facet_wrap(~wavelength,nrow=1) +
  geom_raster() +
  coord_equal() +
  geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0,0),breaks=seq(-10,10,by=10)) +
  scale_y_continuous(expand = c(0,0)) +
  # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
  scale_fill_distiller(palette = 'PiYG') +
  labs(x='x /nm', y='y /nm',
       fill='LDOC',linetype = 'region') +
  theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
        panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
        panel.spacing.x = unit(0,'mm'))

p2

## ----animated----

library(gifski)
library(gganimate)
library(patchwork)
library(ggforce)

tmp <- tempdir()
steps <- unique(d1$wavelength)
# all(unique(d1$wavelength) == unique(d2$wavelength))

for (i in seq_along(steps)) {
  p_step1 <- ggplot(d1 %>% filter(wavelength == steps[i]), aes(x, y, fill=Ey)) +
    # facet_wrap(~wavelength,nrow=1) +
    geom_raster() +
    coord_equal() +
    geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill=NA, lty=2, col='grey98', inherit.aes = FALSE) +
    scale_x_continuous(expand = c(0,0),breaks=seq(-10,10,by=10)) +
    scale_y_continuous(expand = c(0,0)) +
    # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
    scale_fill_viridis_c(option = 'C', lim = c(0, max(d1$Ey))) +
    labs(x='x /nm', y='y /nm',
         fill=expression(E^2),linetype = 'region') +
    theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
          panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
          panel.spacing.x = unit(0,'mm'))
  
  
  p_step2 <- ggplot(d2 %>% filter(wavelength == steps[i]), aes(x, y, fill=C)) +
    # facet_wrap(~wavelength,nrow=1) +
    geom_raster() +
    coord_equal() +
    geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, inherit.aes = FALSE) +
    scale_x_continuous(expand = c(0,0),breaks=seq(-10,10,by=10)) +
    scale_y_continuous(expand = c(0,0)) +
    # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
    scale_fill_distiller(palette = 'PiYG', , lim = range(d2$C)) +
    labs(x='x /nm', y='y /nm',
         fill='LDOC',linetype = 'region') +
    theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
          panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3),
          panel.spacing.x = unit(0,'mm'))
  
  ggsave(glue::glue('{tmp}/{i}.png'), 
         p_step1 + p_step2 + plot_annotation(title=glue::glue('lambda={steps[i]}nm')), 
         width=8,height=6)
}

png_files <- glue::glue('{tmp}/{seq_along(steps)}.png')
gif_file <- 'anim.gif'
gifski(png_files, gif_file, delay = 0.2, progress = FALSE, width = 800, height=600)
unlink(png_files)
