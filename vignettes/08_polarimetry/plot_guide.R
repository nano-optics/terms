library(here)
setwd(here('vignettes/08_polarimetry'))

## ----load----
library(terms)
theme_set(egg::theme_grey(base_size = 10, base_family = 'Source Sans Pro'))

library(rhdf5)
lf <- h5ls('farfield.h5')
lf

ld <- h5read('farfield.h5', "Polarimetry")

# Sca_angle(alpha, beta, gamma)  lambda(nm)     I             Q              U              V
d <- setNames(data.frame(ld$Stokes_Sca_Vec), c('alpha', 'beta', 'gamma', 'wavelength','I', 'Q', 'U','V'))
head(d)
d <- d %>% mutate(docp = V/I)


p <- ggplot(d, aes(wavelength,beta*180/pi, fill=docp))+
  geom_raster()+
  # coord_cartesian()+
  scale_x_continuous(expand=c(0,0),lim=c(200,1000), breaks=seq(200,1000,by=200))+
  scale_y_continuous(expand=c(0,0), breaks = seq(0,180,by=45))+
  labs(x = 'wavelength /nm', y = 'scattering angle /deg', fill='Deg. of CP') +
  # scale_fill_distiller(palette='RdBu') #BrBG+
# scale_fill_distiller(palette='BrBG') #
scale_fill_distiller(palette='PRGn', lim=c(-1,1)) + #
  theme(aspect.ratio = 1, legend.background = element_rect(fill = '#FFFFFFBB'), 
        legend.text = element_text(size=8), legend.title = element_text(size=9),
        legend.position = c(0.83,0.79), legend.margin = margin(1,1,2,1,'mm')) +
  guides(fill=guide_colourbar(barwidth = unit(3,'mm'), ticks.colour = 'black'))
# p

g <- egg::set_panel_size(p, width=unit(2.2,'in'), height=unit(2.2,'in'))
ggsave('plot_guide.pdf',g,width=2.9,height=2.7, device = cairo_pdf)
# system('open .')

# g <- egg::set_panel_size(p, width=unit(2.2,'in'), height=unit(2.2,'in'))
# ggsave('map_guide.pdf',g,width=2.8,height=2.6)

