library(here)
setwd(here('vignettes/08_polarimetry'))

## ----load----
library(terms)
theme_set(egg::theme_article())

## ----run----
system("../../build/terms input > log")

## ----read----

library(rhdf5)
lf <- h5ls('farfield.h5')
lf

ld <- h5read('farfield.h5', "Polarimetry")

## ----dsca----

d <- data.frame(ld$diff_Sca_CS)
names(d) <- c('alpha', 'beta', 'gamma',  'wavelength', 'diff_Sca_CS')
head(d)

p <- ggplot(d, aes(wavelength, beta*180/pi, fill = log10(diff_Sca_CS))) +
  geom_raster() +
  scale_fill_distiller(palette = 'PiYG') +
  scale_x_continuous('wavelength /nm', expand=c(0,0)) +
  scale_y_continuous(expression(theta), expand=c(0,0), lim=c(0,180), breaks=seq(0,180,by=45)) +
  labs(fill=expression(partialdiff*sigma[sca]))

p

## ----phasematrix----

d <- data.frame(ld$Stokes_phase_Mat)
# dput(do.call(paste0, cbind('z',expand.grid(1:4,1:4))))
names(d) <- c('alpha', 'beta', 'gamma',  'wavelength', 
              "z11", "z21", "z31", "z41", "z12", "z22", "z32", "z42", "z13", 
              "z23", "z33", "z43", "z14", "z24", "z34", "z44")
# glimpse(d)
head(d)

p <- ggplot(d, aes(wavelength, beta*180/pi, fill = (z41+z42)/(z11+z12))) +
  geom_raster() +
  scale_fill_distiller(palette = 'RdBu') +
  scale_x_continuous('wavelength /nm', expand=c(0,0)) +
  scale_y_continuous('scattering angle /deg.', expand=c(0,0), lim=c(0,180), breaks=seq(0,180,by=45)) +
  labs(fill='(z41+z42)/\n(z11+z12)')

p

## ----stokes----
# Sca_angle(alpha, beta, gamma)  lambda(nm)     I             Q              U              V
d <- setNames(data.frame(ld$Stokes_Sca_Vec), c('alpha', 'beta', 'gamma', 'wavelength','I', 'Q', 'U','V'))
head(d)
d <- d %>% mutate(docp = V/I)

p <- ggplot(d, aes(wavelength,beta*180/pi, fill=docp))+
  geom_raster()+
  # coord_equal()+
  scale_x_continuous(expand=c(0,0),lim=c(200,1000), breaks=seq(200,1000,by=200))+
  scale_y_continuous(expand=c(0,0), breaks = seq(0,180,by=45))+
  labs(x = 'wavelength /nm', y = 'scattering angle /deg.', fill='Deg. of CP') +
  scale_fill_distiller(palette='RdBu')
p



