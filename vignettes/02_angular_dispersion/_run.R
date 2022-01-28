setwd(here::here('vignettes/02_angular_dispersion/'))

## ----load----
library(terms)
library(ggplot2)
theme_set(theme_grey())

## ----run----
system("../../build/terms input > log")

## ----read----
library(rhdf5)
lf <- h5ls('results.h5')
lf


# grab the results and store in hdf5 format
xs <- consolidate_xsec('results.h5')

str(xs)

## ----fol----
glimpse(xs$mLFO)

mm <- pivot_longer(xs$mLFO %>% filter(variable != 'total'), c("polarisation1","polarisation2"))
str(mm)
mm$polarisation <- factor(mm$name, labels=c('X','Y'))
mm$angle <- seq(0,90,length=7)[as.integer(factor(mm$variable))]

p <- ggplot(subset(mm, angle != 'total'), aes(wavelength, value)) +
  facet_grid(polarisation~crosstype, scales='free_y')+
  geom_line(aes(colour=factor(angle))) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(phi/degree)) +
  ggtitle("Fixed-orientation linear polarisation cross-sections")

p


## ----cleanup----
# unlink("*.dat")

