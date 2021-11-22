setwd(here::here('vignettes/02_angular_dispersion/'))

## ----load----
library(terms)
theme_set(egg::theme_article())

## ----run----
system("../../build/terms input > log")

## ----read----
library(rhdf5)
lf <- h5ls('results.h5')
lf
# grab the results and store in hdf5 format
xs <- consolidate_xsec('results.h5')



## ----fol----
glimpse(xs$mLFO)

xs$mLFO$angle <- xs$mLFO$variable
xs$mLFO$variable <- NULL
mm <- melt(xs$mLFO, id=c("wavelength", "crosstype", "angle"), 
           measure.vars = c("polarisation1","polarisation2"))
mm$polarisation <- factor(mm$variable, labels=c('X','Y'))
mm$angle <- seq(0,90,length=7)[as.integer(mm$angle)]

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

