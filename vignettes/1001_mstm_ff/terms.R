setwd(here::here("vignettes/01_dimer_spheres"))
## ----load----
suppressPackageStartupMessages(require(terms))
library(dplyr)
theme_set(theme_grey())

read.table('400.inp', skip=5, nrows=3,sep=',')

## ----run----
system("../../build/terms input > log")

## ----read----
# grab the results and store in compressed form
xs <- store_xsec(out='xsec.rds')

list.files(pattern=".dat")

# plot OA data
## ----oa----
lfOA

glimpse(xs$mCOA)
mCOAt <- subset(xs$mCOA, variable == 'total')
mCOAn <- subset(xs$mCOA, variable != 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average, colour=crosstype)) +
  facet_wrap(~crosstype) +
  guides(colour='none') +
  scale_colour_brewer(palette='Set1') +
  # geom_line(aes(colour=variable), lty=2) +
  geom_line(data=mCOAt, lty=1) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Orientation-averaged cross-sections")

p1

# p2 <- ggplot(mCOAn, aes(wavelength, dichroism)) +
#   facet_wrap(~crosstype)+
#   geom_line(aes(colour=variable), lty=2) +
#   geom_line(data=mCOAt, lty=1) +
#   labs(x = expression("wavelength /nm"), 
#        y = expression("cross-section /"*nm^2),
#        colour = expression(N)) +
#   ggtitle("Orientation-averaged circular dichroism")

# egg::ggarrange(p1,p2,ncol=1)

## plot fixed orientation data

## ----fol----
lfLinear

glimpse(xs$mLFO)

mLFOt <- subset(xs$mLFO, variable == 'total')
mLFOn <- subset(xs$mLFO, variable != 'total')

p3 <- ggplot(mLFOn, aes(wavelength, average, colour=crosstype)) +
  facet_wrap(~crosstype)+
  # geom_line(aes(colour=variable)) +
  guides(colour='none') +
  scale_colour_brewer(palette='Set1') +
  geom_line(data=mLFOt, map=aes(y=polarisation1, lty='X')) +
  geom_line(data=mLFOt, map=aes(y=polarisation2, lty='Y')) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N), linetype='polarisation') +
  ggtitle("Fixed-orientation linear polarisation cross-sections")

p4 <- ggplot(mLFOn, aes(wavelength, dichroism, colour=crosstype)) +
  facet_wrap(~crosstype)+
  # geom_line(aes(colour=variable), lty=2) +
  geom_line(data=mLFOt, lty=1) +
  scale_colour_brewer(palette='Set1') +
  guides(colour='none') +
  labs(x = expression("wavelength /nm"),
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Fixed-orientation linear dichroism")

egg::ggarrange(p3,p4,ncol=1)

## ----foc----
lfCircular

glimpse(xs$mCFO)

mCFOt <- subset(xs$mCFO, variable == 'total')
mCFOn <- subset(xs$mCFO, variable != 'total')

p5 <- ggplot(mCFOt, aes(wavelength, average, colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line() +
  guides(colour='none') +
  scale_colour_brewer(palette='Set1') +
  # geom_line(data=mCFOt, map=aes(y=polarisation1), lty=1) +
  # geom_line(data=mCFOt, map=aes(y=polarisation2), lty=2) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Fixed-orientation circular polarisation cross-sections")

# p6 <- ggplot(mCFOn, aes(wavelength, dichroism)) +
#   facet_wrap(~crosstype)+
#   geom_line(aes(colour=variable), lty=2) +
#   geom_line(data=mCFOt, lty=1) +
#   labs(x = expression("wavelength /nm"), 
#        y = expression("cross-section /"*nm^2),
#        colour = expression(N)) +
#   ggtitle("Fixed-orientation circular dichroism")

p5
# egg::ggarrange(p5,p6,ncol=1)



## ----cleanup----
# unlink("*.dat")

