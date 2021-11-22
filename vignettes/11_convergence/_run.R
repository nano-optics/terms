setwd(here::here('vignettes/11_convergence/'))
## ----load---
library(terms)
library(glue)
library(purrr)
library(ggplot2)
library(egg)
theme_set(egg::theme_article())


## ----tpl---

tpl <- "ModeAndScheme 2 3
MultipoleCutoff 6
Wavelength 400 700 150
Medium 1.7689 # epsilon of water
OutputFormat HDF5 cross_sections_{nmax}

MultipoleSelections 1
EE1:{nmax}_EM1:{nmax}_ME1:{nmax}_MM1:{nmax} blocks

Scatterers 4
Au_S1 100.0000  0.00000 -30 25.0
Au_S1  80.9017 58.77853 -10 25.0
Au_S1  30.9017 95.10565  10 25.0
Au_S1 -30.9017 95.10565  30 25.0
"

cat(tpl)

## ----run---

for(n in 1:5){
  nmax <- n
  cat(glue(tpl),'\n', file=glue("input_{nmax}"))
  system(glue("../../build/terms input_{nmax} > log_{nmax}"))
}

## ----read----

xs <- purrr::map_df(1:5, function(i) terms::consolidate_xsec(glue('cross_sections_{i}.h5'))$mCOA, .id = 'nmax')
glimpse(xs)
## ----oa----

mCOAt <- subset(xs, variable == 'total')
mCOAn <- subset(xs, variable != 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average, colour=factor(nmax))) +
  # facet_wrap(~ crosstype + abs(Scheme),scales = 'free',ncol=3)+
  facet_wrap(~crosstype)+
  # geom_path(data=mCOAn, aes(colour=variable)) +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[max])) +
  ggtitle("Orientation-averaged cross-sections")

p2 <- ggplot(mCOAt, aes(wavelength, dichroism, colour=factor(nmax))) +
  facet_wrap(~crosstype)+
  # facet_wrap(~ crosstype + abs(Scheme),scales = 'free',ncol=3)+
  # geom_path(data=mCOAn, aes(colour=variable)) +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[max])) +
  ggtitle("Orientation-averaged circular dichroism")

egg::ggarrange(p1,p2,ncol=1)

## ----naive---

mCOAt <- subset(xs, nmax==5 & variable == 'total')
mCOAn <- subset(xs, nmax==5 & variable != 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAn, aes(colour=variable)) +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[max])) +
  ggtitle("Orientation-averaged cross-sections")

p2 <- ggplot(mCOAt, aes(wavelength, dichroism)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAn, aes(colour=variable)) +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[max])) +
  ggtitle("Orientation-averaged circular dichroism")

egg::ggarrange(p1,p2,ncol=1)


## ----scatterer_centred---


tpl_centred <- "ModeAndScheme 2 1
MultipoleCutoff 6
Wavelength 400 700 150
Medium 1.7689 # epsilon of water
OutputFormat HDF5 cross_sections_{nmax}_centred

MultipoleSelections 1
EE1:{nmax}_EM1:{nmax}_ME1:{nmax}_MM1:{nmax} blocks
#ScattererCentredCrossSections

Scatterers 4
Au_S1 100.0000  0.00000 -30 25.0
Au_S1  80.9017 58.77853 -10 25.0
Au_S1  30.9017 95.10565  10 25.0
Au_S1 -30.9017 95.10565  30 25.0
"

nmax <- 5
cat(glue(tpl_centred), file=glue("input_{nmax}_centred"))
system(glue("../../build/terms input_{nmax}_centred > log_{nmax}_centred"))

xs2 <- terms::consolidate_xsec(glue('cross_sections_5_centred.h5'))$mCOA


mCOAt <- subset(xs2, nmax==5 & variable == 'total')
mCOAn <- subset(xs2, nmax==5 & variable != 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAn, aes(colour=variable)) +
  # geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[max])) +
  ggtitle("Orientation-averaged cross-sections")

p2 <- ggplot(mCOAt, aes(wavelength, dichroism)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAn, aes(colour=variable)) +
  # geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[max])) +
  ggtitle("Orientation-averaged circular dichroism")

egg::ggarrange(p1,p2,ncol=1)


