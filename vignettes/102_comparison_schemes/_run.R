library(here)
setwd(here('vignettes/102_comparison_schemes'))

## ----load---
library(terms)
library(glue)
library(purrr)
library(ggplot2)
library(egg)
theme_set(theme_grey())

## ----tpl----

tpl <- "ModeAndScheme 2 {scheme}
MultipoleCutoff 8
Wavelength 300 900 300
Medium 1.7689 # epsilon of water
Incidence  0.0 0.0 0.0 # default along z
OutputFormat HDF5 xsec_{scheme}
Scatterers 4
Au 100.0000  0.00000 -30 20.0
Au  80.9017 58.77853 -10 20.0
Au  30.9017 95.10565  10 20.0
Au -30.9017 95.10565  30 20.0
"

cat(tpl)

## ----run----
for(scheme in 0:3){
  cat(glue(tpl), "\n", file=glue("input_{scheme}"))
  system(glue("../../build/terms input_{scheme} > log_{scheme} &"))
}

## Read and store the OA data

## ----oa----

loa <- map_df(1:3, function(s) consolidate_xsec(glue("xsec_{s}.h5"))$mCOA, .id = 'scheme')
glimpse(loa)

p1 <- ggplot(loa %>% filter(variable=='total') %>% pivot_longer(c('average','dichroism')), 
             aes(wavelength, value, linetype=factor(scheme), colour=crosstype)) +
  # facet_wrap(~ crosstype + abs(Scheme),scales = 'free',ncol=3)+
  facet_grid(name~crosstype, scales='free')+
  # geom_path(data=mCOAn, aes(colour=variable)) +
  geom_line() + guides(colour='none') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       linetype = "Scheme", colour="") +
  ggtitle("Orientation-averaged cross-sections")

p1

## ----fo----
lfo <- map_df(0:3, function(s) consolidate_xsec(glue("xsec_{s}.h5"))$mLFO, .id = 'scheme')
glimpse(lfo)


p2 <- ggplot(lfo %>% filter(variable=='total') %>% pivot_longer(c('polarisation1','polarisation2')), 
             aes(wavelength, value, linetype=factor(scheme, labels=0:3), colour=crosstype)) +
  facet_grid(name~crosstype, scales='free')+
  geom_line() + guides(colour='none') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       linetype = "Scheme", colour="") +
  ggtitle("Fixed-orientation linear polarisation cross-sections")

p2
