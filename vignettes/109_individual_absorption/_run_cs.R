setwd(here::here('vignettes/109_individual_absorption'))

## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())

## ----tpl----

tpl <- "ModeAndScheme 2 2
MultipoleCutoff 8
Wavelength 300 800 100
Incidence file incidence
Medium 1.7689 # for water
Verbosity 2
OutputFormat HDF5 cross_sections
ScattererCentredCrossSections

Scatterers 2
Au@Ag  0 0 0 {Rcore1+d1}  {Rcore1} 
Au@Ag  {Rcore1 + Rcore2 + d1 + d2 + gap} 0 0 {Rcore2+d2}  {Rcore2} 
"

Rcore1 <- 20
Rcore2 <- 30
d1 <- 1
d2 <- 2
gap <- 2

cat(glue(tpl))
cat(glue(tpl), '\n', file='input')

## ----run----

system("../../build/terms input > log")

## ----xsec----

library(rhdf5)
lf <- h5ls('cross_sections.h5')
lf

# h5readAttributes('cross_sections.h5', "Far-Field/oa_incidence/csAbsOA_split")
# h5closeAll()
# h5readAttributes(file = "map.h5", 
#                  name = "Near-Field/normalised_ldoc")
ld <- h5read('cross_sections.h5', "Far-Field")
glimpse(ld)

# total OA cross-sections
xsec <- consolidate_xsec('cross_sections.h5')
glimpse(xsec)

total <- xsec$mCOA %>% filter(variable=='total')
total$scatterer <- 'both'
total$region <- 'total_avg'
total$calculation <- 'total'

p <- ggplot(total, aes(wavelength, average, colour=crosstype)) +
  geom_line() +
  facet_wrap(~crosstype)

p

# Mackowski's partial shell absorptions, with numerical cubature
xsec_part <- consolidate_partials('cross_sections.h5')
xsec_part$crosstype <- 'Abs'
xsec_part$calculation <- 'parts'
xsec_part <- xsec_part %>% mutate(total = partial_1, core = partial_0, shell = total -core)

glimpse(xsec_part)

## Stout's per-particle OA absorptions
head(ld$oa_incidence$csAbsOA_split)

split <- setNames(data.frame(ld$Wavelengths, ld$oa_incidence$csAbsOA_split), c('wavelength','both','1','2'))
glimpse(split)

split$crosstype <- 'Abs'

ds <- split %>% pivot_longer(c('both','1','2'))
ds$scatterer <- ds$name
ds$region <- 'total_avg'
ds$calculation <- 'split'

mparts <- xsec_part %>% filter(polarisation=='4L') %>% 
  pivot_longer(c('total','core','shell'), names_to = 'region')
glimpse(mparts)

ggplot(total %>% filter(crosstype=='Abs'), aes(wavelength, average, linetype=calculation)) +
  geom_line(lwd=1, alpha=0.5) +
  geom_line(mparts, map=aes(wavelength, value, colour=region)) +
  # facet_grid(scatterer~crosstype, scales='free_y') +
  geom_line(data=ds, map=aes(wavelength, value, colour=region)) +
  facet_wrap(~scatterer, scales='free_y') +
  scale_colour_brewer(palette='Set1')

