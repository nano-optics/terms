setwd(here::here('vignettes/109_individual_absorption'))

## ----load----
library(terms)
theme_set(egg::theme_article())

## ----tpl----


system("../../build/terms input > log")

## ----total----

library(rhdf5)
lf <- h5ls('cross_sections.h5')
# lf

ld <- h5read('cross_sections.h5', "Far-Field")
# glimpse(ld)

# total OA cross-sections
xsec <- consolidate_xsec('cross_sections.h5')

total <- xsec$mCOA %>% filter(variable=='total')
total$scatterer <- 'both'
total$region <- 'total_avg'
total$calculation <- 'total'

glimpse(total)

p <- ggplot(total, aes(wavelength, average, colour=crosstype)) +
  geom_line() +
  facet_wrap(~crosstype)

p


## ----partial----

# Mackowski's partial shell absorptions, with numerical cubature
xsec_part <- consolidate_partials('cross_sections.h5')
glimpse(ld$partial_absorption)

xsec_part$crosstype <- 'Abs'
xsec_part$calculation <- 'parts'
xsec_part <- xsec_part %>% mutate(total = partial_0)

# glimpse(xsec_part)

## ----split----

## Stout's per-particle OA absorptions

glimpse(ld$oa_incidence)

split <- setNames(data.frame(ld$Wavelengths, ld$oa_incidence$csAbsOA_split), c('wavelength','both','1','2'))
# glimpse(split)

split$crosstype <- 'Abs'

ds <- split %>% pivot_longer(c('both','1','2'))
ds$scatterer <- ds$name
ds$region <- 'total_avg'
ds$calculation <- 'split'

mparts <- xsec_part %>% filter(polarisation=='4L') %>% 
  pivot_longer(c('total'), names_to = 'region')

glimpse(mparts)


## ----comparison----

ggplot(total %>% filter(crosstype=='Abs'), aes(wavelength, average, linetype=calculation)) +
  geom_line(lwd=1, alpha=0.5) +
  geom_line(mparts, map=aes(wavelength, value, colour=region),lwd=1) +
  # facet_grid(scatterer~crosstype, scales='free_y') +
  geom_line(data=ds, map=aes(wavelength, value, colour=region),lwd=1) +
  facet_wrap(~scatterer, scales='free_y', labeller = label_both) +
  scale_colour_brewer(palette='Set1')

