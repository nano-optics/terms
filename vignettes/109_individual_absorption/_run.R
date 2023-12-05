setwd(here::here('vignettes/109_individual_absorption'))

## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())

## ----tpl----


# system("../../build/terms input > log")

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
mack <- consolidate_partials('cross_sections.h5')

mack$crosstype <- 'Abs'
mack$calculation <- 'Mackowski'
mack <- mack %>% mutate(total = partial_0) %>% filter(polarisation=='4L') %>% 
  pivot_longer(c('total'), names_to = 'region')
glimpse(mack)

mack_both <- mack %>% pivot_wider(names_from = scatterer, id_cols = c(polarisation, wavelength,calculation, region), values_from = value) %>% 
  mutate(both = `1` + `2`)
# glimpse(mack_both)

# ggplot(mack_both, aes(wavelength, both)) + geom_line()

## ----split----

## Stout's per-particle OA absorptions

# glimpse(ld$oa_incidence)

split <- setNames(data.frame(ld$Wavelengths, ld$oa_incidence$csAbsOA_split), c('wavelength','both','1','2'))
# glimpse(split)

split$crosstype <- 'Abs'

stout <- split %>% pivot_longer(c('both','1','2'))
stout$scatterer <- stout$name
stout$region <- 'total_avg'
stout$calculation <- 'Stout'

glimpse(stout)

## ----byscatterer----

ggplot(map= aes(wavelength, average)) +
  geom_line(mack %>% filter(scatterer != 'both'), map=aes(wavelength, value, colour=calculation,linetype=calculation),lwd=1) +
  geom_line(data=stout %>% filter(scatterer != 'both'), map=aes(wavelength, value, colour=calculation,linetype=calculation),lwd=1) +
  facet_wrap(~scatterer, scales='free_y', labeller = label_both) +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_brewer(palette='Set1')+
  labs(title='Orientation-averaged absorption split by scatterer')


## ----both----

ggplot(total %>% filter(crosstype=='Abs', scatterer == 'both'), aes(wavelength, average)) +
  geom_line(mack_both, map=aes(wavelength, both, colour=calculation,linetype=calculation),lwd=1) +
  geom_line(data=stout %>% filter(scatterer == 'both'), map=aes(wavelength, value, colour=calculation,linetype=calculation),lwd=1) +
  geom_line(aes(linetype=calculation, colour=calculation), lwd=1) +
  # facet_wrap(~scatterer, scales='free_y', labeller = label_both) +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_brewer(palette='Set1') +
  labs(title='Total orientation-averaged absorption')

