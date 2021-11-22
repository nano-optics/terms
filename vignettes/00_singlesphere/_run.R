setwd(here::here("vignettes/00_singlesphere"))
## ----load----
library(terms)
library(dplyr)
theme_set(egg::theme_article())

## ----run----
system("../../build/terms input > log")

## ----read----
# grab the results and store in compressed form
xs <- consolidate_xsec('results.h5')


## ----oa----
p1 <- ggplot(xs$mCOA %>% filter(variable=='total'), aes(wavelength, average, colour=crosstype)) +
  facet_wrap(~crosstype) +
  guides(colour='none') +
  scale_colour_brewer(palette='Set1') +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Orientation-averaged cross-sections")

p1


