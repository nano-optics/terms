library(here)
setwd(here('vignettes/04_chiral_dimer_spheroids'))

## ----load----
library(terms)
theme_set(theme_grey())

## ----run----
system("../../build/terms input > log")

## ----read----
# grab the results
xs <- terms::consolidate_xsec('cross_sections.h5')

# plot OA data
## ----oa----

glimpse(xs$mCOA)
mCOAt <- subset(xs$mCOA, variable == 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAt, lty=1) +
  guides(colour='none')+
  scale_color_brewer(palette = 'Set1') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(type)) +
  ggtitle("Orientation-averaged cross-sections")

p2 <- ggplot(mCOAt, aes(wavelength, dichroism,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAt, lty=1) +
  guides(colour='none')+
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(limits = egg::symmetric_range) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(type)) +
  ggtitle("Orientation-averaged circular dichroism")

egg::ggarrange(p1, p2,ncol=1)


## ----cleanup----
# unlink("*.dat")

