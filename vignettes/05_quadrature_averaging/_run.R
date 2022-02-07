setwd(here::here('vignettes/05_quadrature_averaging/'))
## ----load----
library(terms)
theme_set(theme_grey())

## ----run----
system("../../build/terms input > log")

## ----read----
# grab the results 
xs <- terms::consolidate_xsec('cross_sections.h5')

mCOAt <- subset(xs$mCOA, variable == 'total')

mCFOt <- subset(xs$mCFO, variable == 'total')


## ----comparison----
p7 <- ggplot(mCFOt, aes(wavelength, average, colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line(lty=1) +
  geom_line(data=mCOAt, lty=2) +
  guides(colour='none')+
  scale_color_brewer(palette = 'Set1') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Orientation-averaged circular polarisation cross-sections")

p8 <- ggplot(mCFOt, aes(wavelength, dichroism, colour=crosstype)) +
  facet_wrap(~crosstype)+
  # annotate('hline', yintercept=0,y=0,lty=3) +
  geom_line(lty=1) +
  geom_line(data=mCOAt, lty=2) +
  guides(colour='none')+
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(limits = symmetric_range) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Orientation-averaged circular dichroism")


egg::ggarrange(p7,p8,ncol=1)

