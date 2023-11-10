setwd("~/Documents/nano-optics/terms/test/epsilon")
library(terms)
theme_set(theme_grey())

## ----run----
system("../../build/terms input_Au1 > log_Au1")

system("../../build/terms input_Au2 > log_Au2")

## ----read----
xs1 <- terms::consolidate_xsec('Au_default.h5')
xs2 <- terms::consolidate_xsec('Au2.h5')

## ----oa----
mCOA <- rbind(cbind(xs1$mCOA, epsilon = 'Au built-in'),
              cbind(xs2$mCOA, epsilon = 'Au Raschke'))

mCOAt <- subset(mCOA, variable == 'total')

p <- ggplot(mCOA, aes(wavelength, average, linetype = epsilon,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAt) +
  guides(colour='none') +
  scale_colour_brewer(palette='Set1') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N)) +
  ggtitle("Orientation-averaged cross-sections")

p
  
