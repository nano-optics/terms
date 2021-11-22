setwd(here::here('vignettes/06_dielectric_function/'))
## ----load----
library(terms)
theme_set(egg::theme_article())

## ----run----
system("../../build/terms input1 > log1")

system("../../build/terms input2 > log2")

## ----read----
xs1 <- terms::consolidate_xsec('Ag_default.h5')
xs2 <- terms::consolidate_xsec('Ag_PRB.h5')

## ----oa----
mCOA <- rbind(cbind(xs1$mCOA, epsilon = 'Ag built-in'),
              cbind(xs2$mCOA, epsilon = 'Ag PRB'))

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

