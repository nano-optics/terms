setwd(here::here('vignettes/19_dispersive_medium/'))
## ----load----
library(terms)
library(dielectric)
theme_set(theme_grey())


d <- read.table('Ag_PRB.txt')
d
d <- d %>% mutate(V1 = V1/1.33, V2 = V2/1.33^2, V3 = V3/1.33^2)
write.table(d, file="Ag_scaled.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
## ----run----
system("../../build/terms input1 > log1")

system("../../build/terms input2 > log2")

## ----read----
xs1 <- terms::consolidate_xsec('water.h5')
xs2 <- terms::consolidate_xsec('vacuum.h5')
d2 <- xs2$mCOA %>% mutate(wavelength = wavelength*1.33)

## ----oa----
mCOA <- rbind(cbind(xs1$mCOA, epsilon = 'Ag in water'),
              cbind(d2, epsilon = 'Ag equivalent in vacuum'))

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
  
