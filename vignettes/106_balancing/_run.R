library(here)
setwd(here('vignettes/106_balancing'))


## ----load----
library(terms)
theme_set(theme_grey())

# ----tpl----

tpl <- "ModeAndScheme 2 2
MultipoleCutoff {Lmax}
Wavelength 300 900 100
Verbosity 1
Medium 1.7689 

Incidence 0 0 0 1
OutputFormat HDF5 {out}
{comment}DisableStoutBalancing

Scatterers 2
Ag 0.0  0.0  0.0 50.0
Ag 0.0  101.0 0.0 50.0
"


cat(tpl)

## ----run----
comment <- ''
Lmax <- 20
out <- 'xsec_1'
cat(glue(tpl), file = 'input1')
system("../../build/terms input1 > log1")

comment <- '#'
out <- 'xsec_2'
cat(glue(tpl), file = 'input2')
system("../../build/terms input2 > log2")


## ----read----

xs1 <- consolidate_xsec('xsec_1.h5')
xs2 <- consolidate_xsec('xsec_2.h5')

## ----oa----
mCOA <- rbind(cbind(xs1$mCOA, balancing = 'no balancing'),
              cbind(xs2$mCOA, balancing = 'balancing'))

mCOAt <- subset(mCOA, variable == 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average, linetype=balancing, colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line() +
  coord_cartesian(ylim=c(0, 9e4)) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = '', linetype='') +
  ggtitle("Orientation-averaged cross-sections")

p1


