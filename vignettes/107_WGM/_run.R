library(here)
setwd(here('vignettes/107_WGM'))


## ----load----
library(terms)
theme_set(egg::theme_article())

# ----tpl----

tpl <- "ModeAndScheme 2 1
MultipoleCutoff 20
Wavelength 600 800 100
Verbosity 1
Medium 1.7689 

OutputFormat HDF5 {out}
{comment}DisableStoutBalancing

Scatterers 2
Si 0.0  0.0 0.0 2000.0
Si 0.0  4005 0.0 2000.0
"

cat(tpl)

## ----run----
comment <- ''
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
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = '', linetype='') +
  ggtitle("Orientation-averaged cross-sections")

p1


