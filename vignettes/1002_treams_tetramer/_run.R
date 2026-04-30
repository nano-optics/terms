setwd(here::here('vignettes/1002_treams_tetramer/'))
## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())
library(patchwork)


## ----run----
system("../../build/terms input_tetrahedron_external > log1")

system("../../build/terms input_tetrahedron_internal > log2")

## ----read----
xs1 <- terms::consolidate_xsec('tetrahedron_external.h5')
xs2 <- terms::consolidate_xsec('tetrahedron_internal.h5')

## ----plot----
mCOA <- rbind(cbind(xs1$mCOA, tmatrix = 'External'),
              cbind(xs2$mCOA, tmatrix = 'Internal')) 
mCOA$crosstype <- factor(mCOA$crosstype, levels = c("Ext", "Abs", "Sca"),
                         labels = c("Extinction", "Absorption", "Scattering")
)
mCOAt <- subset(mCOA, variable == 'total' & crosstype == "Extinction")

p1 <- ggplot(mCOA, aes(wavelength, average,
                       linetype = tmatrix,colour=tmatrix)) +
  geom_line(data=mCOAt) +
  scale_alpha_manual(values=c(0.8,1)) +
  guides(colour='none') +
  scale_x_continuous(expand=c(0,0))+ 
  scale_colour_brewer(palette='Set1') +
  # theme(legend.position	="inside", legend.position.inside = c(0.915,0.7))+
  labs(x = expression("wavelength /nm"), 
       y = expression("avg. extinction cross-sec. /"*nm^2),
       colour = "T-matrix", linetype="T-matrix", pch="",alpha="T-matrix") 

p1

