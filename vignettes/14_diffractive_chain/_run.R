
## ----load---
suppressPackageStartupMessages(require(terms))
library(reshape2)
library(purrr)
library(ggplot2)
library(egg)
theme_set(theme_grey())

## ----cluster----

cl <- cluster_chain(20, pitch = 450)
radius <- 50.0
# tpl

tpl <- "ModeAndScheme 2 2
MultipoleCutoff 5
MultipoleSelections 1
EE1:1  blocks
OutputFormat HDF5 cross_sections_chain
ScattererCentredCrossSections
Wavelength 300 700 100
Medium 1.7689 # epsilon of water
Incidence  0.0 0.0 0.0 
Scatterers {N}"

tpl_ref <- "ModeAndScheme 2 2
MultipoleCutoff 5
MultipoleSelections 1
EE1:1  blocks
OutputFormat HDF5 cross_sections_ref
ScattererCentredCrossSections
Wavelength 300 700 100
Medium 1.7689 # epsilon of water
Incidence  0.0 0.0 0.0 
Scatterers 1
Au_S1 0   0     0     50  
"
cat(glue(tpl_ref), "\n", file='input_ref')

N <- ncol(cl$positions)
# scatterers <- 
scatterers <-  format(cbind("Au_S1", t(cl$positions), radius), mode = "double")
scatterers
cat(glue(tpl), "\n", file='input_chain')
write.table(scatterers, file = 'input_chain', append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)


## ----run----
system("../../build/terms input_ref > log_ref")
system("../../build/terms input_chain > log_chain")

## ----read----
xs <- terms::consolidate_xsec('cross_sections_chain.h5')
glimpse(xs)
ref <- terms::consolidate_xsec('cross_sections_ref.h5')


## ----fo----
glimpse(xs$mLFO %>% pivot_longer(c(polarisation1, polarisation2)))

m1 <- xs$mLFO %>% filter(variable=='total') %>% pivot_longer(c(polarisation1, polarisation2))
m2 <- ref$mLFO %>% filter(variable=='total') %>% pivot_longer(c(polarisation1, polarisation2))

p3 <- ggplot(m1, aes(wavelength, value)) +
  facet_grid(.~crosstype, scales='free_y')+
  geom_line(aes(linetype=name)) +
  geom_line(aes(y=value*20, colour='reference'), 
            data=m2) +
  # geom_path(aes(colour=variable), lty=3) +
  # geom_path() +
  scale_colour_manual(values=blues9[5])+
  scale_linetype(labels=c("x","y"))+
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       linetype = expression(polarisation),colour='',
  title="Fixed-orientation linear polarisation cross-sections",
  subtitle = 'linear chain of 20 spheres along x')

p3
# ggsave("plot_fo.pdf", p3, width=8, height=3)
# ggsave("plot_fo_stout.pdf", p3, width=6, height=4)
