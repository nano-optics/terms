library(here)
setwd(here('vignettes/105_factorised_RTR'))

## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())

## ----tpl----

tpl <- "ModeAndScheme 2 3
Verbosity 3
MultipoleCutoff 8
Wavelength 300 700 100
Medium 1.7689 # epsilon of water
OutputFormat HDF5 {file}
"

cat(tpl)

## ----run----

N <- 10
file <- 'xsec_1'
cat(glue(tpl), file = 'input1', append = FALSE)
cat("\nScatterers ", N, "\n", file = 'input1', append = TRUE)
cluster_positions(N=N, 
                  cl_fun=cluster_chain,
                  pitch = 50, rot = rotation_euler_passive(0,pi/2,0),
                  radius = 20, label="Au", out='input1', digits=8)


system("../../build/terms input1 > log1")

file <- 'xsec_2'
cat(glue(tpl), file = 'input2', append = FALSE)
cat("\nDisableRTR \n", file = 'input2', append = TRUE)
cat("\nScatterers ", N, "\n", file = 'input2', append = TRUE)
cluster_positions(N=N, 
                  cl_fun=cluster_chain,
                  pitch = 50, rot = rotation_euler_passive(0,pi/2,0),
                  radius = 20, label="Au", out='input2', digits=8)

system("../../build/terms input2 > log2")


## ----read----

xs1 <- consolidate_xsec('xsec_1.h5')
xs2 <- consolidate_xsec('xsec_2.h5')

r1 <- extract_time('log1')
r2 <- extract_time('log2')


## ----oa----

mCOA <- rbind(cbind(xs1$mCOA, RTR = 'yes'),
              cbind(xs2$mCOA, RTR = 'no'))

mCOAt <- subset(mCOA, variable == 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average, linetype=RTR,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = '') +
  scale_colour_brewer(palette='Set1') +
  ggtitle("Orientation-averaged cross-sections")


p1

