setwd(here::here('test/epsilon/'))
library(terms)
library(dielectric)
library(glue)

d <- read.table('ph5004237_si_002.txt')
str(d)
# Wavelength(nm) Ag(e1) Ag(e2) Au(e1) Au(e2) Cu(e1) Cu(e2)

Ag <- data.frame(wavelength = d[,1], real = d[,2], imag = d[,3])
Au <- data.frame(wavelength = d[,1], real = d[,4], imag = d[,5])
Cu <- data.frame(wavelength = d[,1], real = d[,6], imag = d[,7])

write.table(Ag, "epsilon_Ag.txt", row.names = FALSE, col.names = FALSE)
write.table(Au, "epsilon_Au.txt", row.names = FALSE, col.names = FALSE)
write.table(Cu, "epsilon_Cu.txt", row.names = FALSE, col.names = FALSE)
