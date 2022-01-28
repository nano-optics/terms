#setwd("~/home/Atefeh/New_TERMS/tex_pdf_files/Figs")
library(egg)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tikzDevice)
library(tinytex)
a <- setNames(read.table('/home/Atefeh/New_TERMS/tex_pdf_files/Figs/Fig1/orAveNOC_singlepoint.txt')[,c(1,2)],c('wavelength','c'))
a$material <- 'Orientation averaged, LCP'; a$k <- '$\\lldoc$'
b <- setNames(read.table('/home/Atefeh/New_TERMS/tex_pdf_files/Figs/Fig1/N_OC_pointA_LCP.txt')[,c(1,5)],c('wavelength','c'))
b$material <- 'NOC, LCP'; b$k <-  '$k_z$, LCP' #expression(k_z)#'$$, LCP'
c <- setNames(read.table('/home/Atefeh/New_TERMS/tex_pdf_files/Figs/Fig1/N_OC_pointA_RCP.txt')[,c(1,5)],c('wavelength','c'))
c$material <- 'NOC, RCP'; c$k <- '$k_z$, RCP'
#d <- setNames(read.table('/home/Atefeh/New_TERMS/tex_pdf_files/Figs/Fig1/N_OC_pointB_LCP.txt')[,c(1,5)],c('wavelength','x','y','z','c'))
#e <- setNames(read.table('/home/Atefeh/New_TERMS/tex_pdf_files/Figs/Fig1/N_OC_pointB_RCP.txt')[,c(1,5)],c('wavelength','x','y','z','c'))
f <- a
f$c <- -a$c
f$material <- 'Orientation averaged, RCP'; f$k <- '$\\lldoc$'
data_a <- rbind(a,f, b, c)
pa <- ggplot(data_a, aes( wavelength, c, color= material, linetype=k))+
  geom_line()+
  scale_x_continuous(expand=c(0,0), lim=c(400, 1500))+
  scale_y_continuous(expand=c(0,0),lim=c(-22,20))+
  labs(x='wavelength/ nm', y= 'Normalized optical chirality', map='C')+
  scale_colour_brewer(palette='Set1') +
  theme_grey(10) +
  theme()
  print(pa)
  fig <- pa
  ggsave('/home/Atefeh/New_TERMS/tex_pdf_files/fig1.pdf',fig, width=8, height=4.5)

