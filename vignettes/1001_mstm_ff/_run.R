setwd(here::here("vignettes/1001_mstm_ff"))
## ----load----
library(terms)
library(dplyr)
library(ggplot2)
library(dplyr)
library(glue)
library(purrr)
library(tidyr)

theme_set(theme_grey())

## ----run----
# system("../../build/terms input_EE1 > log")


## ----mstm----

read_one <- function(f = '650.dat', N=8){
  
  m <- read.table(pipe(glue('grep -A{N+1} " sphere unpolarized extinction" {f}')), skip=1, header = T)
  d <- setNames(data.frame(t(colSums(as.matrix(m[,-c(1,4)])))), c('Ext', 'Abs'))
  d$Sca <- d$Ext - d$Abs
  d$wavelength <- as.numeric(gsub('.dat', '', f))
  d
}

# read_one()

d <- map_df(glue('{seq(400,800,by=50)}.dat'), read_one)

# saveRDS(d, 'mstm.rds')


## ----comparison----

d <- readRDS('mstm.rds')
terms <- consolidate_xsec('terms.h5')


ggplot(d %>% pivot_longer(-wavelength), aes(wavelength, value*pi*50^2/1e6,colour=name)) +
  geom_point(aes(pch='MSTM')) +
  geom_line(data=terms$mCOA %>% filter(variable=='total'), aes(wavelength, average/1e6, colour=crosstype)) +
  scale_colour_brewer(palette='Set1') +
  scale_x_continuous(lim=c(380,820),expand=c(0,0)) +
  guides(color = guide_legend(override.aes = list(shape=NA) ) ) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*mu*m^2),
       colour = '', pch='') 

