setwd("~/Documents/nano-optics/terms2/build")
## ----load----
library(terms)
library(egg)
# library(reshape2)
# library(purrr)
library(ggplot2)
# library(egg)
# library(grid)
# library(readr)



## ----collective, message=FALSE, echo=1:2----

f <- "tmat_Au20x50_Nmax4_lambda520.tmat"
f2 <- "collective_tmat_EE1.tmat"
f3 <- "alpha_col.txt"

v <- read_tmat(f)
v2 <- read_tmat(f2)
v3 <- read_tmat(f3)
# glimpse(v2)

s <- subset(v,  n <=2 & np <=2)
s2 <- subset(v2, n <=2 & np <=2)
s3 <- subset(v3, n <=2 & np <=2)

g <- egg::ggarrange(display_tmat(s) + ggtitle('spheroid T-matrix'),
display_tmat(s2) + ggtitle('dimer T-matrix'),
display_tmat(s3) + ggtitle('dimer alpha-tensor'),ncol=3)



