setwd(here::here("vignettes/110_alpha_tensor/"))
## ----load----
suppressPackageStartupMessages(require(terms))
library(egg)
# library(reshape2)
# library(purrr)
library(ggplot2)
library(egg)
# library(grid)
# library(readr)


## ----run----
# system("../../build/terms input_EE1 > log")

## ----collective, message=FALSE, echo=1:2----

f <- "tmat_Au20x50_Nmax4_lambda520.tmat"
f2 <- "collective_tmat_EE1.tmat"
f3 <- "alpha_col.txt"

v <- read_tmat(f)
v2 <- read_tmat(f2)
v3 <- read_tmat(f3)
# glimpse(v2)

# s <- subset(v,  n <=2 & np <=2)
# s2 <- subset(v2, n <=2 & np <=2)
# s3 <- subset(v3, n <=2 & np <=2)

s <- subset(v,   n <=3 & np <=3)
s2 <- subset(v2, n <=3 & np <=3)
s3 <- subset(v3, n <=3 & np <=3)

display_tmat(s) + ggtitle('spheroid T-matrix') +
  display_tmat(s2) + ggtitle('dimer T-matrix') +
  display_tmat(s3) + ggtitle('dimer alpha-tensor')


