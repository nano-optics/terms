setwd(here::here('vignettes/13_masking'))

## ----load----
library(terms)
library(patchwork)
# library(reshape2)
# library(purrr)
# library(ggplot2)
# library(egg)
# library(grid)
# library(readr)


## ----run----
# unlink("*.dat")
# system("../../build/terms input > log1")
# system("../../build/terms input_EE1 > log2")
# system("../../build/terms input_EE2_MM0 > log3")

## ----staged, message=FALSE, echo=1:2----
af1 <- 'all_orders/prestagedA.txt'
af2 <- 'all_orders/stagedA_bal.txt'

display_amat(read_amat(af1)) + labs(title='pre-staged') +
  display_amat(read_amat(af2)) + labs(title='staged')

## ----collective, message=FALSE, echo=1:2----

f <- "tmat_Au10x20_Nmax4_lambda650.tmat"
f2 <- "all_orders/collective_all_orders.tmat"

v <- read_tmat(f)
v2 <- read_tmat(f2)
# glimpse(v2)

s <- subset(v, wavelength == 650 & n <=3 & np <=3)
s2 <- subset(v2, wavelength == 650 & n <=3 & np <=3)

display_tmat(s) + ggtitle('spheroid T-matrix') +
  display_tmat(s2) + ggtitle('dimer T-matrix')

## ----staged2, message=FALSE, echo=1:2----
af1 <- 'EE1/prestagedA.txt'
af2 <- 'EE1/stagedA_bal.txt'

# read.table(af1, header = F, nrows = 1)[2]


egg::ggarrange(display_amat(read_amat(af1)) + labs(title='pre-staged'),
               display_amat(read_amat(af2)) + labs(title='staged'), nrow=1)

## ----collective2, message=FALSE, echo=1:2----

f <- "tmat_Au10x20_Nmax4_lambda650.tmat"
f2 <- "EE1/collective_tmat_EE1.tmat"

v <- read_tmat(f)
v2 <- read_tmat(f2)
# glimpse(v2)

s <- subset(v, wavelength == 650 & n <=3 & np <=3)
s2 <- subset(v2, wavelength == 650 & n <=3 & np <=3)

display_tmat(s) + ggtitle('spheroid T-matrix') +
  display_tmat(s2) + ggtitle('dimer T-matrix')



## ----staged3, message=FALSE, echo=1:2----
af1 <- 'EE2_MM0/prestagedA.txt'
af2 <- 'EE2_MM0/stagedA_bal.txt'

display_amat(read_amat(af1)) + labs(title='pre-staged') +
  display_amat(read_amat(af2)) + labs(title='staged')

## ----collective3, message=FALSE, echo=1:2----

f <- "tmat_Au10x20_Nmax4_lambda650.tmat"
f2 <- "EE2_MM0/collective_tmat_EE2_MM0.tmat"

v <- read_tmat(f)
v2 <- read_tmat(f2)
# glimpse(v2)

s <- subset(v, wavelength == 650 & n <=3 & np <=3)
s2 <- subset(v2, wavelength == 650 & n <=3 & np <=3)

display_tmat(s) + ggtitle('spheroid T-matrix') +
  display_tmat(s2) + ggtitle('dimer T-matrix')


