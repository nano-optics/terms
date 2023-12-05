library(rhdf5)
library(glue)
library(R.matlab)
library(dplyr)
library(arrow)

setwd("~/Documents/nano-optics/terms/vignettes/1004_scuff_tmatrix")
d <- read.table('../../test/tmatrix/tmat_Au20x50_Nmax4.tmat')
# s sp n np m mp Tr Ti | a= 20 c= 50
# lambda= 400 nelements= 136 epsIn= -1.649657+5.771763j
names(d) <- c('s', 'sp', 'n', 'np', 'm', 'mp', 'Tr', 'Ti')
str(d)

## subsetting
d |> filter(n < 2, np < 2)



## transformation

d$prefactor <- with(d, (-1)^m * sqrt((2*n + 1) / (4*pi*n*(n+1))))
head(d)

## save as arrow
write_feather(d, 'long.arrow') # 346kb, from 1.9Mb in plain text


