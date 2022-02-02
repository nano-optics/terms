setwd(here::here('test/incidence/'))
library(terms)
library(cubs)
library(glue)

for (N in lebedev_table$N[1:5]) terms::export_cubature(cubs(N, "lebedev"), out = glue::glue("incidence_leb_{N}"))
