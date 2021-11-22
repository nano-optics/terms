library(cubs)
data("lebedev")

write_cubature <- function(cubature, N, out=''){
  q <- cubs(N, cubature)
  cat(nrow(q),'\n',file=out)
  write.table(format(cbind(q[,1:2],0,q[,3]),digits = 15), file = out, quote=FALSE,
              row.names = FALSE,col.names = FALSE, append=TRUE )
}

for (N in lebedev_table$N[1:5]) write_cubature("lebedev", N, out=glue("incidence_leb_{N}"))
