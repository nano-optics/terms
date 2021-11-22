library(cubs)

q <- cubs(N = 18, cubature = 'lebedev')

cat(nrow(q), "\n", file = 'incidence')
write.table(format(cbind(q[,1],q[,2],0,q[,3]), mode='double'), file = 'incidence', 
            col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE)
