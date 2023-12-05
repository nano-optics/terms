library(here)
setwd(here('vignettes/111_polarised_nf'))

## ----load----
suppressPackageStartupMessages(require(terms))
library(ggplot2)
theme_set(theme_grey())

tetramer <- function(radius=30, gap = 5, lab  = 'Au'){
  
  d <- 2*radius
  ## o    
  ## |
  ## o – o     
  ##      \
  ##       o   
  p1 <- c(0, d+gap, 0)
  p2 <- c(0, 0, 0)
  p3 <- c(d+gap, 0, 0)
  p4 <- c(d+gap, 0, d+gap)
  res <- cbind(lab, rbind(p1,p2,p3,p4), radius)
  write.table(res, quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(res)
}

# tetramer()
ge <- get_geometry('input_L')

## ----run----

system("../../build/terms input_u > log")
system("../../build/terms input_p > log")
system("../../build/terms input_L > log")


## ----readL----

library(rhdf5)
ld <- h5read('map_L.h5', "Near-Field")
d1 <- data.frame(ld$mapOaQuantity)
names(d1) <- c('wavelength', 'x','y','z','E^2', 'B^2', 'C^2')
glimpse(d1)

## ----readU----

ld <- h5read('map_unpolarised.h5', "Near-Field")
d2 <- data.frame(ld$mapOaQuantity)
names(d2) <- c('wavelength', 'x','y','z','E^2', 'B^2') # no C, because formula not available: use 1/2(L+R)
glimpse(d2)

## ----readP----

ld <- h5read('map_polarised.h5', "Near-Field")
d3 <- data.frame(ld$mapOaQuantity)
names(d3) <- c('wavelength', 'x','y','z','E^2[L]','E^2[R]', 'B^2[L]','B^2[R]','C[L]','C[R]')
glimpse(d3)

combined <- left_join(d2,d3,by=c('wavelength', 'x','y','z'))
combined <- combined %>% mutate(Ediff = (`E^2[L]` - `E^2[R]`),
                                Bdiff = (`B^2[L]` - `B^2[R]`),
                                Cdiff = (`C[L]` - `C[R]`)) #%>% pivot_longer(c(Ediff,Bdiff,Cdiff)) 

## ----plot----

pal <- RColorBrewer::brewer.pal(9,'PiYG')
pal2 <- RColorBrewer::brewer.pal(9,'PRGn')

p1 <- 
ggplot(combined, aes(x,y, fill=Ediff)) +
  geom_raster() + coord_equal() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, lwd=0.1, inherit.aes = FALSE) +
  scale_fill_gradient2(low = pal[1],mid = pal[5],high = pal[9], midpoint = 0) + 
  labs(x='x /nm', y='y /nm', fill=expression(E[L]^2-E[R]^2))

c <- 299792458
p2 <- p1 + aes(fill = c^2*Bdiff)+ labs(fill=expression(c^2*B[L]^2-c^2*B[R]^2))
p3 <- p1 + aes(fill = Cdiff)+ labs(fill=expression(C[L]-C[R]))
p4 <- 
  ggplot(combined, aes(x,y, fill=log10(`E^2`))) +
  geom_raster() + coord_equal() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  geom_circle(data=ge, aes(x0=x,y0=y,r=r), fill='grey90', lty=2, lwd=0.1, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = 'A') + labs(x='x /nm', y='y /nm', subtitle = 'Electric field intensity', fill=expression(log10(E^2)))
p5 <- p4 + scale_fill_viridis_c(option = 'D') + aes(fill = c^2*`B^2`)+ 
  labs(subtitle = 'Magnetic field intensity', fill=expression(c^2*B^2))
p6 <- p4 + aes(fill = `C[R]` + `C[L]`) + labs(subtitle = 'Degree of optical chirality', fill=expression(C[L]+C[R])) +
  scale_fill_gradient2(low = pal2[1],mid = pal2[5],high = pal2[9], midpoint = 0) 

egg::ggarrange(p4,p1,p5,p2,p6,p3,ncol=2)

