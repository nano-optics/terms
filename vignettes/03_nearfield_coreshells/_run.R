setwd(here::here('vignettes/03_nearfield_coreshells/'))
## ----load----
library(terms)
library(rhdf5)
theme_set(theme_grey())


## ----run----
system('../../build/terms input_ff > log_ff')
system('../../build/terms input_nf > log_nf')
system('../../build/terms input_nfcloseup > log_nfcloseup')

## ----ff----
ff <- terms::consolidate_xsec('cross_sections.h5')
# ld <- h5read('ff/cross_sections.h5', "Far-Field")

m <- subset(ff$mCOA, variable == 'total')

ggplot(m, aes(wavelength, average, colour=crosstype)) +
  geom_line()+
  scale_colour_brewer(palette = 'Set1') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = '') +
  ggtitle("Orientation-averaged cross-sections")

## ----nf----
d <- h5read('map.h5', "Near-Field")
map <- data.frame(d$map_E)
names(map) <- c('lambda', 'x', 'y', 'z', 'scatID', 'volID',  'E2avg', 'E2X', 'E2Y')
glimpse(map)

geometry <- get_geometry('input_nf')

library(tidyr)
m <- pivot_longer(map, c(E2X,E2Y))

p <- ggplot(m, aes(x,y)) +
  geom_raster(aes(fill=log10(value))) +
  facet_wrap(~name, ncol=2) +
  coord_equal() +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.2,lty=2) +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=alpha,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.2,lty=2) +
  scale_fill_viridis_c(option = 'plasma') +
  theme_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="x /nm", y="y /nm", fill="log(I)", 
       title="Near-field intensity map at 680 nm",
       subtitle = "incidence along z; x and y polarisations")
p

## ----nfcloseup----

d <- h5read('map_closeup.h5', "Near-Field")
map <- data.frame(d$map_E)
names(map) <- c('lambda', 'x', 'y', 'z', 'scatID', 'volID', 'E2')

indexed_epsilon <- function(volID, scatID, wavelength){
  case_when(scatID != 0 & volID == -1 ~ Im(dielectric::epsAu(wavelength)$epsilon),
            scatID != 0 & volID == 0 ~ Im(dielectric::epsAg(wavelength)$epsilon),
            # scatID ==0 ~ NA_real_,
            TRUE ~ 0)
}

map <- map %>% mutate(absorption = E2 * indexed_epsilon(volID, scatID, lambda))
glimpse(map)
geometry <- get_geometry('input_nfcloseup')

library(tidyr)
p <- ggplot(map, aes(x,y)) +
  geom_raster(aes(fill=(absorption))) +
  coord_equal(xlim=c(-10,10),ylim=c(-10,10)) +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=r,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.2,lty=2) +
  geom_circle(data=geometry, aes(x0=x,y0=y,r=alpha,group=label),
              inherit.aes=FALSE,
              colour='white', alpha=0.5, lwd=0.2,lty=2) +
  scale_fill_viridis_c(option = 'plasma') +
  theme_grey() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="x /nm", y="y /nm", fill=expression(Im(epsilon)*"|E|"^2), 
       title="Local absorption map at 680 nm",
       subtitle = "incidence along z; x polarisation")
p

# ggsave('nf.pdf',p,width=8,height=4)
