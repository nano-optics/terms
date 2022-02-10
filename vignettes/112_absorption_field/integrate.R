setwd(here::here('vignettes/109_individual_absorption'))

## ----load----
library(terms)
theme_set(theme_grey())

## ----tpl----


system("../../build/terms input > log")

## ----total----

library(rhdf5)
lf <- h5ls('cross_sections.h5')
# lf

ld <- h5read('cross_sections.h5', "Far-Field")
# glimpse(ld)

# total OA cross-sections
xsec <- consolidate_xsec('cross_sections.h5')

total <- xsec$mCOA %>% filter(variable=='total')
total$scatterer <- 'both'
total$region <- 'total_avg'
total$calculation <- 'total'

glimpse(total)

p <- ggplot(total, aes(wavelength, average, colour=crosstype)) +
  geom_line() +
  facet_wrap(~crosstype)

p


## ----partial----

# Mackowski's partial shell absorptions, with numerical cubature
xsec_part <- consolidate_partials('cross_sections.h5')
glimpse(ld$partial_absorption)

xsec_part$crosstype <- 'Abs'
xsec_part$calculation <- 'parts'
xsec_part <- xsec_part %>% mutate(total = partial_0)

# glimpse(xsec_part)

## ----split----

## Stout's per-particle OA absorptions

glimpse(ld$oa_incidence)

split <- setNames(data.frame(ld$Wavelengths, ld$oa_incidence$csAbsOA_split), c('wavelength','both','1','2'))
# glimpse(split)

split$crosstype <- 'Abs'

ds <- split %>% pivot_longer(c('both','1','2'))
ds$scatterer <- ds$name
ds$region <- 'total_avg'
ds$calculation <- 'split'

mparts <- xsec_part %>% filter(polarisation=='4L') %>% 
  pivot_longer(c('total'), names_to = 'region')

glimpse(mparts)


## ----integration----

file_v <- 'positions_volume'
# centre of satellite to integrate
x0 <- 0
y0 <- 0
z0 <- 0

Nr <- 10
Ntp <- 30


# inc <- cubs::cubs(N=1500, cubature = 'lebedev')
# inc <- data.frame()

Nb <- 300
GL_cbeta <- statmod::gauss.quad(Nb)
beta = acos(GL_cbeta$nodes)

# grid of angles
nodes <- data.frame(alpha=0, beta=beta)
# corresponding weights for 2D cubature
weights <- data.frame(alpha=1, beta=GL_cbeta$weights)
# combine the weights and divide by 2 
# (1/4pi for the average, but * 2pi from range of alpha)
weights <- 1/2 * weights$alpha * weights$beta
inc <- setNames(data.frame(as.matrix(nodes), weights), c('phi','theta','weight'))
export_cubature(inc, out='incidence_fine')

Nr <- 10
Ntp <- 50
library(cubs)
quadrature <- cubs::cubs(N=Ntp, cubature = 'sphericaldesigns')
rquad <- statmod::gauss.quad(Nr) #[-1,1]
rmin <- 0.
rmax <- 20
thickness <- rmax - rmin
rnodes <- thickness * (rquad$nodes + 1 )/ 2 + rmin
rweights <- thickness / 2 * rquad$weights * rnodes^2

ptr <- expand.grid(phi = quadrature[,1],
                   r = rnodes)
ptr$theta <- rep(quadrature[,2], Nr)

positions_volume <- transmute(ptr,
                              x=r*sin(theta)*cos(phi) + x0,
                              y=r*sin(theta)*sin(phi) + y0,
                              z=r*cos(theta) + z0
)

cat(nrow(positions_volume), "\n", file = file_v)
write.table(format(positions_volume, mode='double'),
            file = file_v, append = TRUE,
            col.names = FALSE, row.names = FALSE,
            quote = FALSE)


system("../../build/terms input_integrate > log")
# 

dm <- rhdf5::h5read('map.h5', 'Near-Field')

ncol(dm$map_E)
names_MapQuantity <- c('wavelength','x','y','z', 'ScatID', 'volID', 'I')
ddv <- setNames(data.frame(dm$map_E[,1:7]), names_MapQuantity)

unique(ddv$ScatID)
unique(ddv$volID)
head(dm$map_E)
ptr$rf <- factor(ptr$r)
ddv$lf <- factor(ddv$wavelength)
levels(ptr$rf)
lambda <- unique(ddv$wavelength)
all_pos <- ptr[rep(seq_len(nrow(ptr)),length(lambda)),]
str(all_pos)
comb2 <- cbind(all_pos,ddv) 
str(comb2)

comb2 %>% group_by(rf,lf) %>% summarise(n())

s <- comb2 %>% filter(lf == ddv$lf[1])
unique(s$r)
unique(s$wavelength)
s  %>% group_by(r,wavelength) %>% summarise(n())

str(comb2)
tpweights <- quadrature[,3]
dr <- comb2 %>% group_by(r,wavelength) %>%  #group_by(rf,lf) %>% 
  summarise(integrand = 4*pi*sum(I*tpweights)) %>% #*sin(theta)
  ungroup() 
dr

# eps0 <- 8.854e-12
# c <- 299792458*1e9
n <- 1.33

integrated <- dr %>% group_by(wavelength) %>% 
  summarise(integral =  sum(rweights*integrand)) %>% ungroup() %>% 
  mutate(absorption = 2*pi/wavelength/n*integral*Im(dielectric::epsAu(wavelength)$epsilon),
         calculation = 'NF', scatterer=1)

# plot(integrated$wavelength,integrated$absorption)

## ----comparison----

ggplot(total %>% filter(crosstype=='Abs'), aes(wavelength, average)) +
  geom_line(aes(linetype=calculation), lwd=1, alpha=0.5) +
  geom_line(mparts, map=aes(wavelength, value, colour=region,linetype=calculation),lwd=1) +
  # geom_point(data=integrated, aes(wavelength, absorption,pch=calculation)) +
  # facet_grid(scatterer~crosstype, scales='free_y') +
  geom_line(data=ds, map=aes(wavelength, value, colour=region,linetype=calculation),lwd=1) +
  facet_wrap(~scatterer, scales='free_y', labeller = label_both) +
  scale_colour_brewer(palette='Set1')

