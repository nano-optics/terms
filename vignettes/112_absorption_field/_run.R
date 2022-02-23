setwd(here::here('vignettes/112_absorption_field'))

## ----load----
library(terms)
theme_set(theme_grey())

## ----runff----

# system("../../build/terms input > log")

## ----total----

library(rhdf5)
lf <- h5ls('cross_sections.h5')
# lf

ld <- h5read('cross_sections.h5', "Far-Field")
# glimpse(ld)

# total OA cross-sections
xsec <- consolidate_xsec('cross_sections.h5')

total <- xsec$mCFO %>% filter(variable=='total')
total$scatterer <- 'both'
total$region <- 'total_avg'
total$calculation <- 'total'

# glimpse(total)

p <- ggplot(total, aes(wavelength, average, colour=crosstype)) +
  geom_line() +
  facet_wrap(~crosstype) +
  # scale_x_continuous(expand=c(0,0)) +
  scale_colour_brewer(palette='Set1') +
  labs(title='Far-field cross-sections', colour='',
       x=expression(wavelength/nm),y=expression(sigma/nm^2)) +
  theme(legend.position = 'none')

  
p


## ----partial----

# Mackowski's partial shell absorptions
mack <- consolidate_partials('cross_sections.h5')
# glimpse(mack)
mack$crosstype <- 'Abs'
mack$calculation <- 'Mackowski'
mack <- mack %>% mutate(core=partial_0, total = partial_1, shell=partial_1-partial_0) %>% 
  filter(polarisation=='4L') %>% 
  pivot_longer(c('total','core','shell'), names_to = 'region')
# glimpse(mack)

mack_both <- mack %>% pivot_wider(names_from = scatterer, id_cols = c(polarisation, wavelength,calculation, region), values_from = value) %>% 
  mutate(both = `1` + `2`)

ggplot(total %>% filter(crosstype=='Abs', scatterer == 'both'), 
       aes(wavelength, polarisation1)) +
  geom_area(aes(linetype=calculation, lwd='total'),alpha=0.5, fill='grey98') +
  geom_line(mack_both %>% filter(region == 'total'), 
            map=aes(wavelength, both, lty=region)) +
  geom_line(mack, map=aes(wavelength, value,lty=region, colour=factor(scatterer))) +
  # facet_wrap(~scatterer, scales='free_y', labeller = label_both) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),lim=c(0,12000)) +
  scale_colour_brewer(palette='Set1') +
  scale_fill_brewer(palette='Pastel1') +
  scale_size_manual(values=2)+
  labs(title='Absorption, divided by regions', colour='scatterer',linetype='', lwd='',
       x=expression(wavelength/nm),y=expression(sigma[abs]/nm^2))


## ----volumes----
volume_cubature <- function(rmin, rmax, r0 = c(0,0,0), 
                            Nr=5, Ntp=10, method='sphericaldesigns', out=NULL){
  
  cubature <- cubs::cubs(N=Ntp, cubature = method)
  rquad <- statmod::gauss.quad(Nr) #[-1,1]
  thickness <- rmax - rmin
  rnodes <- thickness * (rquad$nodes + 1 )/ 2 + rmin
  rweights <- thickness / 2 * rquad$weights * rnodes^2
  # quadrature <- data.frame(nodes=rnodes, weights = rweights)
  
  ptr <- expand.grid(phi = cubature[,1],
                     r = rnodes)
  ptr$theta <- rep(cubature[,2], Nr)
  
  positions_volume <- transmute(ptr,
                                x=r*sin(theta)*cos(phi) + r0[1],
                                y=r*sin(theta)*sin(phi) + r0[2],
                                z=r*cos(theta) + r0[3]
  )
  
  if(!is.null(out)){ # export directly to positions file
    cat(nrow(positions_volume), "\n", file = out)
    write.table(format(positions_volume, mode='double'),
                file = out, append = TRUE,
                col.names = FALSE, row.names = FALSE,
                quote = FALSE)
    
  }
  # return all information for post-processing
  invisible(list(positions = positions_volume,
                 ptr=ptr, rweights=rweights, tpweights=cubature[,3]))
}

volumes_list <- list(core1 = volume_cubature(rmin=0, rmax=18, r0 = c(0,0,0), Nr=15, Ntp=36),
                     shell1 = volume_cubature(rmin=18, rmax=21, r0 = c(0,0,0), Nr=15, Ntp=36),
                     core2 = volume_cubature(rmin=0, rmax=30, r0 = c(55,0,0), Nr=15, Ntp=36),
                     shell2 = volume_cubature(rmin=30, rmax=32, r0 = c(55,0,0), Nr=15, Ntp=36))

positions_list <- lapply(volumes_list, '[[', "positions")
# glimpse(volumes_list, max.level = 3)
# glimpse(positions_list)

positions <- do.call(rbind, positions_list)
pos_file <- 'positions'
cat(nrow(positions), "\n", file = pos_file)
write.table(format(positions, mode='double'),
            file = pos_file, append = TRUE,
            col.names = FALSE, row.names = FALSE,
            quote = FALSE)

# look at one set of points
ptr <- volumes_list$shell2$ptr

ggplot(ptr %>% mutate(r=factor(round(r,3))), 
       aes(phi, theta)) +
  geom_point() + #
  # theme_void() + 
  coord_equal() +
  facet_wrap(~r, labeller = label_both) +
  theme(panel.grid.major = element_line(colour = 'grey30', 
                                        linetype =1, size = 0.1),
        plot.margin = margin(4,4,4,4))+
  scale_x_continuous(expression(varphi),expand=c(0,0),lim=c(-pi,pi),breaks=seq(-2,2,by=1/2)*pi,
                     labels = c('-2π','-3π/2','-π','-π/2','0','π/2','π','3π/2','2π'))+
  scale_y_continuous(expression(theta),expand=c(0,0),lim=c(0,pi),breaks=seq(0,1,by=1/2)*pi,
                     labels = c('0','π/2','π'))+
  guides(size='none', colour='none') +
  labs(title='Sampled points',
       subtitle = 'across the shell thickness of scatterer 1')

## ----runnf----

# system("../../build/terms input_map > log_map")

## ----readnf----

dm <- rhdf5::h5read('map.h5', 'Near-Field')
names_MapQuantity <- c('wavelength','x','y','z', 'ScatID', 'volID', 'I')
ddv <- setNames(data.frame(dm$map_E[,1:7]), names_MapQuantity)
# str(ddv)
# unique(ddv$ScatID)
# unique(ddv$volID)

dl <- split(ddv, interaction(ddv$volID,ddv$ScatID))

## ----processnf----

ptr <- volumes_list$shell2$ptr
map <- dl[[2]]
Nlambda <- length(unique(map$wavelength))
all_pos <- ptr[rep(seq_len(nrow(ptr)), Nlambda),]
test <- cbind(all_pos, map) 
glimpse(test)

munchy <- data.frame(rbind(cbind(phi=rep(-90,10), 
                                 theta=seq(-90,90, length=10)),
                           cbind(phi=seq(-90,90, length=10), 
                                 theta=rep(90,10)),
                           cbind(phi=rep(90,10), 
                                 theta=seq(90,-90, length=10)),
                           cbind(phi=seq(90,-90, length=10), 
                                 theta=rep(-90,10))))
ggplot(test %>% filter(wavelength == 550) %>% mutate(r=factor(round(r,3))), 
       aes(phi*180/pi, 90-theta*180/pi, 
           colour=log10(I))) +
  geom_polygon(data=munchy, map=aes(phi, theta),fill='grey90',
               alpha=0.3, colour = NA, inherit.aes = FALSE) +
  geom_point() + #
  theme_void() + 
  facet_wrap(~r, labeller = label_both) +
  theme(panel.grid.major = element_line(colour = 'grey30', 
                                        linetype =1, size = 0.1),
        plot.margin = margin(4,4,4,4))+
  scale_x_continuous(breaks = seq(-180,180, by = 30),
                     minor_breaks = seq(-165,165, by = 30))+
  scale_y_continuous(breaks = seq(-180,180, by = 30),
                     minor_breaks = seq(-165,165, by = 30)) +
  coord_map("orthographic",orientation = c(0,0,0), 
            clip='off', xlim = c(-180, 180)) +
  scale_size(range = c(0.2, 1.6)) +
  scale_colour_viridis_c() +
  guides(size='none', colour='none') +
  labs(title='Sampled near-field points',
       subtitle = 'across the shell thickness of scatterer 1')

map_as_absorption <- function(map, ptr, rweights, tpweights,
                              n=1.33, fun_epsilon = dielectric::epsAu, ...){
  
  
  ptr$rf <- factor(ptr$r)
  map$lf <- factor(map$wavelength)
  lambda <- unique(map$wavelength)
  Nlambda <- length(lambda)
  
  # check size compatibility: map positions x number of wavelengths
  stopifnot(nrow(map) / Nlambda == nrow(ptr))
  
  # repeat positions for each wavelength and merge
  # note: not using left_join by x,y,z since those values 
  # are subject to numerical accuracy - with rbind we trust the order was maintained
  all_pos <- ptr[rep(seq_len(nrow(ptr)), Nlambda),]
  comb2 <- cbind(all_pos, map) 
  
  # integrate for a shell at constant r
  Ir <- comb2 %>% group_by(r, wavelength) %>%  
    summarise(integrand = 4*pi*sum(I*tpweights)) %>%
    ungroup() 
  
  # integrate all shells over radius
  integrated <- Ir %>% group_by(wavelength) %>% 
    summarise(integral =  sum(rweights*integrand)) %>% ungroup() %>% 
    mutate(absorption = 2*pi/wavelength/n*integral*Im(fun_epsilon(wavelength)$epsilon), 
           ...)
  
  integrated
}

# names(dl) # scatID, volID
core1 <- c(volumes_list[[1]], list(map = dl[[1]]))
shell1 <- c(volumes_list[[2]], list(map = dl[[2]]))
core2 <- c(volumes_list[[3]], list(map = dl[[3]]))
shell2 <- c(volumes_list[[4]], list(map = dl[[4]]))
# glimpse(shell1)

library(dielectric)
epsPd <- function(wavelength){
  
  raw_pd <- readRDS('pd.rds')
  
  epsr <- approx(raw_pd$wavelength,raw_pd$re_epsilon,wavelength)$y
  epsi <- approx(raw_pd$wavelength,raw_pd$im_epsilon,wavelength)$y
  data.frame(wavelength=wavelength, 
             epsilon = epsr + 1i*epsi)
  
}

partials <- list(with(core1, map_as_absorption(map, ptr, rweights, tpweights,
                                               n=1.33, fun_epsilon = dielectric::epsAu, 
                                               scatterer = '1', region='core')),
                 with(shell1, map_as_absorption(map, ptr, rweights, tpweights,
                                               n=1.33, fun_epsilon = dielectric::epsAu, 
                                               scatterer = '1', region='shell')),
                 with(core2, map_as_absorption(map, ptr, rweights, tpweights,
                                               n=1.33, fun_epsilon = dielectric::epsAu, 
                                               scatterer = '2', region='core')),
                 with(shell2, map_as_absorption(map, ptr, rweights, tpweights,
                                               n=1.33, fun_epsilon = epsPd, 
                                               scatterer = '2', region='shell')))


all_partials <- do.call(rbind, partials)

## ----comparison----

ggplot(mack, aes(wavelength, average, colour=region)) +
  # geom_line(mack_both, map=aes(wavelength, both, colour=calculation),lwd=1) +
  # geom_line(mack, map=aes(wavelength, core,lty=factor(scatterer)),lwd=1) +
  geom_line(mack, map=aes(wavelength, value, linetype=factor(scatterer))) +
  # geom_line(aes(linetype=calculation, colour=calculation), lwd=1) +
  geom_point(data=all_partials, aes(wavelength, absorption, colour=region,pch='integral(E^2)*dV')) +
  # facet_wrap(~scatterer, scales='free_y', labeller = label_both) +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_brewer(palette='Set1') +
  scale_shape_manual(values=21,labels=scales::label_parse())+
labs(title='Absorption, divided by regions', colour='region',
     linetype='scatterer', lwd='', shape='via near-field',
     x=expression(wavelength/nm),y=expression(sigma[abs]/nm^2))





