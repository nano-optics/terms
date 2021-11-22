
## ----load---
library(terms)
library(glue)
library(purrr)
library(ggplot2)
library(egg)
theme_set(egg::theme_article())


## ----tpl---

tpl <- "ModeAndScheme 2 3
Wavelength 400 800 200
Medium 1.7689 # epsilon of water
MultipoleSelections 1
EE1:4_EM1:4_ME1:4_MM1:4 blocks
OutputFormat HDF5 cross_sections_{n1}_{n2}

MultipoleCutoff {n1} {n2}

Scatterers 4
Au_S1 200 0 -75 25
Au_S1 161.803398874989 117.557050458495 -25 25
Au_S1 61.8033988749895 190.211303259031 25 25
Au_S1 -61.8033988749895 190.211303259031 75 25
"

cat(tpl)

params <- crossing(n1 = 4:8, n2=4:10) %>% filter(n2 >= n1 & n2 <= n1+2)

## ----run---

for(ii in 1:nrow(params)){
  cat(glue_data(params[ii,], tpl),'\n', file=glue_data(params[ii,], "input_{n1}_{n2}"))
  system(glue_data(params[ii,], "../../build/terms input_{n1}_{n2} > log_{n1}_{n2}"))
}

## ----read----

xs <- purrr::pmap_df(params, function(n1, n2) terms::consolidate_xsec(glue('cross_sections_{n1}_{n2}.h5'))$mCOA, .id = 'id')
params$id <- as.character(1:nrow(params))
xs <- left_join(params, xs, by='id')
glimpse(xs)

## ----oa----

mCOAt <- subset(xs, variable == 'total')
mCOAn <- subset(xs, variable != 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average, colour=factor(n2), linetype=factor(n1))) +
  # facet_wrap(~ crosstype + abs(Scheme),scales = 'free',ncol=3)+
  facet_wrap(~crosstype)+
  # geom_path(data=mCOAn, aes(colour=variable)) +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[2]),
       linetype = expression(N[1])) +
  ggtitle("Orientation-averaged cross-sections")


p2 <- ggplot(mCOAt, aes(wavelength, dichroism, colour=factor(n2), linetype=factor(n1))) +
  facet_wrap(~crosstype)+
  # facet_wrap(~ crosstype + abs(Scheme),scales = 'free',ncol=3)+
  # geom_path(data=mCOAn, aes(colour=variable)) +
  geom_line() +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(N[2]),
       linetype = expression(N[1])) +
  ggtitle("Orientation-averaged circular dichroism")

library(patchwork)
p1 + p2 + plot_layout(ncol=1, guides = 'collect')


