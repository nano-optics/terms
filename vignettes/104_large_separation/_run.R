library(here)
setwd(here('vignettes/104_large_separation'))

## ----load---
library(terms)
theme_set(egg::theme_article())

# ----tpl----

tpl <- "ModeAndScheme 2 {scheme}
MultipoleCutoff 8
Wavelength 300 700 100
Medium 1.7689 # epsilon of water
Incidence  0.0 0.0 0.0 
# only dipoles, similar outcome
# MultipoleSelections 1
# EE1:1_EM2:1_ME2:1_MM2:1  blocks
OutputFormat HDF5 xsec_{separation}_{comment_centred}_{scheme}

{comment}ScattererCentredCrossSections

Scatterers 4  
Au  0 0    {-3*separation/2}  30
Au  0 0    {-separation/2}  30
Au  0 0    {+separation/2}  30
Au  0 0    {+3*separation/2}  30
"

cat(tpl)

params <- crossing(separation = c(200, 500, 1000), comment_centred=c(0, 1), scheme=0:3)

## ----run----

for(ii in 1:nrow(params)){
  current <- params[ii,]
  current$comment <- if(current$comment_centred) '#' else ''
  input <- glue_data(current, "input_{separation}_{comment_centred}_{scheme}")
  log <- glue_data(current, "log_{separation}_{comment_centred}_{scheme}")
  cat(glue_data(current, tpl), file=input)
  if(!file.exists(glue_data(params[ii,], "xsec_{separation}_{comment_centred}_{scheme}.h5")))
  system(glue("../../build/terms {input} > {log} &"))
}


## ----read----

lf <- glue_data(params, "xsec_{separation}_{comment_centred}_{scheme}.h5")

m <- map_df(lf, function(f) consolidate_xsec(f)$mLFO, .id='file')
params$file <- as.character(seq(1,nrow(params)))
combined <- left_join(params, m, by='file')
glimpse(combined)

## ----fo----

ground_truth <- combined %>% filter(variable=='total', scheme == 0, !as.logical(comment_centred))
ground_truth$scheme <- NULL
combined$comment <- ifelse(!combined$comment_centred, "ScattererCentred", "ClusterCentred")

p <- ggplot(combined %>% filter(variable=='total'), 
            aes(wavelength, average, colour=factor(crosstype))) +
  geom_line(data = ground_truth, alpha=0.5,lwd=1.2, show.legend = c(TRUE,FALSE)) +
  geom_line(aes(linetype=comment)) +
  facet_grid(separation~scheme, scales='free_y')+
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = '', linetype = '') +
  scale_linetype_manual(values=c(2,1)) +
  ggtitle("Fixed-orientation linear polarisation cross-sections") +
  theme(legend.position = 'top', legend.direction = 'horizontal')

p

