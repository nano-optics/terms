setwd(here::here("vignettes/10_shell_absorption"))

## ----load----
library(terms)
theme_set(egg::theme_article())

## ----run----
system("../../build/terms input_AuAg > log1")
system("../../build/terms input_AuAu > log2")

## ----xsec----

library(rhdf5)
lf <- h5ls('cross_sections_AuAg.h5')
lf

ld_AuAu <- h5read('cross_sections_AuAu.h5', "Far-Field")
ld_AuAg <- h5read('cross_sections_AuAg.h5', "Far-Field")


# grab the total results
xsec_AuAu <- consolidate_xsec('cross_sections_AuAu.h5')
xsec_AuAg <- consolidate_xsec('cross_sections_AuAg.h5')

full <- rbind(cbind(xsec_AuAu$mCOA, case = "Au@Au"),
              cbind(xsec_AuAg$mCOA, case = "Au@Ag"))

total <- subset(full, variable == 'total')
glimpse(total)

## ----coats----

part_AuAu <- consolidate_partials('cross_sections_AuAu.h5')
part_AuAg <- consolidate_partials('cross_sections_AuAg.h5')

glimpse(part_AuAu)

process_partials <- function(part){
  part %>% mutate(core = partial_0, total = partial_1, shell = total - core) %>%
    select(-partial_0, -partial_1) %>%
    ungroup() %>% group_by(scatterer, wavelength) %>%
    filter(polarisation %in% c('3R', '4L')) %>% # take both circular
    summarise(core_avg = mean(core), shell_avg = mean(shell), total_avg = mean(total), .groups = 'drop') %>%
    pivot_longer(c('total_avg','shell_avg','core_avg'), names_to = 'region') %>%
    ungroup()
}


coats <- rbind(cbind(process_partials(part_AuAu), case = "Au@Au"),
               cbind(process_partials(part_AuAg), case = "Au@Ag"))

coats$crosstype <- 'Abs'
glimpse(coats)


## ----comparison----

all <- ggplot(coats %>% filter(scatterer==1), aes(wavelength, value, colour=crosstype)) +
  facet_grid(~case)+
  geom_line(aes(linetype=region)) +
  geom_line(data=total, aes(wavelength, average/2)) +
  scale_linetype_manual(values=c(3,2,1)) +
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(lim=c(0,19000),expand=c(0,0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x='wavelength /nm', y=expression("per-part. cross-section /"*nm^2),
       colour='',linetype = 'region') +
  theme_article() +
  theme(panel.grid.major.y = element_line(colour = 'grey80',size = 0.2,linetype=3),
        panel.grid.minor.y = element_line(colour = 'grey80',size = 0.1,linetype=3))

all


