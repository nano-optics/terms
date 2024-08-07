library(here)
setwd(here('vignettes/108_timing_nf/'))
getwd()

## ----params----

suppressPackageStartupMessages(require(terms))
library(dplyr)

# easy case, check scaling with Npart
par1 <- crossing(Npart = c(1, 2, 5, seq(10, 40, by = 10)),
                 Nmax =  c(8, 16),
                 Ninc = 1,
                 OA = TRUE)

# few particles, check scaling with Nmax
par2 <- crossing(Npart = c(3, 10),
                 Nmax =  seq(5, 20, by = 5),
                 Ninc = 1,
                 OA = TRUE)

# medium difficulty, check scaling with Ninc, no OA
par3 <- crossing(Npart = c(4, 8),
                 Nmax =  8,
                 Ninc = c(1, 50, seq(100, 300, by=100)),
                 OA = FALSE)

library(DT)
params <- crossing(rbind(par1, par2, par3), scheme=0:3) %>% arrange(Npart, -Ninc)
datatable(params, rownames = FALSE,
          options = list(pageLength = 5, scrollX=T) )

## ----tpl----

tpl <- "ModeAndScheme 1 {scheme} 
MultipoleCutoff {Nmax}
ScattererCentredCrossSections
Incidence {Ninc} 0.0 0.0 1 
OutputFormat HDF5 map_tmp
Wavelength 633
Medium 1.7689 # epsilon of water
Verbosity 2 # show timings

SpacePoints 32 0 0 0 0 0 0 0 0 # picking a point
MapQuantity 2 E C
{comment_oa}MapOaQuantity E C

Scatterers  {Npart}
"

input_file <- function(Npart=15, Ninc=1, Nmax=8, 
                       scheme = 0, OA = FALSE, 
                       radius = 30,
                       label="Au", ...){
  
  comment_oa <- if(OA) '' else '#'
  out <- glue("runs/input_{Npart}_{Ninc}_{Nmax}_{scheme}")
  c <- cluster_chain(N=Npart, pitch = 100, a = radius, b = radius, c = radius, 
                     rot = rotation_euler_passive(pi/4, pi/4, 0)) # arbitrary rotation
  
  if(Ninc == 1) Ninc <- 0 else if(Ninc > 1) Ninc <- paste0("-", Ninc) # discretise Ninc angles
  cat(glue(tpl), "\n", file = out, append = FALSE)
  write.table(format(cbind(label, t(c$positions), radius), mode='double'), 
              file = out, 
              col.names = FALSE, append = TRUE, 
              row.names = FALSE, quote = FALSE)
  
}


## ----run---

bashrun <- 'run.sh'

fs::dir_create('runs')  
cat('#!/bin/bash \n\n', file = bashrun, append = FALSE)
for(ii in 1:nrow(params)){
  
  do.call(input_file, params[ii,])
  
  input <- glue_data(params[ii,], 'runs/input_{Npart}_{Ninc}_{Nmax}_{scheme}')
  log <- glue_data(params[ii,], 'runs/log_{Npart}_{Ninc}_{Nmax}_{scheme}')
  cat(glue("../../build/terms {input} > {log}"),
      '\n', 'rm map_tmp.h5\n', sep = '', file = bashrun, append = TRUE)
  
}


## ----read---
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())

old_extract_time <- function(log = 'log'){
  command <- paste0("grep '(CPU & real in s)' ", log)
  timestrings <- readLines(pipe(command))
  splits <- strsplit(timestrings, ">")
  routine <- sapply(splits, "[", 1)
  action <- sapply(splits, "[", 2)
  stplits2 <- strsplit(timestrings, " ")
  cpu <- sapply(stplits2, function(v) as.numeric(v[length(v)-2]))
  total <- sapply(stplits2, function(v) as.numeric(v[length(v)]))
  
  tibble(routine = routine, cpu = cpu, total=total)
}



extract_time <- function(log = 'log'){
  command <- paste0("grep '(CPU & real in s)' ", log)
  timestrings <- readLines(pipe(command))
  calculation <- strsplit(timestrings, ">")
  routine <- sapply(calculation, "[", 1)
  detail <- stringr::str_extract(string = sapply(calculation, "[", 2),
                                 pattern = "(?<=\\[).*(?=\\])")
  
  cpu <- sapply(strsplit(timestrings, " "), 
                function(v) as.numeric(v[length(v)-2]))
  total <- sapply(strsplit(timestrings, " "), 
                  function(v) as.numeric(v[length(v)]))
  
  tibble(routine = routine, detail = detail, cpu = cpu, total=total)
}

if(!file.exists('timings.rda')){
  
  params2 <- params[1:100,]
  all_logs <- glue_data(params2,"runs/log_{Npart}_{Ninc}_{Nmax}_{scheme}")
  times <- map_df(all_logs, extract_time, .id='log')
  params2$log <- as.character(seq(1,nrow(params2)))
  
  combined <- inner_join(params2, times, by='log')
  saveRDS(combined,'timings.rda')
  
} else {
  combined <- readRDS('timings.rda')
}

d1 <- inner_join(par1, combined %>% filter(routine == 'termsProgram') %>% select(-detail), by = c("Npart", "Nmax", "Ninc"))

p1 <- ggplot(d1,  aes(Npart, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Nmax ~ routine, labeller = label_both, scales='free_y')+
  scale_y_log10(
    # breaks = 10^seq(-3,3,by=1),lim=10^c(-3,3.7),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2,10^3)
    labels=scales::label_number()
  ) + 
  scale_x_continuous(breaks = unique(d1$Npart)[-2]) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs number of particles',
       subtitle = 'Ninc = 1, Nmax = 8 and 16') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

# p1


# d2 <- subset(combined, Npart %in% c(3,10) & Ninc == 1 & !(Nmax %in% c(8,16)))
d2 <- inner_join(par2, combined %>% filter(routine == 'termsProgram') %>% select(-detail), by = c("Npart", "Nmax", "Ninc"))
p2 <- ggplot(d2, aes(Nmax, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Npart ~ routine, labeller = label_both, scales='free_y')+
  scale_y_log10(
    # breaks = 10^seq(-3,3,by=1),lim=10^c(-3,3.7),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2,10^3)
    labels=scales::label_number()
  ) + 
  scale_x_continuous(breaks = unique(d2$Nmax)) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs multipole order',
       subtitle = 'Ninc = 1, Npart = 3 and 10')  +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

# p2

d3 <- inner_join(par3, combined %>% filter(routine == 'termsProgram') %>% select(-detail), by = c("Npart", "Nmax", "Ninc"))
p3 <- ggplot(d3, aes(Ninc, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Npart ~ routine, labeller = label_both, scales='free_y')+
  scale_y_log10(
    labels=scales::label_number()
    # breaks = 10^seq(-3,2,by=1), lim=10^c(-3,2),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2)
  ) + 
  scale_x_continuous(breaks = sort(unique(d3$Ninc))[-2]) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs incident angles',
       subtitle = 'Nmax = 8, Npart = 4 and 8') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

# p3


## ----fulltimings---

smaller_titles <- theme(plot.title = element_text(size = 10), 
                        plot.subtitle = element_text(size = 8))
egg::ggarrange(p1 + smaller_titles, 
               p2 + smaller_titles, 
               p3  + smaller_titles, nrow=1)


## ----detailedtimings_Npart---

dd1 <- inner_join(par1, combined %>% filter(routine != 'termsProgram') %>% 
                    select(-detail), by = c("Npart", "Nmax", "Ninc"))

pp1 <- ggplot(dd1,  aes(Npart, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Nmax ~ routine, labeller = labeller(.rows = label_both), scales='free_y')+
  scale_y_log10(
    # breaks = 10^seq(-3,3,by=1),lim=10^c(-3,3.7),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2,10^3)
    labels=scales::label_number()
  ) + 
  scale_x_continuous(breaks = unique(d1$Npart)[-2]) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs number of particles',
       subtitle = 'Ninc = 1, Nmax = 8 and 16') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp1

## ----detailedtimings_Npart_2---

dd1p <- inner_join(par1, combined %>% filter(routine == 'solve', !is.na(detail)),
                   by = c("Npart", "Nmax", "Ninc"))

pp1p <- ggplot(dd1p,  aes(Npart, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Nmax ~ detail, labeller = labeller(.rows = label_both), scales='free_y')+
  scale_y_log10(
    # breaks = 10^seq(-3,3,by=1),lim=10^c(-3,3.7),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2,10^3)
    labels=scales::label_number()
  ) + 
  scale_x_continuous(breaks = unique(d1$Npart)[-2]) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs number of particles',
       subtitle = 'Ninc = 1, Nmax = 8 and 16') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp1p

## ----detailedtimings_Nmax---

dd2 <- inner_join(par2, combined %>% filter(routine != 'termsProgram') %>% select(-detail), by = c("Npart", "Nmax", "Ninc"))
pp2 <- ggplot(dd2, aes(Nmax, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Npart ~ routine, labeller = labeller(.rows = label_both), scales='free_y')+
  scale_y_log10(
    # breaks = 10^seq(-3,3,by=1),lim=10^c(-3,3.7),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2,10^3)
    labels=scales::label_number()
  ) + 
  scale_x_continuous(breaks = unique(d2$Nmax)) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs multipole order',
       subtitle = 'Ninc = 1, Npart = 3 and 10')  +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp2


## ----detailedtimings_Nmax_2---


dd2p <- inner_join(par2, combined %>% filter(routine == 'solve', !is.na(detail)), by = c("Npart", "Nmax", "Ninc"))
pp2p <- ggplot(dd2p, aes(Nmax, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Npart ~ detail, labeller = labeller(.rows = label_both), scales='free_y')+
  scale_y_log10(
    # breaks = 10^seq(-3,3,by=1),lim=10^c(-3,3.7),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2,10^3)
    labels=scales::label_number()
  ) + 
  scale_x_continuous(breaks = unique(d2$Nmax)) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs multipole order',
       subtitle = 'Ninc = 1, Npart = 3 and 10')  +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp2p


## ----detailedtimings_Ninc---

dd3 <- inner_join(par3, combined %>% filter(routine != 'termsProgram') %>% select(-detail), by = c("Npart", "Nmax", "Ninc"))
pp3 <- ggplot(dd3, aes(Ninc, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_grid(Npart ~ routine, labeller = labeller(.rows = label_both), scales='free_y')+
  scale_y_log10(
    labels=scales::label_number()
    # breaks = 10^seq(-3,2,by=1), lim=10^c(-3,2),
    # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2)
  ) + 
  scale_x_continuous(breaks = sort(unique(d3$Ninc))[-2]) +
  labs(colour='scheme', y='log(time) /s',
       title = 'Timing vs incident angles',
       subtitle = 'Nmax = 8, Npart = 4 and 8') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp3


## ----detailedtimings_Ninc_2---

dd3p <- inner_join(par3, combined %>% 
                     filter(routine == 'solve', !is.na(detail)), by = c("Npart", "Nmax", "Ninc"))
pp3p <- ggplot(dd3p, aes(Ninc, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_wrap(Npart ~ detail, nrow=2, labeller = labeller(.rows = label_both), scales='free_y')+
  # scale_y_log10(
  # labels=scales::label_number()
  # breaks = 10^seq(-3,2,by=1), lim=10^c(-3,2),
  # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2)
  # ) + 
  scale_x_continuous(breaks = sort(unique(d3$Ninc))[-2]) +
  labs(colour='scheme', y='time /s',
       title = 'Timing vs incident angles',
       subtitle = 'Nmax = 8, Npart = 4 and 8') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp3p


## ----detailedtimings_Ninc_3---

dd4p <- inner_join(par3, combined %>% filter(routine == 'mapNF', !is.na(detail)), by = c("Npart", "Nmax", "Ninc"))
pp4p <- ggplot(dd4p, aes(Ninc, total, colour=factor(scheme))) +
  geom_point() +
  geom_line(stat = "summary", fun.data=function(x) tibble(ymin=max(x),y=max(x),ymax=max(x))) +
  facet_wrap(Npart ~ detail, nrow=2, scales='free_y')+
  # scale_y_log10(
  # labels=scales::label_number()
  # breaks = 10^seq(-3,2,by=1), lim=10^c(-3,2),
  # labels = expression(10^-3,10^-2,10^-1,1,10^1,10^2)
  # ) + 
  scale_x_continuous(breaks = sort(unique(d3$Ninc))[-2]) +
  labs(colour='scheme', y='time /s',
       title = 'Timing vs incident angles',
       subtitle = 'Nmax = 8, Npart = 4 and 8') +
  theme(legend.position = 'top', plot.margin = margin(20,5,20,5))

pp4p

