setwd(here::here('vignettes/1003_cda'))

## ----load----
library(terms)
library(ggplot2)
library(R.matlab)
theme_set(theme_grey())

## ----run----
system("../../build/terms input > log")
system("../../build/terms input_dip > log")


## ----read----
## grab CDA results

dmat <- readMat('fingers.mat')
cda <- data.frame(wavelength = rep(as.vector(dmat$wavelength), 3),
                  average = 2*c(dmat$xsec[,1],dmat$xsec[,2],dmat$xsec[,3]),
                  dichroism = -2*c(dmat$xsec[,5],dmat$xsec[,6],dmat$xsec[,7]),
                  crosstype = rep(c('Ext','Abs','Sca'), each = length(dmat$wavelength)))

# grab the TERMS results
xs <- terms::consolidate_xsec('cross_sections.h5')
xs_dip <- terms::consolidate_xsec('cross_sections_dip.h5')

# plot OA data
## ----oa----

# glimpse(xs$mCOA)
mCOAt <- subset(xs$mCOA, variable == 'total')
mCOAt_dip <- subset(xs_dip$mCOA, variable == 'total')

library(egg)
library(gggrid)
# library(png)
# png <- png::readPNG('fingers_dimer.png')

p1 <- ggplot(mCOAt, aes(wavelength, average,colour=crosstype)) +
  facet_wrap(~crosstype)+
  # grid_panel(data=tibble(crosstype='Abs',wavelength=400,average=1e4),grid::rasterGrob(png, width = unit(3,'cm'), x=0,hjust=0)) +
  geom_area(data=mCOAt, aes(fill=crosstype),colour=NA) +
  geom_line(data=mCOAt_dip, aes(lty='L[max]==1'),lwd=0.5) +
  geom_line(data=cda, aes(lty='CDA'),lwd=0.8) +
  guides(colour='none')+
  scale_y_continuous(lim=c(0,7700),expand=c(0,0)) +
  scale_x_continuous(lim=c(460,740),expand=c(0,0)) +
  scale_linetype_manual(values=c(2,1), labels = scales::parse_format()) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Pastel1') +
  guides(linetype=guide_legend(label.hjust = 0))+
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(type), fill='exact', linetype='',
       subtitle = "Orientation-averaged cross-sections")

# p1

p2 <- ggplot(mCOAt, aes(wavelength, dichroism,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_area(data=mCOAt, aes(fill=crosstype),colour=NA) +
  geom_line(data=mCOAt_dip, aes(lty='L[max]==1'),lwd=0.5) +
  geom_line(data=cda, aes(lty='CDA'),lwd=0.8) +
  guides(colour='none')+
  scale_x_continuous(lim=c(460,740),expand=c(0,0)) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Pastel1') +
  guides(linetype=guide_legend(label.hjust = 0))+
  scale_y_continuous(limits = c(-740,260), expand=c(0,0)) + #limits = symmetric_range
  scale_linetype_manual(values=c(2,1), labels = scales::parse_format()) +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(type), fill='exact', linetype='',
       subtitle = "Orientation-averaged circular dichroism")

g <- egg::ggarrange(p1, p2,ncol=1, draw = F)

grid.draw(g)
# ggsave('comparison-cda.pdf', g, width=8, height=6)




## ----cleanup----
# unlink("*.dat")

