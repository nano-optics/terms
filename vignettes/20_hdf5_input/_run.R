setwd(here::here('vignettes/20_hdf5_input/'))
## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())
library(patchwork)


## ----run----
system("../../build/terms input_spheroid1_dimer > log1")

system("../../build/terms input_spheroid2_dimer > log2")

## ----read----
xs1 <- terms::consolidate_xsec('JCMsuite_dimer.h5')
xs2 <- terms::consolidate_xsec('smarties_dimer.h5')

## ----plot----
mCOA <- rbind(cbind(xs1$mCOA, tmatrix = 'SMARTIES'),
              cbind(xs2$mCOA, tmatrix = 'JCMsuite')) 
mCOA$crosstype <- factor(mCOA$crosstype, levels = c("Ext", "Abs", "Sca"),
                         labels = c("Extinction", "Absorption", "Scattering")
)
mCOAt <- subset(mCOA, variable == 'total')

p1 <- ggplot(mCOA, aes(wavelength, average, alpha=tmatrix,
                       linetype = tmatrix,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAt |> filter(tmatrix=='SMARTIES',wavelength>420, 
                                 wavelength <760),
            lwd=0.8,lty=1) +
  geom_point(data=mCOAt |> filter(tmatrix=='JCMsuite')) +
  scale_alpha_manual(values=c(0.8,1)) +
  guides(colour='none') +
  scale_x_continuous(expand=c(0,0))+ 
  scale_colour_brewer(palette='Set1') +
  theme(legend.position	="inside", legend.position.inside = c(0.915,0.7))+
  labs(x = expression("wavelength /nm"), 
       y = expression("avg. cross-sec. /"*nm^2),
       colour = expression(N), linetype="", pch="",alpha="T-matrix") 

p2 <- ggplot(mCOA, aes(wavelength, dichroism, alpha=tmatrix,
                       linetype = tmatrix,colour=crosstype)) +
  facet_wrap(~crosstype)+
  annotate("segment", x = -Inf, xend=Inf, y=0,yend=0, lty=3, col='grey')+
  geom_line(data=mCOAt |> filter(tmatrix=='SMARTIES',wavelength>420, wavelength <760),
            lwd=0.8,lty=1) +
  geom_point(data=mCOAt |> filter(tmatrix=='JCMsuite')) +
  scale_alpha_manual(values=c(0.8,1)) +
  scale_x_continuous(expand=c(0,0))+ 
  guides(colour='none') +
  scale_colour_brewer(palette='Set1') +
  labs(x = expression("wavelength /nm"), 
       y = expression("dich. cross-sec. /"*nm^2),
       colour = expression(N), linetype="", pch="") +
  scale_y_continuous(limits = symmetric_range) +
  theme(legend.position = "none")

print(p1/p2)

