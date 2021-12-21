# terms 1.0.0

## Website

- new vignette for alpha tensor

## Bug fix

- LR vs RL error in calculation of scattering circular dichroism

## User-facing changes

- MapOaQuantity [p] new format

## New functionality

- implemented conversion of multipoles $l<=3$ into cartesian "alpha tensor"

# terms 0.9.9


## Website

- updated examples for release
- new version of pkgdown, some css tweaks

## R code

- some helper functions for hdf5 format
- vignette schematics in x3d format

## Fortran

- added support for hdf5
- added computation of local degree of optical chirality
- added computation of Stokes vectors and Mueller matrices
- multi-wavelength calculations of near-fields
- near field maps have more options (E, B, C) and include information about material regions


## compilation

- Cmake is now used for building terms


# terms 0.9.0

## Website

- initial website using pkgdown

## R code

- basic package structure
- examples in Rmd format ("vignettes")
- basic utilities in `R/`
