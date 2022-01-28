# terms 1.0.0

## Website

New vignettes:

- Polarised orientation-averaged near fields
- Reproducing cartesian `Higher-Order Polarizability Tensors' for a dimer
- Validation against MSTM v4 for orientation-averaged cross-sections
- Validation against coupled-dipole approximation for far-field circular dichroism

## Bug fix

- Block matrices LR and RL were swapped, resulting in a (small) error when computing scattering circular dichroism

## User-facing changes

- `MapOaQuantity [p]` has a new format (polarised or unpolarised); if specific polarisation(s) are requested, they are taken from the `Incidence` keyword

## New functionality

- implemented conversion of multipoles $l<=3$ into cartesian "alpha tensor"
- implemented polarised (L or R) orientation-averaged near-fields E2, B2, and LDOC (note: the results are incorrect inside spheres)

## Low level changes

- `calcOaExtField` renamed `calcOaNFUnpol`
- `calcOaLDOC` renamed `calcOaNF`; now combines calculation of LDOC and field intensities (many terms are common)
- new routine `alphaTensor`

## R utilities

- `export_cubature()` utility function to produce an incidence file suitable for TERMS from the `cubs` package.
- removed dependency on `reshape2`; using `tidyr` instead

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
