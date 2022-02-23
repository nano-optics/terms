
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TERMS <img src="man/figures/logo.png" width="120" align="right" />

TERMS stands for **T**-matrix for **E**lectromagnetic **R**adiation with
**M**ultiple **S**catterers — it is a Fortran program to simulate the
near-field and far-field optical properties of collections of particles.
TERMS solves rigorously the Maxwell equations via the superposition
*T*-matrix method, where incident and scattered fields are decomposed
into a basis of multipolar electric and magnetic spherical waves.

In a multiple-scattering problem the net field exciting a given particle
is composed of the incident field plus the scattering contribution from
neighbouring particles, restulting in a coupled system of equations to
be solved for the total fields. TERMS implements several algorithms to
describe the self-consistent electromagnetic interaction between
multiple scatterers, and from there compute optical properties such as
absorption, scattering, extinction, circular dichroism, as well as
near-field intensities and the local degree of optical chirality.

By describing the incident and scattered fields in a basis of spherical
waves the *T*-matrix framework lends itself to analytical formulas for
orientation-averaged quantities such as far-field cross-sections and
near-field intensities, greatly reducing the computational time needed
to simulate particles and systems of particles in random orientation.

### Features

The possible computations are divided into three main modes:

-   Far-field quantities (absorption, scattering, extinction, circular
    dichroism) for multiple wavelengths and angles of incidence, as well
    as orientation-averages
-   Near-field calculations for multiple wavelengths and incident
    angles, also computing the local degree of chirality, as well as
    orientation-averages
-   Stokes parameters and differential scattering cross-sections for
    multiple incidence or scattering angles

The computational cost scales with the size of the linear system,
proportional to the number of particles Np, and to the square of the
maximum multipolar order Nmax. On a typical PC we may treat up to \~500
particles with Nmax=1, and a dimer with Nmax up to \~60.

Notable features of TERMS include:

-   Incident plane waves along arbitrary directions, with linear or
    circular polarisation
-   Built-in calculation of individual *T*-matrices for coated spheres;
    import of general *T*-matrices from other programs
    (e.g. [SMARTIES](https://www.victoria.ac.nz/scps/research/research-groups/raman-lab/numerical-tools/smarties))
-   Built-in dielectric functions for common materials such as Au, Ag,
    Al, Cr, Pt, Pd, Si, and Water, or from tabulated values
-   Per-layer absorption in layered spheres
-   Orientation-averaging of far-field cross-sections, as well as linear
    and circular dichroism
-   Near-field maps of electric and magnetic field components, E^2, E^4,
    local degree of chirality
-   Calculation of the global cluster *T*-matrix
-   “Masking” of specific multipolar orders
-   Calculation of Stokes parameters, phase matrix, differential
    scattering
-   Plain text or HDF5 I/O format
-   Possible compilation in quad-precision

### System requirements

-   Fortran 90 compiler
-   Cmake
-   (optional) HDF5 library
-   (optional) LAPACK

The electromagnetic field is expanded in the basis of vector spherical
waves, with the Bessel/Hankel functions computed using
[TOMS644](http://www.netlib.org/toms-2014-06-10/644) library (source
included in **TERMS**). Determining the collective scattering amounts to
either solving or inverting a large linear system, which is done using
[LAPACK](http://www.netlib.org/lapack/). All the relevant LAPACK
routines are included in **TERMS**, but it is recommended to link with
your system installed BLAS/LAPACK at compilation stage, because it can
enable multi-threading during runtime.  
Results are output in plain text files, or, alternatively, in `HDF5`
data format, which requires suitable `hdf5` libraries to be available on
your system.

### Compilation

We recommend using the cross platform compilation tool
[cmake](https://cmake.org/), to resolve all dependencies for your
system. From within the `build/` directory, type

    > cmake ..
    > make

to produce the executable `terms`. Note: if Cmake doesn’t find hdf5 (or
another library path) by default, you may need to export it explicitly
beforehand. For example on MacOS, the following has proved useful:

    # brew install hdf5
    export HDF5_ROOT=/usr/local/Cellar/hdf5/1.12.0_3/

Alternatively, a minimal build script is available in the `build/`
directory; the program can be compiled by executing `bash buildTERMS.sh`
from a Linux terminal with [bash](https://www.gnu.org/software/bash/).
Edit ‘buildTERMS.sh’ to specify a compiler other than
[GFortran](https://gcc.gnu.org/wiki/GFortran).

#### Downloading the code

We recommend downloading the [latest release
here](https://github.com/nano-optics/terms/releases)
\[`terms_code_1.0.0(.zip|.tar)`\]. You can also browse/clone/fork the
[entire repository](https://github.com/nano-optics/terms), but note that
it contains many files used to generate the website, which are not
relevant for using TERMS.

#### Repository structure

For developer’s convenience, we’ve structured this repository as a R
package, to benefit from R’s comprehensive documentation ecosystem and
define some helper functions for data post-processing/plotting in R.
This means that several sub-folders are not relevant for TERMS users,
especially if you don’t use R for data analysis/plotting. The directory
structure is as follows:

`build/`  
contains minimal compilation scripts, as an alternative to Cmake

`src/`  
contains the TERMS Fortran code, and external libraries

The following directories are not needed for users of the program.

`docs/`  
this directory is auto-generated by a script, and serves the [published
website accessible at
nano-optics.ac.nz/terms](http://nano-optics.ac.nz/terms)

`inst`  
extra files needed for generating the website

`man/`  
documentation for R convenience functions

`pkgdown/`  
extra files needed for generating the website

`R/`  
R helper functions, which can be loaded with `library(terms)` after
installing this package locally

`test/`  
Example `input` files to test that the program works correctly. Also
includes a few incidence files for Lebedev cubature.

`vignettes/`  
a collection of standalone examples illustrating of the different
features of the program, as well as timing and convergence tests. These
examples are run automatically and the output is [published on the
website](http://nano-optics.ac.nz/terms/).

### The input file

When running the stand-alone executable, main input parameters are read
from a plain text input file (line by line and from left to right). Each
line is interpreted as a sentence and split into space-separated words.
The first (left-most) word is interpreted as a keyword, and the
subsequent words as arguments for that keyword. In each sentence, text
from the first word starting with the hash character (`#`) is
interpreted as a human-readable comment and thus ignored by the program.
All the supported keywords and corresponding arguments are [documented
on this page](http://nano-optics.ac.nz/terms/articles/Keywords.html).
The order of keywords generally doesn’t matter, with just two
exceptions: `ModeAndScheme` must be the first keyword, and `Scatterers`
must be the last. Two examples of input files are provided in the
`/test` directory.

## Citing TERMS

If you use TERMS, please cite the following user guide, as well as other
publications listed below if relevant:

<div style="display: none;">

<sup>1</sup>,,<sup>2</sup>,<sup>3</sup>,<sup>4</sup>,<sup>5</sup>,<sup>6</sup>,<sup>7</sup><sup>8</sup>

</div>

<div id="refs" class="references csl-bib-body">

<div id="ref-Schebarchov:2022wc" class="csl-entry">

<span class="csl-left-margin">(1) </span><span
class="csl-right-inline">Schebarchov, D.; Fazel-Najafabadi, A.; Le Ru,
E. C.; Auguié, B. Multiple Scattering of Light in Nanoparticle
Assemblies: User Guide for the Terms Program. *Journal of Quantitative
Spectroscopy and Radiative Transfer* **2022**, 108131.
https://doi.org/<https://doi.org/10.1016/j.jqsrt.2022.108131>.</span>

</div>

<div id="ref-Schebarchov:2021ut" class="csl-entry">

<span class="csl-left-margin">(2) </span><span
class="csl-right-inline">Schebarchov, D.; Fazel-Najafabadi, A.; Le Ru,
E. C.; Auguié, B. *TERMS Website*; 2021.
<https://doi.org/10.5281/zenodo.5703291>.</span>

</div>

<div id="ref-Somerville:2016aa" class="csl-entry">

<span class="csl-left-margin">(3) </span><span
class="csl-right-inline">Somerville, W. R. C.; Auguié, B.; Le Ru, E. C.
SMARTIES: User-Friendly Codes for Fast and Accurate Calculations of
Light Scattering by Spheroids. *J. Quant. Spectrosc. Ra.* **2016**,
*174*, 39–55. <https://doi.org/10.1016/j.jqsrt.2016.01.005>.</span>

</div>

<div id="ref-schebarchov2019mind" class="csl-entry">

<span class="csl-left-margin">(4) </span><span
class="csl-right-inline">Schebarchov, D.; Le Ru, E. C.; Grand, J.;
Auguié, B. Mind the Gap: Testing the Rayleigh Hypothesis in *T*-Matrix
Calculations with Adjacent Spheroids. *Optics express* **2019**, *27*
(24), 35750–35760. <https://doi.org/10.1364/OE.27.035750>.</span>

</div>

<div id="ref-Lee:2020aa" class="csl-entry">

<span class="csl-left-margin">(5) </span><span
class="csl-right-inline">Lee, S.; Hwang, H.; Lee, W.; Schebarchov, D.;
Wy, Y.; Grand, J.; Augui’e, B.; Wi, D. H.; Cort’es, E.; Han, S. W.
Core–Shell Bimetallic Nanoparticle Trimers for Efficient
Light-to-Chemical Energy Conversion. *ACS Energy Letters* **2020**, *5*
(12), 3881–3890. <https://doi.org/10.1021/acsenergylett.0c02110>.</span>

</div>

<div id="ref-Fazel-Najafabadi:2021uq" class="csl-entry">

<span class="csl-left-margin">(6) </span><span
class="csl-right-inline">Fazel-Najafabadi, A.; Schuster, S.; Auguié, B.
Orientation Averaging of Optical Chirality Near Nanoparticles and
Aggregates. *Physical Review B* **2021**, *103* (11), 115405.
<https://doi.org/10.1103/PhysRevB.103.115405>.</span>

</div>

<div id="ref-Fazel-Najafabadi:2022ud" class="csl-entry">

<span class="csl-left-margin">(7) </span><span
class="csl-right-inline">Fazel-Najafabadi, A.; Auguié, B. Orientation
Dependence of Optical Activity in Light Scattering by Nanoparticle
Clusters. *Mater. Adv.* **2022**, –.
<https://doi.org/10.1039/D1MA00869B>.</span>

</div>

<div id="ref-Fazel-Najafabadi:2022aa" class="csl-entry">

<span class="csl-left-margin">(8) </span><span
class="csl-right-inline">Fazel-Najafabadi, A.; Auguié, B.
Orientation-Averaged Light Scattering by Nanoparticle Clusters:
Far-Field and Near-Field Benchmarks of Numerical Cubature Methods,
2022.</span>

</div>

</div>
