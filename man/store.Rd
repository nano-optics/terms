% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/store.R
\name{store_xsec}
\alias{store_xsec}
\alias{consolidate_xsec}
\alias{consolidate_partials}
\title{Reshape cross-section results into a convenient format for post-processing and plotting}
\usage{
store_xsec(..., out = "xsec.rds")

consolidate_xsec(hdf5, verbose = TRUE, ...)

consolidate_partials(hdf5, verbose = TRUE)
}
\arguments{
\item{...}{extra arguments passed to the final list}

\item{out}{output Rds filename}

\item{hdf5}{filename}

\item{verbose}{logical: print attributes}
}
\value{
returns a list containing data.frames in long format
}
\description{
Read and store plain text cross-sections

Extracts commonly-used information from a HDF5 file storing far-field cross-sections (Mode=2), and reshapes the data into long-format data.frames suitable for plotting

Extracts partial absorption cross-sections from a HDF5 file storing far-field cross-sections (Mode=2), and reshapes the data into long-format data.frames suitable for plotting
}
\section{Functions}{
\itemize{
\item \code{store_xsec()}: store plain text cross-sections

\item \code{consolidate_xsec()}: consolidate cross-sections

\item \code{consolidate_partials()}: consolidate partial absorption cross-sections for multilayered spheres

}}
