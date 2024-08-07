% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualise.R
\name{x3d_scene}
\alias{x3d_scene}
\title{Interactive display of cluster geometries}
\usage{
x3d_scene(
  cl,
  viewpoint = c(0, 0, 100),
  orientation = NULL,
  width = "300px",
  height = "300px",
  scale = 100,
  ...
)
}
\arguments{
\item{cl}{cluster}

\item{viewpoint}{viewpoint position (3-vector)}

\item{orientation}{optional viewpoint orientation (axis-angle vector: x,y,z,angle)}

\item{width}{display width}

\item{height}{display height}

\item{scale}{size of axes}

\item{...}{extra arguments passed to cluster_to_x3d}
}
\value{
returns X3D object to embed in a html document with suitable X3D support
}
\description{
Displays a cluster in X3D format
}
