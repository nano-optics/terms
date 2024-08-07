##
## Functions for creating specific geometries, and utilities
##


##' @title Clusters of particles
##' @description Defines various cluster geometries and exports in a format suitable for TERMS
##' @describeIn clusters write cluster positions to input file
##' @param N number of particles
##' @param cl_fun cluster function
##' @param radius particle radius
##' @param label particle material label
##' @param out filename
##' @param digits accuracy
##' @param ... extra arguments passed to cl_fun
##' @return returns scatterers positions and sizes for an input file
##' @importFrom stats setNames
##' @importFrom utils read.table strcapture write.table 
##' @import tibble
##' @export
##' @examples 
##' cluster_positions()
cluster_positions <- function(N=5, cl_fun=cluster_chain, 
                              radius = 50, label="Au", ..., out='', digits = 8){
  c <- cl_fun(N=N, ...)

  write.table(cbind(label, round(cbind(t(c$positions), radius), digits=digits)), 
              file = out, 
              col.names = FALSE, append = TRUE, 
              row.names = FALSE, quote = FALSE)
  
}


##' @description Helix of particles
##' @describeIn clusters helical cluster
##' @param a semi-axis
##' @param b semi-axis
##' @param c semi-axis
##' @param R0 helix radius
##' @param pitch helix pitch
##' @param delta helix angle step
##' @param delta0 helix start angle
##' @param hand helix handedness
##' @export
##' @examples 
##' cluster_helix()
cluster_helix <- function (N = 5, a = 10, b = 10, c = 50, R0 = 100,
                           pitch = 200, 
                           delta = pi/5, delta0 = 0, hand = 1,  ...) 
{
  
  
  ph = seq(from = delta0, by = delta, length = N)
  x = R0 * cos(ph)
  y = R0 * sin(ph)
  z = hand * ph * pitch / (2*pi)
  positions <- rbind(x, y, z=z-mean(z))
  
  # angles calculation
  xp =  - y
  yp =   x
  zp =  hand * pitch / (2*pi)
  n = sqrt(xp^2+yp^2+zp^2)
  
  phi =  atan2(yp, xp);
  theta =  acos(zp/n);
  psi = theta*0; # don't care for now
  angles <- rbind(phi, theta, psi)
  
  sizes <- matrix(c(a=a, b=b, c=c), nrow=3, ncol=N)
  structure(list(positions = positions, 
                 sizes = sizes, 
                 angles = as.matrix(angles), 
                 R0 = R0), class = "cluster")
}


##' @description Chain of particles
##' @describeIn clusters linear chain cluster
##' @param a semi-axis
##' @param b semi-axis
##' @param c semi-axis
##' @param pitch chain pitch
##' @param rot rotation matrix applied to each particle
##' @export
##' @examples 
##' cluster_chain()
cluster_chain <- function (N=5, pitch = 500, 
                           a = 50, b = 30, c = b, 
                           rot=rotation_euler_passive(0,0,0)) 
{
  positions <- rbind(x = (seq_len(N) - (N+1)%/%2) * pitch, y = 0, z = 0)
  positions <- t(rot) %*% positions
                     
  sizes <- equal_sizes(a = a, b = b, c = c, N = N)
  angles <- equal_angles(0, 0, 0, N = N)
  structure(list(positions = positions, sizes = sizes, angles = angles), 
            class = "cluster")
}



##' @description Core-satellite cluster of spheres
##' @describeIn clusters core-satellite cluster
##' @param Rcore core radius
##' @param Rsat satellite radius
##' @param gap gap distance
##' @param exclusion minimum exclusion distance for hc positions
##' @export
##' @examples 
##' cluster_satellite()
cluster_satellite <- function(N = 30, Rcore = 30, Rsat = 4,gap = 0.1, 
                              exclusion = 10, ...){
  
  R <- Rcore + Rsat + gap
  
  positions <- R * sample_fibonacci(N)
  sizes <- Rsat + 0*positions
  angles <- 0*positions
  
  structure(list(positions = positions, 
                 sizes = sizes, 
                 angles = angles, 
                 R0 = Rcore), class = "cluster")
}



##' @description Matrix of equal sizes
##' @describeIn utility equal sizes
##' @param a semi-axis
##' @param b semi-axis
##' @param c semi-axis
##' @param N number of particles
##' @export
equal_sizes <- function (a, b, c, N) 
{
  rbind(a = rep(a, N), b = rep(b, N), c = rep(c, N))
}


##' @description Matrix of equal angles
##' @describeIn utility equal angles
##' @param phi Euler angle
##' @param theta Euler angle
##' @param gamma Euler angle
##' @param N number of particles
##' @export
equal_angles <- function (phi, theta, gamma, N) 
{
  rbind(phi = rep(phi, N), theta = rep(theta, N), gamma = rep(gamma, 
                                                              N))
}

##' Fibonacci coverage of a sphere
##'
##' Produces a set of points that covers rather uniformly the unit sphere with N points
##' with a spiral-like pattern based on a Fibonacci sequence
##' @title sample_fibonacci
##' @param N number of points
##' @export
##' @family low_level sample fibonacci sampling of a sphere
sample_fibonacci <- function(N=301){
  N0 <- N
  if(N%%2 == 1) N0 <- N+1
  P <- (N0-1)/2
  ii <- seq(-P,P,by=1)
  # note: uses latitude (internally), not colatitude
  # but we don't use angles, just xyz in the end
  lat <- asin(2*ii/N0)
  Phi <- (1+sqrt(5))/2
  long <- 2*pi*ii/Phi
  
  # return N xyz positions
  rbind(x = cos(long)*cos(lat), 
        y = sin(long)*cos(lat), 
        z = sin(lat))[,1:N]
}

