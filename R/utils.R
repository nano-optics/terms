

##' @title Utility functions
##' @description Generate an incidence file for spherical cubature
##' @describeIn utility export a spherical cubature
##' @param q data.frame with angles and weights, from cubs::cubs()
##' @param out filename
##' @importFrom cubs cubs
##' @export
export_cubature <- function(q = cubs::cubs(N = 10, cubature = 'lebedev'), out=''){
  
  cat(nrow(q),'\n',file=out)
  write.table(format(cbind(q[,1:2], 0, q[,3]), digits = 15), file = out, quote=FALSE,
              row.names = FALSE, col.names = FALSE, append = TRUE)
  
}


##' @title Utility functions
##' @description Generate a dielectric function in suitable format for TERMS
##' @describeIn utility export a dielectric function
##' @param m data.frame with wavelength and epsilon, e.g. from dielectric::epsAu()
##' @param out filename
##' @importFrom dielectric epsAu
##' @export
export_dielectric <- function(m = dielectric::epsAu(seq(400,800)), out=''){
  
  wavelength <- m[["wavelength"]]
  epsilon <- m[["epsilon"]]
  write.table(format(cbind(wavelength, Re(epsilon), Im(epsilon)), digits = 15), file = out, quote=FALSE,
              row.names = FALSE, col.names = FALSE, append = FALSE)
  
}

##' @description Euler rotation matrix
##' @describeIn utility passive rotation matrix
##' @param phi Euler angle
##' @param theta Euler angle
##' @param psi Euler angle
##' @export
rotation_euler_passive <- function(phi, theta, psi){
  
  cosphi = cos(phi); cospsi = cos(psi); costheta = cos(theta)
  sinphi = sin(phi); sinpsi = sin(psi); sintheta = sin(theta)
  
  Rot = matrix(NA,3,3)
  Rot[1,1] = cosphi*costheta*cospsi - sinphi*sinpsi
  Rot[1,2] = sinphi*costheta*cospsi + cosphi*sinpsi
  Rot[1,3] = -sintheta*cospsi
  
  Rot[2,1] = -cosphi*costheta*sinpsi - sinphi*cospsi
  Rot[2,2] = -sinphi*costheta*sinpsi + cosphi*cospsi
  Rot[2,3] = sintheta*sinpsi
  
  Rot[3,1] = cosphi*sintheta
  Rot[3,2] = sinphi*sintheta
  Rot[3,3] = costheta
  
  Rot
}




##' @description Euler rotation matrix
##' @describeIn utility active rotation matrix
##' @param phi Euler angle
##' @param theta Euler angle
##' @param psi Euler angle
##' @export
rotation_euler_active <- function(phi, theta, psi){
  cosphi = cos(phi); cospsi = cos(psi); costheta = cos(theta);
  sinphi = sin(phi); sinpsi = sin(psi); sintheta = sin(theta);
  
  R=matrix(0, 3,3);
  
  R[1,1] = cosphi*costheta*cospsi - sinphi*sinpsi;
  R[2,1] = sinphi*costheta*cospsi + cosphi*sinpsi;
  R[3,1] = -sintheta*cospsi;
  
  R[1,2] = -cosphi*costheta*sinpsi - sinphi*cospsi;
  R[2,2] = -sinphi*costheta*sinpsi + cosphi*cospsi;
  R[3,2] = sintheta*sinpsi;
  
  R[1,3] = cosphi*sintheta;
  R[2,3] = sinphi*sintheta;
  R[3,3] = costheta;
  R
}



##' @description Axis-angle rotation from Euler angles
##' @describeIn utility axis-angle rotation
##' @param a Euler angle
##' @param b Euler angle
##' @param c Euler angle
##' @export
euler_to_axisangle <- function(a,b,c){
  
  M = rotation_euler_active(a,b,c)
  theta = acos((M[1, 1] + M[2, 2] + M[3, 3] - 1) / 2)+1e-5
  e1 = (M[3, 2] - M[2, 3]) / (2*sin(theta))
  e2 = (M[1, 3] - M[3, 1]) / (2*sin(theta))
  e3 = (M[2, 1] - M[1, 2]) / (2*sin(theta))
  
  c(e1, e2, e3, theta)
  
}


##' @description Extend a range symmetrically about 0
##' @describeIn utility symmetric range
##' @param range range (2-vector)
##' @export
symmetric_range <- function(range) 
{
  max_abs <- max(abs(range))
  c(-max_abs, max_abs)
}

##' @noRd
match_closest <- function(x, y, tol=1e-5){
  
  breaks <- sort(c(y-0.5*tol, y+0.5*tol))
  f <- cut(x, breaks=breaks, include.lowest = TRUE)
  rep(sort(y),each=2)[as.numeric(f)]
  
}

