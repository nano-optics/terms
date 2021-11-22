
##' @noRd
#' @export
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


##' @noRd
#' @export
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

##' @noRd
#' @export
euler_to_axisangle <- function(a,b,c){
  
  M = rotation_euler_active(a,b,c)
  theta = acos((M[1, 1] + M[2, 2] + M[3, 3] - 1) / 2)+1e-5
  e1 = (M[3, 2] - M[2, 3]) / (2*sin(theta))
  e2 = (M[1, 3] - M[3, 1]) / (2*sin(theta))
  e3 = (M[2, 1] - M[1, 2]) / (2*sin(theta))
  
  c(e1, e2, e3, theta)
  
}


##' @noRd
#' @export
match_closest <- function(x, y, tol=1e-5){
  
  breaks <- sort(c(y-0.5*tol, y+0.5*tol))
  f <- cut(x, breaks=breaks, include.lowest = TRUE)
  rep(sort(y),each=2)[as.numeric(f)]
  
}

