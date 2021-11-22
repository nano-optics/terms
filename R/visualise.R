
##' @title Interactive display of cluster geometries
##' @description Displays a cluster in X3D format
##' @param cl cluster 
##' @param viewpoint optional viewpoint (3-vector)
##' @param width display width
##' @param height display height
##' @param scale scaling of axes
##' @param ... extra arguments passed to cluster_to_x3d
##' @return returns X3D object to embed in a html document with suitable X3D support
##' @import glue
##' @export
x3d_scene <- function(cl, viewpoint=NULL, width="300px", height="300px", scale=100, ...){
  
  if(is.null(viewpoint)) viewpoint <- c(0,0,1) * 1.1 * max(apply(cl$positions,1,max))
  
  header <- glue('<X3D width="{width}" height="{height}" class="x3d_scene"><scene>')
  footer <- glue('<Shape><IndexedLineSet colorPerVertex="false" colorIndex="0 1 2" coordIndex="0 1 -1 0 2 -1 0 3 -1">
<Coordinate point="0 0 0 {scale} 0 0 0 {scale} 0 0 0 {scale}"/><Color color="1 0 0 0 1 0 0.2 0.2 1"/></IndexedLineSet></Shape> 
<viewpoint centerOfRotation="0 0 0" position="{paste(viewpoint, collapse=" ")}" orientation="0 0 1 0" />
</scene>
</X3D>')
  
  cat(header, cluster_to_x3d(cl, ...), footer, sep="\n") 
  
}

#' @noRd
#' @export
particle_to_x3d <- function(position, size, angle, material = 'gold'){
  
  if(material == 'gold'){
    diffuseColor <- "1.0000000 0.8431373 0.0000000"
    specularColor<- "0.8549020 0.6470588 0.1254902"
  } else if(material == 'silver'){
    diffuseColor <- "0.7529412 0.7529412 0.7529412"
    specularColor<- "0.827451 0.827451 0.827451"
  } else {
    diffuseColor <- "0.3921569 0.5843137 0.9294118"
    specularColor<- "0.5294118 0.8078431 0.9803922"
  } 

  
  rotation <- euler_to_axisangle(angle[1],angle[2],angle[3])
  scale <- size
  translation <- position
  
  glue::glue('<transform scale="{paste(scale,collapse=" ")}" rotation="{paste(rotation,collapse=" ")}" translation="{paste(translation,collapse=" ")}">
      <shape>
        <appearance>
          <material ambientIntensity="0.4"  diffuseColor="{diffuseColor}" shininess="0.16" specularColor="{specularColor}"/>
        </appearance>
        <sphere radius="1" />
      </shape>
    </transform>
    ')
}


#' @noRd
#' @export
cluster_to_x3d <- function(cl, default_material = 'dielectric'){
  N <- ncol(cl$sizes)
  mean_pos <- apply(cl$positions, 1, mean)
  material <- rep(default_material, N)
  gold <- grepl('Au', cl$materials)
  silver <- grepl('Ag', cl$materials)
  if(any(gold)) material[gold] <- 'gold'
  if(any(silver)) material[silver] <- 'silver'
  paste(lapply(1:N, function(ip) particle_to_x3d(cl$positions[,ip]-mean_pos,cl$sizes[,ip],cl$angles[,ip], material[ip])), sep='\n', collapse='\n')
}




#' @noRd
#' @export
get_geometry <- function(input = 'input'){
  
  # this assumes all scatterers of the same kind
  # should use readLines instead, since different lengths
  # g <- read.table(pipe(paste("sed '1,/Scatterers/d'", input)))
  
  # first, find the Scatterers line
  rl <- readLines(input)
  scat <- which(grepl('^Scatterers', rl))
  # now re-read as table
  g <- read.table(input, skip = scat, header = FALSE)
  
  # cs <- grepl("@", g)
  if(ncol(g) < 9) # spherical scatterers, complete with dummy
    g <- cbind(g, 
               matrix(0, nrow=nrow(g), ncol=8-ncol(g)), 
               matrix(1, nrow=nrow(g), ncol=1))
  
  names(g) <- c("tag", 'x','y','z','r', 'alpha','beta','gamma','chi') 
  
  g$label <- paste0(g$tag, seq(1,nrow(g)))
  
  g
}

#' @noRd
#' @export
cluster_geometry <- function(ge){
  positions <- t(ge[,c("x","y","z")])
  sizes <- t(ge[,c("r","r","r")] * cbind(1/ge[,c("chi")],1/ge[,c("chi")],1))
  angles <- t(ge[,c("alpha","beta","gamma")])
  materials <- ge$tag
  structure(list(positions=positions,angles=angles,sizes=sizes, materials=materials),class='cluster')
}

#' @noRd
#' @export
visualise_rgl <- function(cl, outfile=NULL, show_core=FALSE, ...){
  rgl.ellipsoids(cl$positions, cl$sizes, cl$angles, ...)
  
  if("R0" %in% names(cl) && show_core) 
    rgl::rgl.spheres(0,0,0, radius=cl$R0, col="grey")
  
  if(!is.null(outfile))
    rgl::rgl.snapshot( outfile, fmt = "png", top = TRUE )
}


#' @noRd
#' @export
rgl.ellipsoid <- function (x=0, y=0, z=0, a = 1, b=1, c=3, phi=0, theta=0, psi=0,
                           subdivide = 2, ...)
{
  # sigma <- diag(c(1,1,3))
  # o <- ellipse3d(sigma, centre = c(0,0,0), subdivide=2)
  sphere <- rgl::subdivision3d(rgl::icosahedron3d(...), subdivide)
  class(sphere) <- c("mesh3d","shape3d")
  
  norm <- sqrt(sphere$vb[1, ]^2 + 
                 sphere$vb[2, ]^2 + 
                 sphere$vb[3, ]^2 )
  for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ]/norm
  sphere$vb[4, ] <- 1
  sphere$normals <- sphere$vb
  result <- rgl::scale3d(sphere, a,b,c)
  
  # shade3d(result, col='red')
  
  rotM <- rotation_euler_passive(phi,theta,psi)
  result <- rgl::rotate3d(result,matrix=rotM)
  
  # shade3d(result, col='gold')
  
  result <- rgl::translate3d(result, x,y,z)
  invisible(result)
}


#' @noRd
#' @export
rgl.ellipsoids <- function(positions, sizes, angles, colour = "red", ...){
  
  N <- NCOL(positions)
  colours <- rep(colour, length.out=N)
  ll <- lapply(seq(1,N), function(ii)
    rgl.ellipsoid(positions[1,ii],positions[2,ii],positions[3,ii],
                  sizes[1,ii],sizes[2,ii],sizes[3,ii],
                  angles[1,ii],angles[2,ii],angles[3,ii], col = colours[ii], ...))
  
  rgl::shapelist3d(ll)
  
}

#' @noRd
#' @export
rgl_annotate <- function(){
  rgl::axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
  rgl::axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
  rgl::axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
  rgl::axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
  rgl::title3d('','','x axis','y axis','z axis')
}

##' @noRd
#' @export
generate_scad <- function(cl, out='', res=10,
                          x=cl$positions[1,],
                          y=cl$positions[2,],
                          z=cl$positions[3,],
                          a=cl$sizes[1,],
                          b=cl$sizes[2,],
                          c=cl$sizes[3,],
                          r=1+0*cl$sizes[3,],
                          phi=cl$angles[1,]*180/pi,
                          theta=cl$angles[2,]*180/pi,
                          psi=cl$angles[3,]*180/pi){
  
  npart <- length(x)
  
  tpl <- "
  x = [{{x}}];
  y = [{{y}}];
  z = [{{z}}];
  a = [{{a}}];
  b = [{{b}}];
  c = [{{c}}];
  r = [{{r}}];
  phi = [{{phi}}];
  theta = [{{theta}}];
  psi = [{{psi}}];
  
  for (i=[0:{{npart}}])
   translate([x[i], y[i], z[i]]) multmatrix([ [cos(phi[i])*cos(theta[i])*cos(psi[i]) - sin(phi[i])*sin(psi[i]), 
              -cos(phi[i])*cos(theta[i])*sin(psi[i]) - sin(phi[i])*cos(psi[i]),
              cos(phi[i])*sin(theta[i]), 0],
             [sin(phi[i])*cos(theta[i])*cos(psi[i]) + cos(phi[i])*sin(psi[i]),
              -sin(phi[i])*cos(theta[i])*sin(psi[i]) + cos(phi[i])*cos(psi[i]),
              sin(phi[i])*sin(theta[i]), 0],
             [-sin(theta[i])*cos(psi[i]), 
              sin(theta[i])*sin(psi[i]), 
              cos(theta[i]), 0],
             [0,0,0,1]
] )  scale([a[i], b[i], c[i]]) sphere(r = r[i], $fn={{res}});
  "
  
  st <- glue::glue(tpl, .open = "{{", .close = "}}", 
                   x=paste0(x,collapse=","), 
                   y=paste0(y,collapse=","),
                   z=paste0(z,collapse=","),
                   a=paste0(a,collapse=","), 
                   b=paste0(b,collapse=","),
                   c=paste0(c,collapse=","), 
                   r=paste0(r,collapse=","), 
                   phi=paste0(phi,collapse=","),
                   theta=paste0(theta,collapse=","),
                   psi=paste0(psi,collapse=","),
                   npart = npart-1, res=res)
  
  cat(st, file=out)
  invisible(st)
}





#' get_geometry <- function(input = 'input'){
#'   
#'   g <- read.table(pipe(paste("sed '1,/Scatterers/d'", input)))
#'   
#'   if(ncol(g) < 9) # basic scatterers, complete with dummy
#'     g <- cbind(g, matrix(NA, nrow=nrow(g), ncol=9 - ncol(g)))
#'   
#'   names(g) <- c("tag", 'x','y','z','r', 'alpha','beta','gamma','chi') 
#'   
#'   g$label <- paste0(g$tag, seq(1,nrow(g)))
#'   
#'   g
#' }
#' 
#' 
#' cluster_geometry <- function(ge){
#'   positions <- t(ge[,c("x","y","z")])
#'   sizes <- t(ge[,c("r","r","r")] * cbind(1/ge[,c("chi")],1/ge[,c("chi")],1))
#'   angles <- t(ge[,c("alpha","beta","gamma")])
#'   
#'   structure(list(positions=positions,angles=angles,sizes=sizes),class='cluster')
#' }


# 
# input <- readLines('input')
# rg <- paste(input[(grep("Scatterer",input)+1):(length(input)-1)], collapse="\n")
# 
# add_geometry <- function(descriptor){
#   geometry <- setNames(read.table(textConnection(descriptor), 
#                                   colClasses = c("character","double","double","double","double")), 
#                        c("tag", 'x','y','z','r'))
#   geometry$label <- paste0(geometry$tag, seq(1,nrow(geometry)))
#   
#   geometry
# }
# 
# ge <- add_geometry(rg)
# 
