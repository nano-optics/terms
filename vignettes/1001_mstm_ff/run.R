library(glue)

tpl <- "number_spheres
{N}
print_sphere_data
.true.
sphere_data
{scatterers}
end_of_sphere_data
length_scale_factor
{k0}
random_orientation
.true.
mie_epsilon
-{Lmax}
output_file
{results}
end_of_options
"

xyz <- matrix(1:9,3,3, byrow=T)
epsilon <- dielectric::epsAu(633)

spheres <- function(positions, radii, epsilon){
  # sphere_data
  # 0.d0,0.d0,0.d0,1.d0,(1.5d0,0.d0)
  # 0.d0,0.d0,2.d0,1.d0,(1.5d0,0.d0)
  # end_of_sphere_data
  
  m <- sqrt(epsilon)
  # glue('{positions[,1]}d0,{positions[,2]}d0,{positions[,3]}d0,{radii}d0,({Re(m)}d0,{Im(m)}d0)')
  glue('{positions[,1]},{positions[,2]},{positions[,3]},{radii},({Re(m)},{Im(m)})')

}

# cl <- cda::cluster_helix(N = 8, a = 50, b=50, c=50, R0 = 150, pitch = 400, delta = pi/4)

cl <- cda::cluster_helix(N = 3, a = 50, b=50, c=50, R0 = 120, pitch = 400, delta = pi/4)
# cda::visualise(cl, show_core = FALSE)

input_file <- function(cluster = cl, 
                       wavelength = 633, Lmax=12, out = glue('{wavelength}.inp')){
  
  k0 <- 2*pi/wavelength
  epsilon <- dielectric::epsAu(wavelength)$epsilon
  scatterers <- glue_collapse(spheres(t(cluster$positions), cluster$sizes[1,], epsilon), sep = '\n')
  N <- ncol(cluster$positions)
  results <- gsub('.inp','.dat', out)
  cat(glue(tpl), file = out)
}




for(lambda in seq(400, 800, by=50)){
  input_file(wavelength = lambda)
  system(glue('./mstm {lambda}.inp'))
}


