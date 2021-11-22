
## common functions

#' @noRd
#' @export
lfOA <- c("cdAbsOA.dat","csAbsOA.dat",
          "cdExtOA.dat", "csExtOA.dat",
          "cdScaOA.dat", "csScaOA.dat")

#' @noRd
#' @export
lfLinear <- c("csAbs1X.dat", "csAbs2Y.dat",  
              "csExt1X.dat", "csExt2Y.dat",  
              "csSca1X.dat", "csSca2Y.dat")

#' @noRd
#' @export
lfCircular <- c("csAbs3R.dat", "csAbs4L.dat",
                "csExt3R.dat", "csExt4L.dat", 
                "csSca3R.dat", "csSca4L.dat")


## Read and store the OA data

#' @noRd
#'@export
store_xsec <- function(..., out = 'xsec.rds'){
  
  results <- list()
  
  if(all(file.exists(lfOA))){
    
    lOA <- lapply(lfOA, read_file)
    lCD <- lapply(lOA[c(1,3,5)], reshape2::melt, 
                  id=c("wavelength","crosstype"))
    lCS <- lapply(lOA[c(2,4,6)], reshape2::melt, 
                  id=c("wavelength","crosstype"))
    mCOA <- purrr::map2_df(lCD, lCS, process_OA)
    
    results <- c(results, list(mCOA=mCOA))
    
  } 
  
  ## Read in fixed orientation data
  
  if(all(file.exists(lfLinear))){
    
    lLFO <- lapply(lfLinear, read_file)
    lX <- lapply(lLFO[c(1,3,5)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    lY <- lapply(lLFO[c(2,4,6)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    mLFO <- purrr::map2_df(lX, lY, process_fixed)
    
    results <- c(results, list(mLFO = mLFO))
    
  }
  
  if(all(file.exists(lfCircular))){
    
    lCFO <- lapply(lfCircular, read_file)
    lL <- lapply(lCFO[c(1,3,5)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    lR <- lapply(lCFO[c(2,4,6)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    mCFO <- purrr::map2_df(lL, lR, process_fixed)
    
    results <- c(results, list(mCFO=mCFO))
  }
  
  # add extra metadata if needed
  results <- c(results, ...)
  
  saveRDS(results, file=out)
  invisible(results)
}


##' @title Reshape cross-section results into a convenient format for post-processing and plotting
##' @description Extracts commonly-used information from a HDF5 file storing far-field cross-sections (Mode=2), and reshapes the data into long-format data.frames suitable for plotting
##' @describeIn store consolidate cross-sections
##' @param hdf5 filename
##' @param verbose logical: print attributes
##' @param ... extra arguments passed to the final list
##' @return returns a list containing data.frames in long format
##' @import tibble
##' @importFrom rhdf5 h5read h5readAttributes
##' @export
consolidate_xsec <- function(hdf5, verbose = TRUE, ...){
  
  ld <- rhdf5::h5read(hdf5, "Far-Field")
  
  results <- list()
  
  if(("oa_incidence" %in% names(ld)) && !is.null(ld$oa_incidence)){ # Scheme=0: no OA
    
    att <- rhdf5::h5readAttributes(hdf5, "Far-Field/oa_incidence")
    if(!is.null(att) && verbose) message(att)
    
    nms <- c("cdAbsOA", "cdExtOA", "cdScaOA", "csAbsOA", "csExtOA", "csScaOA")
    # present <- intersect(nms, names(ld$oa_incidence)) # sometimes only cdExtOA, others not implemented
    lOA <- lapply(nms, function(id) {
      if(!(id %in% names(ld$oa_incidence))) {
        d <- data.frame(wavelength = ld$Wavelengths, NA) } else {
          d <- data.frame(wavelength = ld$Wavelengths, ld$oa_incidence[[id]])
        }
      
      if(ncol(d) > 2) # decomposition in 1..nmax orders
        names(d) = c('wavelength', 'total', paste0('I', seq_len(ncol(d)-2))) else names(d) = c('wavelength', 'total')
        
        d$crosstype <- substr(id,3,5)
        d
    })
    
    names(lOA) <- nms
    lCD <- lapply(lOA[c("cdAbsOA", "cdExtOA", "cdScaOA")], reshape2::melt, 
                  id=c("wavelength","crosstype"))
    lCS <- lapply(lOA[c("csAbsOA", "csExtOA", "csScaOA")], reshape2::melt, 
                  id=c("wavelength","crosstype"))
    mCOA <- purrr::map2_df(lCD, lCS, process_OA)
    
    results <- c(results, list(mCOA=mCOA))
    
    # OA split by particles
    if(("csAbsOA_split" %in% names(ld$oa_incidence)) && !is.null(ld$oa_incidence$csAbsOA_split)){ 
      d <- data.frame(ld$oa_incidence$csAbsOA_split)
      splitOA <- cbind(wavelength = ld$Wavelengths, setNames(d, c('total', paste0('scatterer', seq_len(ncol(d)-1)))))
      
      results <- c(results, list(splitOA=splitOA))
      }
    
    results
  } 
  
  ## Read in fixed orientation data
  
  if("fixed_incidence" %in% names(ld)){
    
    att <- rhdf5::h5readAttributes(hdf5, "Far-Field/fixed_incidence")
    if(!is.null(att) && verbose) message(att)
    
    # circular
    lfCircular <- c("csAbs3R", "csAbs4L", "csExt3R", "csExt4L", "csSca3R", "csSca4L")
    
    lCFO <- lapply(lfCircular, function(n) {
      d <- data.frame(wavelength = ld$Wavelengths, ld$fixed_incidence[[n]])
      names(d) = c('wavelength', 'total', paste0('I', seq_len(ncol(d)-2)))
      d$crosstype <- substr(n,3,5)
      d
    })
    names(lCFO) <- lfCircular
    lR <- lapply(lCFO[c(1,3,5)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    lL <- lapply(lCFO[c(2,4,6)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    mCFO <- purrr::map2_df(lR, lL, process_fixed)
    
    # linear
    lfLinear <- c("csAbs1X", "csAbs2Y", "csExt1X", "csExt2Y", "csSca1X", "csSca2Y")
    
    lLFO <- lapply(lfLinear, function(n) {
      d <- data.frame(wavelength = ld$Wavelengths, ld$fixed_incidence[[n]])
      names(d) = c('wavelength', 'total', paste0('I', seq_len(ncol(d)-2)))
      d$crosstype <- substr(n,3,5)
      d
    })
    names(lLFO) <- lfLinear
    
    lX <- lapply(lLFO[c(1,3,5)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    lY <- lapply(lLFO[c(2,4,6)], reshape2::melt, 
                 id=c("wavelength","crosstype"))
    mLFO <- purrr::map2_df(lX, lY, process_fixed)
    
    results <- c(results, list(mCFO = mCFO, mLFO = mLFO))
    
  }
  
  # add extra metadata if needed
  results <- c(results, ...)
  
}

#' @noRd
#' @export
read_file <- function(x) {
  d = read.table(x)
  names(d) = c('wavelength', 'total', paste0('I', seq_len(ncol(d)-2)))
  d$crosstype <- substr(x,3,5)
  d
}

#' @noRd
#' @export
process_OA <- function(m1, m2){
  m <- m1
  m$dichroism <- m1$value
  m$average <- m2$value
  m$value <- NULL
  m
}

#' @noRd
#' @export
process_fixed <- function(m1, m2){
  m <- m1
  m$polarisation1 <- m1$value
  m$polarisation2 <- m2$value
  m$dichroism <- m1$value - m2$value
  m$average <- 0.5*(m1$value + m2$value)
  m$value <- NULL
  m
}


##' @description Extracts partial absorption cross-sections from a HDF5 file storing far-field cross-sections (Mode=2), and reshapes the data into long-format data.frames suitable for plotting
##' @describeIn store consolidate partial absorption cross-sections for multilayered spheres
##' @param hdf5 filename
##' @param verbose logical: print attributes
##' @param ... extra arguments passed to the final list
##' @export
consolidate_partials <- function(hdf5, verbose = TRUE){
  ld <- rhdf5::h5read(hdf5, "Far-Field")
  
  att <- rhdf5::h5readAttributes(hdf5, "Far-Field/partial_absorption")
  if(!is.null(att) && verbose) message(att)
  
  partials <- ld$partial_absorption
  nms <- names(partials)
  pattern <- 'csAbs(1X|2Y|3R|4L)_scat([0-9]+)coat([0-9]+)'
  proto <- data.frame(polarisation = NA_character_, scatterer = NA_integer_, coat = NA_integer_)
  info <- strcapture(pattern, nms, proto = proto)
  
  info$total <- lapply(partials, function(p) data.frame(wavelength = ld$Wavelengths, total = p[,1]))
  
  # complete <- unnest(info, cols='total') %>% group_by(scatterer, wavelength, polarisation) %>% 
  #   pivot_wider(names_from = coat, values_from = total, names_prefix = 'partial_') %>% 
  #   mutate(core = partial_0, total = partial_1, shell = total - core) %>% 
  #   select(-partial_0, -partial_1) %>% 
  #   pivot_longer(c('total','shell','core'), names_to = 'region') %>% 
  #   ungroup()
  
  average <- unnest(info, cols='total') %>% group_by(scatterer, wavelength, polarisation) %>% 
    pivot_wider(names_from = coat, values_from = total, names_prefix = 'partial_') 
  
  # note: this fails for particles with not exactly 1 shell... so let users postprocess it
  
  # mutate(core = partial_0, total = partial_1, shell = total - core) %>% 
  # select(-partial_0, -partial_1) %>% 
  # ungroup() %>% group_by(scatterer, wavelength) %>% 
  # filter(polarisation %in% c('3R', '4L')) %>% # take both circular
  # summarise(core_avg = mean(core), shell_avg = mean(shell), total_avg = mean(total), .groups = 'drop') %>% 
  # pivot_longer(c('total_avg','shell_avg','core_avg'), names_to = 'region') %>%
  # ungroup()
  
  return(average)
}