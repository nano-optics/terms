## display functions

# n * (n+1) + n = pmax
# n^2 + 2n - pmax = 0
# (-b+sqrt(b^2-4*a*c)) / (2*a)
# pmax=264
# (-2+sqrt(2+4*pmax/2)) / (2)


##' @description Combined T-matrix index
##' @describeIn utility p-index
##' @param in1 index
##' @param in2 index
##' @export
p_index <- function(in1,in2){
  # % (n,m) -> p = n * (n+1) + m
  in1 * (in1 + 1) + in2
}

##' @description T-matrix indices
##' @describeIn utility indices
##' @param n_max maximum single-particle index
##' @param n_part number of particles
##' @export
indices <- function(n_max=3, n_part=1){
  un <- seq(1,n_max)
  
  n = do.call(c, lapply(un, function(.n)rep(.n, 2*.n+1)))
  m = do.call(c, lapply(un, function(.n)seq(-.n,.n)))
  
  nn <- expand.grid(n = n, np = n)
  mm <- expand.grid(m = m, mp = m)
  d <- rbind(cbind(nn, mm, q=1,qp=1),
             cbind(nn, mm, q=1,qp=2),
             cbind(nn, mm, q=2,qp=1),
             cbind(nn, mm, q=2,qp=2))
  
  dplyr::mutate(d, 
                p = n*(n+1)+m, 
                l = (q - 1)* max(p) + p,
                i = (n_part-1) * max(l) + l,
                pp = np*(np+1)+mp,
                lp = (qp - 1)* max(pp) + pp,
                j = (n_part-1) * max(lp) + lp)
}


##' @description Unpack T-matrix indices
##' @describeIn utility unpack indices
##' @param n_max maximum single-particle index
##' @param j_max size of collective T-matrix
##' @export
unpack_indices <- function(n_max=3, j_max=2){
  
  lmax <- 2*(n_max*(n_max+1)+n_max)
  imax <- j_max*lmax
  
  dplyr::mutate(data.frame(i=seq(1, imax)),
                j = (i-1) %/% (imax/j_max) + 1,
                l = i - (j-1)*imax %/% j_max,
                q = (2*(l-1))%/%lmax + 1,
                p=l - (q-1)*lmax/2,
                n = floor(sqrt(p)),
                m = p - n*(n+1) )
}



##' @description Unpack T-matrix indices
##' @describeIn utility unpack indices
##' @param f filename
##' @export
read_amat <- function(f){
  size <- read.table(f, header = F, nrows = 1)
  v <- read.table(f, header = F, skip=1)
  names(v) <- c('i', 'ip', 'Tr', 'Ti')
  
  invisible(list(v = v, size=unlist(size)))
}


##' @description Wrap staged matrices into a list of T-matrix like objects
##' @describeIn utility wrap staged matrices
##' @param a A-matrix object (from read_amat)
##' @param n_max maximum single-particle index
##' @param n_part number of particles
##' @export
amat_to_tmatlist <- function(a, n_max=3, n_part=2){
  rows <- unpack_indices(n_max, j_max=n_part)
  cols <- setNames(rows, paste0(names(rows),'p'))
  tmatbig <- merge(merge(a, rows, by = 'i'), cols, by = 'ip')
  split(tmatbig, f = tmatbig$j)
}



##' @description Read T-matrix into long-format data.frame
##' @describeIn utility read T-matrix
##' @param f filename
##' @param save store result as Rds file
##' @export
read_tmat <- function(f, save=FALSE){
  con <- pipe(paste("grep -e ^#", f))
  comments <- readLines(con)
  close(con)
  v <- read.table(f, header = F, comment.char = '#')
  names(v) <- strsplit(split = ' ', comments[1])[[1]][2:9]
  rowcom <- strsplit(comments[-1], '\\s+')
  wavelength <- as.numeric(lapply(rowcom, "[",3))
  Nb <- as.numeric(lapply(rowcom, "[",5))
  epsilon <- as.complex(gsub("j","i",lapply(rowcom, "[",7)))
  meta <- setNames(data.frame(wavelength, Nb, epsilon), 
                   c('wavelength',"Nb","epsilon"))
  meta$wavelength <- as.numeric(meta$wavelength)
  
  replambda <- rep(meta$wavelength, each = unique(meta$Nb))
  
  stopifnot(length(replambda) == nrow(v))
  
  v$wavelength <- replambda
  v$p <- p_index(v$n, v$m)
  v$pp <- p_index(v$np, v$mp)
  
  if(save) saveRDS(v, file=gsub('.txt$','.Rds',f))
  invisible(v)
}




##' @description Display staged matrix
##' @describeIn utility display staged matrix
##' @param l list of T-matrices, from amat_to_tmatlist
##' @import ggplot2 dplyr purrr tidyr
##' @export
display_amat <- function(l){
  
  s <- l$v
  s$modT <- Mod(complex(s$Tr,s$Ti))
  
  maxi <- l$size[1] #max(s$i)
  maxip <- l$size[2]  #max(s$ip)
  dd <- expand.grid(ip=1:maxip,i=1:maxi)
  
  
  ggplot(s, aes(ip,i, fill=modT)) +
    geom_tile(data=dd, fill=NA) +
    geom_tile() +
    scale_y_reverse(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0)) + 
    coord_equal()+
    # scale_fill_continuous(low=blues9[1],high=blues9[9],
    #                       lim=c(0, max(s$modT)), expand=c(0,0)) +
    scale_fill_viridis_c(option='magma',direction = -1,lim=c(0, max(s$modT)), expand=c(0,0)) +
    theme_void() + guides(fill='none') +
    theme(panel.border = element_rect(fill=NA, colour='white'),
          panel.background = element_rect(fill='grey96', colour=NA),
          plot.margin = margin(4,4,4,4),
          # plot.background = element_rect(fill='grey70'),
          # panel.grid.major = element_line(color = 'white', size = 0.5),
          # panel.grid.minor = element_line(color = 'white', size = 0.2),
          strip.text = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5), panel.ontop = FALSE)
  
}


##' @description Display T-matrix
##' @describeIn utility display T-matrix
##' @param s T-matrices, from read_tmat
##' @export
display_tmat <- function(s){
  
  # currently using s and not q in smarties
  if(!("q" %in% names(s))){
    s$q <- s$s
    s$qp <- s$sp
  }
  
  s$qqp <- paste0(s$q,s$qp)
  
  s$modT <- Mod(complex(s$Tr,s$Ti))
  # dput(RColorBrewer::brewer.pal(4,"Pastel1"))
  dummy <- data.frame(qqp = c(11,12,21,22),fill=c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4"))
  
  maxp <- max(s$p)
  maxpp <- max(s$pp)
  dd <- expand.grid(p=1:maxp, pp=1:maxpp)
  # print(str(dd))
  nn <- unique(s$n)
  major <- cumsum(sapply(nn, function(n) length(seq(-n,n)))) + 0.5
  minor <- seq(0, maxp-1) +0.5
  
  ggplot(s, aes(p,pp, fill=modT)) +
    facet_wrap(~qqp, ncol=2) +
    geom_tile(data=dd, fill=NA) +
    geom_tile() +
    geom_vline(data=data.frame(x=major), aes(xintercept=x), 
               size=0.6, colour='white')+
    geom_vline(data=data.frame(x=minor), aes(xintercept=x), 
               size=0.2, colour='white',lty=1)+
    geom_hline(data=data.frame(y=major), aes(yintercept=y), 
               size=0.6, colour='white')+
    geom_hline(data=data.frame(y=minor), aes(yintercept=y), 
               size=0.2, colour='white',lty=1)+
    scale_y_reverse(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0)) + 
    coord_equal()+
    # scale_fill_continuous(low=blues9[1],high=blues9[9],
    #                       lim=c(0, max(s$modT)), expand=c(0,0)) +
    scale_fill_viridis_c(option='magma',direction = -1,lim=c(0, max(s$modT)), expand=c(0,0)) +
    theme_void() + guides(fill='none') +
    theme(panel.border = element_rect(fill=NA, colour='white'),
          panel.background = element_rect(fill='grey96', colour=NA),
          plot.margin = margin(4,4,4,4),
          # plot.background = element_rect(fill='grey70'),
          # panel.grid.major = element_line(color = 'white', size = 0.5),
          # panel.grid.minor = element_line(color = 'white', size = 0.2),
          strip.text = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5), panel.ontop = FALSE)
  
}



##' @description Display prestaged matrix
##' @describeIn utility display prestaged matrix
##' @param a prestaged matrix, from read_amat
##' @param n_max maximum order
##' @param n_part number of particles
##' @param draw logical, draw output
##' @export
display_prestaged <- function(a, n_max=3, n_part=2, draw=TRUE){
  
  tmat_list <- amat_to_tmatlist(a, n_max, n_part)
  lay <- matrix(NA, length(tmat_list), length(tmat_list))
  diag(lay) <- seq_along(tmat_list)
  
  display_one <- function(t){
    display_tmat(t)  + theme(panel.border = element_rect(fill=NA, colour='white'),
                             panel.background = element_rect(fill='grey96', colour=NA),
                             plot.margin = margin(0,0,0,0),
                             axis.title = element_blank(),
                             axis.text = element_blank(),
                             strip.text = element_blank(),
                             plot.title = element_text(hjust = 0.5), 
                             panel.ontop = FALSE,
                             plot.background = element_rect(fill=NA))
  }
  
  # library(grid)
  gl <- lapply(tmat_list, display_one)
  g <- gridExtra::arrangeGrob(grobs=gl, 
                              layout_matrix = lay, 
                              vp=grid::viewport(width=grid::unit(1,"snpc"), 
                                                height=grid::unit(1,"snpc")))
  if(draw) grid::grid.draw(g)
  invisible(g)
}


##' @noRd
read_vtacs <- function(f= 'debug_vtacs_h.txt'){
  d <- read.table(f, skip=1)
  m <- matrix(0, max(d$V1), max(d$V2))
  m[cbind(d$V1,d$V2)] <- complex(real = d$V3, imaginary = d$V4)
  
  m
}