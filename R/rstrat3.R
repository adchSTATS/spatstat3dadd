#' Stratified random sample of point
#'
#' Generates a stratified random point pattern i a window.
#' @param win A window of class \code{\link{box3}} or any data accepted by \code{\link{as.box3}}
#' @param nx The number of cuboids in the x direction, respectively.
#' @param ny The number of cuboids in the y direction, respectively.
#' @param nz The number of cuboids in the z direction, respectively.
#' @param k Number of points generated in each cuboid
#' @return A \code{\link{data.frame}} with 3 columns and \code{nx*ny*nz*k} rows
#' @details The window is divided into a number of cuboid tiles and \code{k} points are then uniformly sampled in each cuboid
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat
#' @export

rstrat3 <- function(win, nx, ny=nx, nz=nx, k=1){
  win <- as.box3(win)

  ## Test certain conditions necessary for the stratification to make sence
  stopifnot(nx >= 1 && ny >= 1 && nz >= 1)
  stopifnot(k >= 1)

  ## Generate centre point of the desired tiles
  centres <- grid3(win, nx, ny, nz)

  ## Generate multiple centre points in each tile if that is desired
  xcent <- rep.int(centres$x, k)
  ycent <- rep.int(centres$y, k)
  zcent <- rep.int(centres$z, k)

  ## Define the dimmensions of the tiles
  width <- diff(win$xrange)/nx
  height <- diff(win$yrange)/ny
  depth <- diff(win$zrange)/nz

  ## Simulate a point in the tile
  ## This is done by simulating a point around the centre point that is contained within the width/height/depth of the tiles
  n <- nx * ny * nz * k
  x <- xcent + stats::runif(n, min = -width/2, max = width/2)
  y <- ycent + stats::runif(n, min = -height/2, max = height/2)
  z <- zcent + stats::runif(n, min = -depth/2, max = depth/2)

  return(pp3(x, y, z, win))
}

