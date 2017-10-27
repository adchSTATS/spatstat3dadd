#' Cuboid grid of points
#'
#' Generate a cuboid grid of points in three dimensions window
#' @param win A window of class \code{\link{box3}} or any data accepted by \code{\link{as.box3}}.
#' @param nx Number of points in the x direction of the grid.
#' @param ny Number of points in the y direction of the grid.
#' @param nz Number of points in the z direction of the grid.
#' @return A list consisting of x, y and z which are numeric vectors constituting the coordinates of the cuboid grid.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat
#' @export
grid3 <- function(win, nx, ny, nz){
  win <- as.box3(win)

  ## Extract the ranges for the 3 dimensions
  xrange <- win$xrange
  yrange <- win$yrange
  zrange <- win$zrange

  ## Define the nx, ny and nz centre values for x, y respectively z
  x <- seq(from = xrange[1], to = xrange[2], length.out = 2 * nx + 1)[2 * (1:nx)]
  y <- seq(from = yrange[1], to = yrange[2], length.out = 2 * ny + 1)[2 * (1:ny)]
  z <- seq(from = zrange[1], to = zrange[2], length.out = 2 * nz + 1)[2 * (1:nz)]

  # # ## Properly replicate the above values to obtain a grid
  x <- rep.int(x, ny*nz)
  y <- rep.int(rep.int(y, rep.int(nx, ny)), nz)
  z <- rep.int(z, rep.int(nx*ny, nz))

  #out <- expand.grid(x, y, z)

  return(data.frame(x = x, y = y, z = z))
}
