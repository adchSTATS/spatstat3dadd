#' Erosion
#'
#' Function for eroding a cuboid by some value r
#' @param win A window of class \code{\link{box3}} or any data accepted by \code{\link{as.box3}}.
#' @param r The erosion factor.
#' @return A window of class \code{\link{box3}}.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat
#' @export

## Function for eroding a window of object class box3 og any data accepted by as.box3
erosion3 <- function(win, r){
  w <- as.box3(win)
  xrange_new <- c(w$xrange[1] + r, w$xrange[2] - r)
  yrange_new <- c(w$yrange[1] + r, w$yrange[2] - r)
  zrange_new <- c(w$zrange[1] + r, w$zrange[2] - r)
  wnew <- box3(xrange_new, yrange_new, zrange_new)
  return(wnew)
}
