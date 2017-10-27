#' Hastings ratio for the Strauss model
#' 
#' Function used for calculating the hastings ratio used for simulation with the Metropolis-Hastings algorithm, \code{\link{rStrauss}}.
#' @param x A point pattern of class \code{\link{pp3}}.
#' @param x_prop A point proposed to be added to the point pattern.
#' @param xi A point in \code{x}.
#' @param h The hardcore parameter.
#' @param beta The intensity parameter.
#' @param gamma The interaction parameter.
#' @param R The interactionradius parameter.
#' @param birth A parameter indicating whether the death or birth Hasting ratio should be calculated.
#' @return A numeric value. The Hastings ratio.
#' @details Used internally.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat

HRStrauss <- function(x, x_prop, xi, h, beta, gamma, R, birth = TRUE){
  Vol <- volume(x$domain)
  np <- npoints.pp3(x)
  if(birth){
    if(np < 2){
      if((np == 1) & (min(nndist(x_prop)) <= h)){
        hr <- 0
      }
      else{
        hr <- beta*gamma^(sum(crossdist.pp3(xi, x) <= R))*Vol/(np + 1)
      }
    }
    else{
      min.dist.x <- min(nndist.pp3(x))
      min.dist.x_prop <- min(nndist.pp3(x_prop))
      not.OK <- ((min.dist.x_prop <= h) & (min.dist.x > h))
      if(not.OK){
        hr <- 0
      }
      else{
        hr <- beta*gamma^(sum(crossdist.pp3(xi, x) <= R))*Vol/(np + 1)
      }
    }
  }
  if(!birth){
    hr <- np / (beta * Vol * gamma^(sum(crossdist.pp3(xi, x_prop) <= R)))
  }
  return(hr)
}
