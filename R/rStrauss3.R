#' Simulate Strauss
#'
#' Function for simulating the Strauss (hardcore) model a number of times.
#' @param win A window in which the simulation should be contained. Should be of class \code{\link{box3}}
#' or any data accepted by \code{\link{as.box3}}.
#' @param beta The intensity parameter (a positive number).
#' @param gamma The interaction parameter (a number between 0 and 1, both inclusive).
#' @param R The interactionradius parameter (a number greater than \code{h}).
#' @param h The hardcore parameter (a non-negative number).
#' @param multi.sim If \code{TRUE} a list of simulations will be returned. 
#' The "distance" between each simulated point pattern will be \code{1000}.
#' @param nit The number of iterations/point pattern to be simulated.
#' @param trace If \code{TRUE} the number of points, number of h-close neighbors, and number of R-close neighbors is printed.
#' The number of points and the number of R-close neighbors are the sufficient statistics.
#' @param Keep.an.eye If \code{TRUE} a number indicating the number of poit patterns so far simulated will be printed.
#' @return A list of three dimensional point patterns
#' @details The initial point pattern is an empty Poisson point pattern.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat
#' @export

rStrauss3 <- function(beta, gamma = 1, R = h, h = 0, win = box3(), multi.sim = FALSE, nit = 10000, trace = TRUE, Keep.an.eye = FALSE){
  Y <- rpoispp3(lambda = 0, domain = win, nsim = 1, drop = TRUE)
  if(multi.sim){
    X <- list()
  }
  when <- seq(1000, nit, by = 1000)
  if(trace){
    np <- npoints.pp3(Y)
    hclose <- sum(nndist.pp3(Y) < h)
    S_R <- (sum(pairdist.pp3(Y) <= R) - npoints(Y))/2
  }
  for(i in 2:nit){
    R1 <- stats::runif(1, 0, 1)
    R2 <- stats::runif(1, 0, 1)
    if(R1 <= 0.5){
      xi <- runifpoint3(1, win)
      if(npoints.pp3(Y) == 0){
        Y_prop <- pp3(x = xi$data$x$x, y = xi$data$y$y, z = xi$data$z$z, Y$domain)
      }
      else{
        tmp_df <- rbind.hyperframe(Y$data, xi$data)
        Y_prop <- pp3(x = tmp_df$x, y = tmp_df$y, z = tmp_df$z, Y$domain)
      }
      hr <- HRStrauss(x = Y, x_prop = Y_prop, xi = xi, h = h, beta = beta, gamma = gamma, R = R, birth = TRUE)
      a <- min(1, hr)
      if(R2 <= a){
        Y <- Y_prop
      }
    } 
    else{
      if(npoints.pp3(Y) == 0){
        if(trace){
          np <- c(np, npoints.pp3(Y))
          hclose <- c(hclose, sum(nndist.pp3(Y) < h))
          S_R <- c(S_R, (sum(pairdist.pp3(Y) <= R) - npoints(Y))/2)
        }
        if(Keep.an.eye){
          cat("It:", i, "\t", "npoints:", np[i], "\t", "nnpoints:",
              hclose[i], "\t", "S_R:", S_R[i], "\n")
        }
        next
      }
      else{
        point.index = as.numeric(sample(seq(1, npoints.pp3(Y), 1), 1))
        tmp_df <- data.frame(Y$data[-point.index, ])
        Y_prop <- pp3(x = tmp_df$x, y = tmp_df$y, z = tmp_df$z, Y$domain)
        tmp_df <- data.frame(Y$data[point.index, ])
        xi <- pp3(x = tmp_df$x, y = tmp_df$y, z = tmp_df$z, Y$domain)
        hr <- HRStrauss(x = Y, x_prop = Y_prop, xi = xi, h = h, beta = beta, gamma = gamma, R = R, birth = FALSE)
        a <- min(1, hr)
        if(R2 <= a){
          Y <- Y_prop
        }
      }
    }
    if(multi.sim & (i %in% when)){
      X[[i / 1000]] <- Y
    }
    if(trace){
      np <- c(np, npoints.pp3(Y))
      hclose <- c(hclose, sum(nndist.pp3(Y) < h))
      S_R <- c(S_R, (sum(pairdist.pp3(Y) <= R) - npoints(Y))/2)
    }
    if(Keep.an.eye){
      cat("It:", i, "\t", "npoints:", np[i], "\t", "nnpoints:",
          hclose[i], "\t", "S_R:", S_R[i], "\n")
    }
  }
  if(multi.sim){
    Y <- X
  }
  out <- list("Three dimensional point pattern" = Y, "Number of points" = np,
              "Number of points NOT fulfilling hardcore" = hclose,
              "Number of R-close pairs of points" = S_R)
  if(length(out) == 1){
    out <- drop(out[[1]])
  }
  return(out)
}
