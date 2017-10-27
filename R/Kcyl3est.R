#' Cylindrical K-function
#'
#' Function for estimating the cylindrical K-function from a three dimensional point pattern.
#' @param X A point pattern of class \code{\link{pp3}}.
#' @param t A numeric vector of heights.
#' @param u A numeric vector determining the direction of the cylinders.
#' @param ... Ignored
#' @param rmax Optional. Maximum value of argument \code{r} for which \code{Kcyl3(r)} will be estimated.
#' @param nrval Optional. Number of values of r for which \code{Kcyl3(r)} will be estimated.
#' A large value of nrval is required to avoid discretisation effects.
#' @param drop A logical value indicating whether or not the object to be returned should be a list or an \code{\link{fv}} oject.
#' Only works when \code{length(t) == 1}.
#' @return A list where each element is an object of class \code{\link{fv}} containing the \code{r} values,
#' theoretical Poisson values and the estimated function values.
#' Each list element corresponds to the different heights of the cylinders.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat genfun
#' @export

Kcyl3est <- function(X, t, u, ..., rmax = NULL, nrval = 128, drop = TRUE){
  inherits(X, "pp3")

  if (is.null(rmax))
    rmax <- diameter(X$domain)/2
  r <- seq(from = 0, to = rmax, length.out = nrval)

  out <- list()
  for(i in 1:length(t)){
    out[[i]] <- fv(x = data.frame(r = r, theo = 2 * pi * r^3 * t[i]),
                   argu = "r",
                   ylab = substitute(Kcyl(r), NULL),
                   valu = "theo",
                   fmla = . ~ r,
                   alim = c(0, max(r)),
                   labl = c("r", "%s[theo](r)"),
                   desc = c("distance argument r", "Theoretical %s for homogenious Poisson"),
                   fname = c("Kcyl", "3"))
  }
  names(out) <- paste("t=", t, sep = "")

  u <- u/sqrt(crossprod(u))

  W <- stats::window(X)$domain
  win <- c(W$xrange[2], W$yrange[2], W$zrange[2])
  win.vol <- volume(W)
  np <- npoints(X)
  rho.sq <- np * (np - 1) / win.vol^2 #find ud af hvordan vi vil estimere rho^2
  dat <- as.data.frame(X$data)

  ind.list <- whichdiffincyl(as.matrix(dat), u, r^2, t)
  nams <- names(ind.list)
  nams <- matrix(gsub("\\.", "", gsub("r", "", unlist(strsplit(nams, split = "t")))), ncol = 2, byrow = TRUE)
  class(nams) <- "numeric"
  colnames(nams) <- c("r", "t")

  empty.ind <- which(sapply(ind.list, function(x){nrow(x) == 0}))

  Kcyl_mat <- matrix(nrow = length(r), ncol = length(t))
  Kcyl_mat[nams[empty.ind, ]] <- 0
  Kcyl_mat <- data.frame(Kcyl_mat)

  if(length(empty.ind) > 0){
    j.ind <- (1:length(ind.list))[-c(empty.ind)]
  }else j.ind <- 1:length(ind.list)

  for(j in j.ind){
    vec <- abs(dat[ind.list[[j]]$index2, ] - dat[ind.list[[j]]$index1, ])
    edge.cor <- apply(sweep(-vec, 2, win, "+"), 1, prod)
    Kcyl_mat[nams[j, 1], nams[j, 2]] <- 2*sum(1/edge.cor)/rho.sq
    colnames(Kcyl_mat)[nams[j, 2]] <- paste("Kcyl", t[nams[j, 2]], sep = "")
  }

  for(i in 1:length(t)){
    tmp_df <- data.frame(Kcyl_mat[paste("Kcyl", t[i], sep = "")])
    names(tmp_df) <- "Kcyl"
    tmp <- bind.fv(x = out[[paste("t=", t[i], sep = "")]],
                   y = tmp_df,
                   labl = paste("hat(%s)[est](r) , ", "t=", t[i], sep = ""),
                   desc = "Estimate of %s",
                   preferred = "Kcyl")
    out[[paste("t=", t[i], sep = "")]] <- tmp
  }

  if(drop == TRUE & length(t) == 1){
    out <- out[[1]]
  }
  return(out)
}
