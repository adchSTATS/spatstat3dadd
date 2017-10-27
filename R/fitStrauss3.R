#' Fit Strauss model
#'
#' Function for fitting the Strauss model using logistic regression to 3D data.
#' @param X Three dimensional point pattern of class \code{\link{pp3}} contained in a window of class \code{\link{box3}}.
#' @param R A vector of interaction radii
#' @param nd For dum.sim.method = "Poisson" nd^3 is the expected number of dummy points in window eroded by.
#' For dum.sim.method = "stratified" nd specifies the grid dimension (grid is on the eroded window).
#' @param erosion.factor A number defining how much the window shall be eroded in the border correction.
#' If \code{0} no correction will be done. The default is \code{max(R)}.
#' @param dum.sim.method A string indicating which model the dummy point should be simulated from.
#' Options \code{"Poisson"} and \code{"stratified"}
#' @return A \code{\link{data.frame}} with \code{length(R)} row with variables \code{R}, \code{beta}, \code{gamma} and \code{logLik}
#' @details The method used is only for Gibbs point process models on the exponential family form.
#' For the Strauss process this means that the interaction radius should be fixed.
#' When these requirements are satisfied, a logistic regression is fitted.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat
#' @export

fitStrauss3 <- function(X, R, nd, erosion.factor = max(R), dum.sim.method = "stratified"){
  stopifnot(inherits(X, "pp3"))
  
  ## Define a new point pattern with the eroded window as domain and subset points to eroded window
  ero.X <- X
  ero.X.domain <- erosion3(X$domain, erosion.factor)
  ero.X$domain <- ero.X.domain
  df_tmp <- ero.X$data
  df_tmp <- df_tmp[((ero.X.domain$xrange[1] <= df_tmp$x) & (ero.X.domain$xrange[2] >= df_tmp$x) &
                    (ero.X.domain$yrange[1] <= df_tmp$y) & (ero.X.domain$yrange[2] >= df_tmp$y) &
                    (ero.X.domain$zrange[1] <= df_tmp$z) & (ero.X.domain$zrange[2] >= df_tmp$z)) ,]
  ero.X <- pp3(df_tmp$x, df_tmp$y, df_tmp$z, ero.X.domain)
  dum.intensity <- nd^3/volume.box3(ero.X.domain)
  ## Simulate the dummy points
  if(dum.sim.method == "Poisson"){
    dum.points <- rpoispp3(dum.intensity, ero.X.domain)
  } else if(dum.sim.method == "stratified"){
    dum.points <- rstrat3(ero.X.domain, nd, nd, nd, k=1)
  } else stop("Method for simulating dummy points not implemented")

  ## Create a pp3 that appends the dummy points to the observed points and indicate the point origin (dummy or observed)
  ## We will use this for fitting the logistic regression
  ## The inero variable will indicate wheather the points are in the eroded window or not
  X$data$resp <- 1
  dum.points$data$resp <- 0
  new.pp <- pp3(c(X$data$x, dum.points$data$x),
                c(X$data$y, dum.points$data$y),
                c(X$data$z, dum.points$data$z),
                X$domain)
  new.pp$data$resp <- c(X$data$resp, dum.points$data$resp)
  new.pp$data$inero <- (new.pp$data$x >= min(ero.X.domain$xrange) & new.pp$data$x <= max(ero.X.domain$xrange) &
                          new.pp$data$y >= min(ero.X.domain$yrange) & new.pp$data$y <= max(ero.X.domain$yrange) &
                          new.pp$data$z >= min(ero.X.domain$zrange) & new.pp$data$z <= max(ero.X.domain$zrange))

  ## Calculate the independent variable
  ## Find all the parwise distances and exclude the columns that does not correspond to a point from the observed point process
  dist.mat <- pairdist.pp3(new.pp)[, (new.pp$data$resp == 1)]

  ## Define the offset to be use in the logistic regression
  new.pp$data$os <- log(1/dum.intensity)

  ## Count for each row the number of columns that satisfies that the value is smaller than R and append to data
  ## Do the logistic regression
  out <- data.frame(R)
  for(i in 1:length(R)){
    new.pp$data$var <- rowSums(dist.mat <= R[i]) - new.pp$data$resp
    fit <- stats::glm(resp ~ 1 + var + offset(os), family = stats::binomial(link = "logit"),
               data = as.data.frame(new.pp$data[new.pp$data$inero, ]))
    out[i, "beta"] <- exp(fit$coefficients[1])
    out[i, "gamma"] <- exp(fit$coefficients[2])
    out[i, "loglik"] <- stats::logLik(fit)
  }

  return(out)
}


