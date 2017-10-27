#' J-function
#'
#' Function for estimating the J-function from a three dimensional point pattern.
#' @param X An observed point pattern of class \code{\link{pp3}}.
#' @param ... Ignored.
#' @return An object of class \code{\link{fv.object}} containing r the vector of radii at which the J-function is to be estimated.
#' @author Andreas Dyreborg Christoffersen \email{andreas@math.aau.dk}
#' @import spatstat genfun
#' @export

J3est <- function(X, ...){
  stopifnot(inherits(X, "pp3"))
  F3est = F3est(X, ...)
  G3est = G3est(X, ...)
  rvals = F3est$r
  rmax = max(rvals)
  J = fv(data.frame(r = rvals, theo = 1), "r", substitute(J(r), NULL),
         "theo", . ~ r, c(0, rmax), c("r", "%s[pois](r)"),
         c("distance argument r", "theoretical Poisson %s"), fname = "J")
  Fnames <- names(F3est)
  Gnames <- names(G3est)
  if ("raw" %in% Gnames && "raw" %in% Fnames) {
    Jun = ratio(1 - G3est$raw, 1 - F3est$raw)
    J = bind.fv(J, data.frame(un = Jun), "hat(%s)[un](r)",
                "uncorrected estimate of %s", "un")
    attr(J, "alim") = range(rvals[F3est$raw <= 0.9])
  }
  if ("rs" %in% Gnames && "rs" %in% Fnames) {
    Jrs = ratio(1 - G3est$rs, 1 - F3est$rs)
    J = bind.fv(J, data.frame(rs = Jrs), "hat(%s)[rs](r)",
                "border corrected estimate of %s", "rs")
    attr(J, "alim") = range(rvals[F3est$rs <= 0.9])
  }
  if ("han" %in% Gnames && "cs" %in% Fnames) {
    Jhan = ratio(1 - G3est$han, 1 - F3est$cs)
    J = bind.fv(J, data.frame(han = Jhan), "hat(%s)[han](r)",
                "Hanisch-style estimate of %s", "han")
    attr(J, "alim") = range(rvals[F3est$cs <= 0.9])
  }
  if ("km" %in% Gnames && "km" %in% Fnames) {
    Jkm = ratio(1 - G3est$km, 1 - F3est$km)
    J = bind.fv(J, data.frame(km = Jkm), "hat(%s)[km](r)",
                "Kaplan-Meier estimate of %s", "km")
    attr(J, "alim") = range(rvals[F3est$km <= 0.9])
  }
  if ("hazard" %in% Gnames && "hazard" %in% Fnames) {
    Jhaz = G3est$hazard - F3est$hazard
    J = bind.fv(J, data.frame(hazard = Jhaz), "hazard(r)",
                "Kaplan-Meier estimate of derivative of log(%s)")
  }
  nama = names(J)
  fvnames(J, ".") = rev(nama[!(nama %in% c("r", "hazard"))])
  attr(J, "F") = F3est
  attr(J, "G") = G3est
  unitname(J) <- unitname(F3est)
  return(J)
}
