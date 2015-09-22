#' Angle constrained MCR-ALS.
#'
#' @param data
#' a matrix with spectral values
#' @param  ncomp
#' number of pure components to identify
#' @param do.angle
#' logical, should the angle constraint be used
#' @param f
#' how much of mean should be added for angle constraint
#' @param  wavelength
#' a vector with wavelength if necessary
#' @param  max.niter
#' maximum number of iterations
#' 
#' @return
#' The method returns results of unmixing as an object of \code{mcrals} class with following fields.
#' \item{spec }{resolved spectra of pure components.}
#' \item{conc }{resolved concentrations of pure components.}
#' \item{ncomp }{number of resolved components.}
#' \item{resvar }{residual variance for each component.}
#' 
#' @details 
#' The method implements a simple MCR-ALS algorithm with non-negativity constraints applied both to
#' spectra and contributions by setting all negative values to zero. In addition an angle constraint
#' can be applied by adding a small portion (default 5%) of mean to spectra or contributions [1].
#' 
#' @seealso 
#' The class has numerous methods for investigating and refininf the unmixing results.
#' \tabular{ll}{
#'  \code{print} \tab prints information about a \code{purity} object.\cr
#'  \code{summary} \tab shows summary statistics for the object.\cr
#'  \code{plot} \tab plot overview of unmixing results.\cr
#'  \code{\link{plotSpectra}} \tab shows a plot with resolved spectra.\cr
#'  \code{\link{plotConcentrations}} \tab shows a plot with resolved concentrations.\cr
#'  \code{\link{plotVariance}} \tab shows a plot with residual variance vs. number of components.\cr
#'  \code{\link{plotSpectraComparison}} \tab compares resolved and ethalon spectra if available.\cr
#'  \code{\link{image.mcrals}} \tab if data is a hypercube, shows resolved concetration map as an image.\cr
#' }
#' 
#' @references 
#' 1. W. Windig, M.R. Keenan, Appl. Spectrosc. 65 (2011) 349–357.
#' 
#' @examples
#' 
#' # unmix hyperspectral data with and without angle constraint and compare the 
#' # results for first two components
#' data(puredata)
#'
#' res1 = mcrals(hsi$D, 3, do.angle = F, wavelength = hsi$wavelength)
#' res2 = mcrals(hsi$D, 3, do.angle = T, wavelength = hsi$wavelength)
#' 
#' summary(res1)
#' summary(res2)
#'
#' # because initial estimates are made with random numbers
#' # the order of components in the results is different  
#' par(mfrow = c(2, 1))
#' plotSpectra(res1)
#' plotSpectra(res2)
#' par(mfrow = c(1, 1))
#' 
#' par(mfrow = c(3, 2))
#' for (i in 1:3)
#' {
#'  image(res1, i)
#'  image(res2, i)
#' }  
#' par(mfrow = c(1, 1))
#' 
#' @export
mcrals = function(D, ncomp, Ct = NULL, St = NULL, do.angle = T, f = 0.05, wavelength = NULL, 
                  max.niter = 500)
{
  nvar = ncol(D)
  
  # generate wavelength if not provided  
  if (is.null(wavelength))
  {
    wavelength = 1:ncol(D)
  }  
  
  # make initial estimate for pure spectra if they are not provided
  if (is.null(St) && is.null(Ct))
  {  
    St = matrix(runif(nvar * ncomp), ncol = ncomp, nrow = nvar)
  }  
  else if (!is.null(Ct))
  {
    St = t(D) %*% Ct %*% pinv(t(Ct) %*% Ct)
    St[St < 0] = 0
  }  
  
  resvarold = 1
  Ct.old = NULL
  St.old = NULL
  
  for (i in 1:max.niter)
  {
    # normalize St
    len = sqrt(colSums(St^2))
    St = sweep(St, 2, len, '/')

    # apply angle constraint    
    if (do.angle)
    {
      m = apply(St, 1, mean)
      St = sweep((1 - f) * St, 1, f * m, '+')
    }  
    
    # calculate Ct and apply non-negativity
    Ct = D %*% St %*% pinv(t(St) %*% St)
    Ct[Ct < 0] = 0
    
    # calculate St and apply non-negativity
    St = t(D) %*% Ct %*% pinv(t(Ct) %*% Ct)
    St[St < 0] = 0
    
    # residual variance
    Dt = Ct %*% t(St)
    res = D - Dt
    resvar = sum(res^2) / sum(D^2)
    
    if (abs(resvar - resvarold) < 0.00001 && i > 5)
    {
      break
    }
    
    resvarold = resvar
    Ct.old = Ct
    St.old = St
  } 
  
  # calculate residual variance for each component
  resvar = rep(0, ncomp)
  for (i in 1:ncomp)
  {  
    Dt = Ct[, 1:i] %*% t(St[, 1:i])
    res = D - Dt
    resvar[i] = sum(res^2) / sum(D^2)
  }
  
  names(resvar) = colnames(Ct.old) = colnames(St.old) = paste('C', 1:ncol(Ct.old), sep = '')
  rownames(Ct.old) = rownames(D)
  rownames(St.old) = colnames(D)

  attr(Ct.old, 'width') = attr(D, 'width')
  attr(Ct.old, 'height') = attr(D, 'height')
  
  res = list(
    spec = St.old,
    conc = Ct.old,
    wavelength = wavelength,
    ncomp = ncomp,
    resvar = resvar,
    do.angle = do.angle,
    f = f
  )
  
  res$call = match.call()   
  class(res) = "mcrals"
  res
}  

#' Image with concentration map
#' 
#' @description 
#' Make an image with map of resolved concentrations 
#' @param obj
#' object with mcrals results
#' 
#' @details 
#' The method is similar to \code{\link{image.purity}}.
#' 
#' @export
image.mcrals = function(obj, ...)
{
  image.purity(obj, ...)
}  


#' Overview plot for mcrals object
#' 
#' @description
#' Shows a set of plots for MCR-ALS results.
#' 
#' @param x
#' a  (object of class \code{mcrals})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{mcrals}} function.
#' 
#' @export
plot.mcrals = function(x, ...)
{   
  layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = T))
  plotSpectra(x)
  plotConcentrations(x)
  plotVariance(x)
  par(mfrow = c(1, 1))
}

#' Print method for mcrals object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' unmixing results (object of class \code{mcrals})
#' @param ...
#' other arguments
#'
#' @export 
print.mcrals = function(x, ...)
{
  cat("\nUnmixing results (class mcrals)\n")
  cat('\nCall:\n')
  print(x$call)
  cat('\nMajor fields:\n')
  cat("$spec - resolved spectra of pure components\n")
  cat("$conc - resolved concentrations of pure components\n")
  cat("$ncomp - number of resolved components\n")
  cat("$resvar - vector with residual variance for resolved components\n")
  cat("$wavelength - vector with wavelength/wavebands/wavenumbers\n")
  cat('\nTry summary(obj) and plot(obj) to see the results.\n')   
}

#' Summary method for mcrals object
#' 
#' @description
#' Shows summary statistics for mcrals objects.
#' 
#' @param object
#' unmixing results (object of class \code{mcrals})
#' @param ...
#' other arguments
#' 
#' @export 
summary.mcrals = function(object, ...)
{
  cat('\nSummary statistics for unmixing results:\n')
  cat(sprintf('Number of resolved components: %d\n', object$ncomp))
  cat(sprintf('Total explained variance: %.1f%%\n', (1 - object$resvar[object$ncomp]) * 100))
  cat(sprintf('Use of angle constraint: %d\n', object$do.angle))
  cat(sprintf('The amount of mean added: %d%%\n', object$f * 100))
  cat('\n')
}