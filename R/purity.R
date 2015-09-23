#' Spectral unmixing with pure variables
#'
#' @description 
#' \code{purity} is a class for resolving pure spectra and contributions from spectra of mixtures
#' using pure variable approach.
#' 
#' @param spectra
#' a matrix with spectral data
#' @param  ncomp
#' number of pure components to identify
#' @param  offset
#' offset in percent for calculation of parameter alpha
#' @param by.spec
#' logical, should the algorithm works with spectra (rows) or with concentrations (columns)
#' @param use.deriv
#' shall a second derivative be used (0 - no, 1 - only for unmixing, 2 - both for locating pure 
#' variables and for unmixing) 
#' @param  savgol
#' vector of parameters for Savitzky-Goley filter: derivative order, width of filter and polynomial order
#' @param exclude
#' vector of column numbers that should be excluded when pure variables are being detected
#' @param  wavelength
#' a vector with wavelength if necessary
#' 
#' @return
#' The method returns an object of \code{purity} class - a list with with following fields:
#' \item{spec }{resolved spectra of pure components.}
#' \item{conc }{resolved concentrations of pure components.}
#' \item{ncomp }{number of resolved components.}
#' \item{purvars }{vector with indices for pure variables.}
#' \item{purityspec }{matrix with purity values for each resolved components.}
#' \item{stdspec }{matrix with weighted standard deviation values.}
#' \item{purity }{vector with purity values for resolved components.}
#' \item{wavelength }{vector with wavelength/wavebands/wavenumbers.}
#' 
#' @details  
#' The unmixing is done with two steps. First, pure variables are identified using methods, described
#' in [1-2]. Then the selected pure variables are used to resolve the pure spectra and contributions
#' with ordinary least squares.
#' 
#' By default method works with spectra, rows of the data matrix, but it is possible to make the
#' unmixing using concentration profiles (with parameter \code{by.spec} set to \code{FALSE}). Also 
#' an inverse second derivative of the data can be used for either unmixing step or for both 
#' identifying the pure variables and for unmixing in case if spectral peaks are overlapped.    
#'  
#' @seealso 
#' The class has numerous methods for exploring and refining the unmixing results.
#' \tabular{ll}{
#'  \code{print} \tab prints information about a \code{purity} object.\cr
#'  \code{summary} \tab shows summary statistics for the object.\cr
#'  \code{plot} \tab plot overview of unmixing results.\cr
#'  \code{\link{plotPurity}} \tab shows purity values and weighted standard deviation as a line plot.\cr
#'  \code{\link{plotSpectra}} \tab shows a plot with resolved spectra.\cr
#'  \code{\link{plotConcentrations}} \tab shows a plot with resolved concentrations.\cr
#'  \code{\link{plotVariance}} \tab shows a plot with residual variance vs. number of components.\cr
#'  \code{\link{plotSpectraComparison}} \tab compares resolved and ethalon spectra if available.\cr
#'  \code{\link{image.purity}} \tab if data is a hypercube, shows resolved concetration map as an image.\cr
#'  \code{\link{explore}} \tab a GUI tool to explore and tune the results (requires `shiny`).\cr
#' }
#' The steps of the algorithms can be also applied separately with functions 
#' \code{\link{getpurevars}} and \code{\link{unmix}}. 
#' 
#' @references 
#' 1. W. Windig, J. Guilment, Anal. Chem. 63 (1991) 1425-1432.
#' 2. W. Windig, Chemom. Intell. Lab. Syst. 23 (1994) 71-86. 
#'  
#' @examples  
#'
#' ## 1. Resolving of Raman spectra of carbonhydrates
#' 
#' # get the data
#' data(puredata)
#' 
#' # do unmixing with standard settings
#' res = purity(carbs$D, 3, wavelength = carbs$wavelength)
#' summary(res)
#' plot(res)
#' 
#' # compare the unmixed spectra with ethalon
#' par(mfrow = c(3, 1))
#' plotSpectraComparison(res, carbs$S, 1)
#' plotSpectraComparison(res, carbs$S, 2)
#' plotSpectraComparison(res, carbs$S, 3)
#' par(mfrow = c(1, 1))
#' 
#' # do unmixing with two different settings and compare purity for first component
#' res1 = purity(carbs$D, 3, offset = 5, use.deriv = FALSE)
#' res2 = purity(carbs$D, 3, offset = 15, use.deriv = TRUE, savgol = c(2, 5, 2))
#' par(mfrow = c(2, 1))
#' plotPurity(res1, 1)
#' plotPurity(res2, 1)
#' par(mfrow = c(1, 1))
#' 
#' ## 2. Hyperspectral image of oil-in-water emulsion
#' 
#' # get the data
#' data(puredata)
#' 
#' # it is important that variable has proper attributes
#' show(attr(hsi$D, 'width'))
#' show(attr(hsi$D, 'height'))
#' 
#' # do unmixing
#' res = purity(hsi$D, 3, offset = 3, wavelength = hsi$wavelength)
#' summary(res)
#' plot(res)
#' 
#' # show concentration maps
#' par(mfrow = c(1, 3))
#' image(res, 1)
#' image(res, 2)
#' image(res, 3)
#' 
#' 
#' @export 
purity = function(spectra, ncomp, offset = 5, by.spec = T, use.deriv = 0, savgol = c(2, 5, 2), 
                  exclude = NULL, wavelength = NULL)
{
  # transpose data if resolving is in contribution space
  if (!by.spec)
    spectra = t(spectra)
  
  # prepare data and take the derivative (if needed)
  D = spectra
  deriv = NULL
  if (use.deriv > 1)
  {  
    deriv = -savgol(D, savgol)
    deriv[deriv < 0] = 0
    D = deriv
  }

  # generate wavelength if not provided  
  if (is.null(wavelength))
  {
    wavelength = 1:ncol(spectra)
  }  
  
  # get pure variables and unmix data
  res1 = getpurevars(D, ncomp, offset = offset, exclude = exclude)
  res2 = unmix(spectra, D[, res1$purevars, drop = F], by.spec = by.spec)

  # combine all values to the final list
  res = c(res1, res2)
  res$by.spec = by.spec 
  res$use.deriv = use.deriv
  res$savgol = savgol
  res$wavelength = wavelength
  
  res$call = match.call()   
  class(res) = "purity"
  res
  
}

#' Identifies pure variables
#'
#' @description 
#' The method identifies indices of pure variables using the SIMPLISMA
#' algorithm.
#' 
#' @param D
#' matrix with the spectra or their inverted second derivative
#' @param ncomp
#' number of pure components 
#' @param  offset
#' offset in percent for calculation of parameter alpha
#' @param exclude      
#' which columns exclude when detecting the pure variables
#' 
#' @return 
#' The function returns a list with with following fields:
#' \item{ncomp }{number of pure components.}
#' \item{purvars }{vector with indices for pure variables.}
#' \item{purityspec }{matrix with purity values for each resolved components.}
#' \item{stdspec }{matrix with weighted standard deviation values.}
#' \item{purity }{vector with purity values for resolved components.}
#'
#' @export
getpurevars = function(D, ncomp, offset = 5, exclude = NULL)
{
  nspec = nrow(D)
  nvar = ncol(D)
  
  # get indices for excluded columns if provided
  colind = rep(TRUE, nvar)
  if (!is.null(exclude))
  {  
    if (max(exclude)> nvar || min(exclude) < 1)
      stop('Wrong values for the parameter "exclude"!')
    colind[exclude] = FALSE
  }
  
  
  # calculate purity spectrum
  mu = apply(D, 2, mean)
  sigma = apply(D, 2, sd) * sqrt((nspec - 1) / nspec)
  alpha = offset / 100 * max(mu)
  purity = sigma / (mu + alpha)
  
  # calculate variance-covariance matrix 
  Dprime = sweep(D, 2, (sqrt(mu^2 + (sigma + alpha)^2)), '/')
  R = 1/nrow(Dprime) * (t(Dprime) %*% Dprime)
  
  # loop to locate the pure variables variables
  purityspec = matrix(0, nrow = ncomp, ncol = nvar)
  stdspec = matrix(0, nrow = ncomp, ncol = nvar)
  purevars = rep(0, ncomp)
  purevals = rep(0, ncomp)
  for (i in 1:ncomp)
  {
    weights = matrix(0, ncol = ncol(D), nrow = 1)
    for (j in 1:ncol(D))
    {
      weights[j] = det(R[c(j, purevars[1:(i-1)]), c(j, purevars[1:(i-1)]), drop = F]);
    }  
    
    purityspec[i, colind] = purity[colind] * weights
    stdspec[i, colind] = sigma * weights
    purevars[i] = which.max(purityspec[i, ])
    purevals[i] = purityspec[i, purevars[i]]
  }  
  
  res = list()
  res$ncomp = ncomp
  res$purity = purevals
  res$purityspec = purityspec
  res$stdspec = stdspec
  res$purevars = purevars
  
  res
}

#' Unmix spectral data with least squares method
#'
#' @description
#' \code{unmix} decomposes the spectral data to a matrix with concentrations and 
#' a matrix with spectra of pure components using set of pure variables.
#' 
#' @param D
#' matrix with the spectral data
#' @param Dr
#' a matrix with concentration estimates (e.g. columns of D with pure variables)
#' @param by.spec
#' logical, should the algorithm works with spectra (rows) or with concentrations (columns)
#'  
#' @return 
#' Returns a list with two matrices, `conc` for concentrations and `spec` for spectra as well
#' as a vector with residual variance values
#'  
#' @export 
unmix = function(D, Dr, by.spec = T)
{
  
  # resolve spectra and concentrations
  St = t(D) %*% Dr %*% solve(t(Dr) %*% Dr)
  Ct = D %*% St %*% solve(t(St) %*% St)
  
  # scale
  #if (savgol[1] > 0)
  #  f = as.matrix(rowSums(deriv))
  #else
  #  f = as.matrix(rowSums(D))
  f = as.matrix(rowSums(D))
  a = solve(t(Ct) %*% Ct) %*% t(Ct) %*% f
  
  if (length(a) == 1)
    A = a
  else
    A = diag(as.vector(a))
  
  Ct = Ct %*% A
  St = St %*% A %*% solve(t(A) %*% A)

  # calculate residual variance
  ncomp = ncol(Ct)
  resvar = rep(0, ncomp)
  for (i in 1:ncomp)
  {  
    exp = Ct[, 1:i, drop = F] %*% t(St[, 1:i, drop = F])
    res = D - exp
    resvar[i] = sum(res^2) / sum(D^2)
  }
 
  # add names 
  names(resvar) = colnames(Ct) = colnames(St) = paste('C', 1:ncomp, sep = '')
  rownames(Ct) = rownames(D)
  rownames(St) = colnames(D)

  # switch spectra and concentrations if unmixing was in contribution space
  if (by.spec)
  {
    spec = St
    conc = Ct
  }
  else
  {
    spec = Ct
    conc = St
  }  
  
  # add attributes for hyperspectral image if any
  attr(conc, 'width') = attr(D, 'width')
  attr(conc, 'height') = attr(D, 'height')

  # return the results 
  res = list(
    conc = conc,
    spec = spec,
    resvar = resvar
  )
  
  res
}

#' Savytzky-Golay filter
#' 
#' @description
#' Applies Savytzky-Golay filter to the rows of data matrix
#' 
#' @param data
#' a matrix with data values
#' @param  params
#' vector with parameters: derivative order, width of filter and polynomial order
#' 
#' @export
savgol = function(data, params)
{
  nobj = nrow(data)
  nvar = ncol(data)
  pdata = matrix(0, ncol = nvar, nrow = nobj)
 
  deriv = params[1]
  width = params[2]
  porder = params[3]
  
  for (i in 1:nobj)
  {
    d = data[i, ]
    
    w = (width - 1)/2                        
    f  = pinv(outer(-w:w, 0:porder, FUN = "^"))  
    d = convolve(d, rev(f[deriv + 1, ]), type = "o")     
    pdata[i, ] = d[(w + 1) : (length(d) - w)] 
  }  
  
  pdata
}

#' Pseudo-inverse matrix
#' 
#' @description
#' Computes pseudo-inverse matrix using SVD
#' 
#' @param data
#' a matrix with data values to compute inverse for
#' 
#' @export
pinv = function(data)
{
  # Calculates pseudo-inverse of data matrix
  s = svd(data)
  s$v %*% diag(1/s$d) %*% t(s$u)
}

#' Image with concentration map
#' 
#' @description 
#' Make an image with map of resolved concentrations 
#' @param x
#' object with purity results
#' @param ncomp
#' which component make the map for
#' @param main
#' main title for the plot
#' @param col
#' a color map for the image
#' @param xaxt
#' x-axis type
#' @param yaxt
#' y-axis type
#' @param ...
#' other image parameters
#' 
#' @details 
#' The method will work only if original spectral data
#' has attributes \code{width} and \code{height} pointing out
#' that it is a hypercube
#' 
#' @export
image.purity = function(x, ncomp, main = NULL, col = getColors(256), xaxt = 'n', yaxt = 'n', ...)
{
  if (ncomp < 1 || ncomp > x$ncomp)
    stop('Wrong value for argument "ncomp"!')
 
  img = x$conc[, ncomp, drop = F]
  width = attr(x$conc, 'width')
  height = attr(x$conc, 'height')
  
  if (is.null(main))
    main = sprintf('Concentration map for %s', colnames(x$conc)[ncomp])
  
  if (is.null(width) || is.null(height))
    stop('The provided spectral data was not a hyperspectral image!')
  
  dim(img) = c(height, width)  
  image(img, col = col, main = main, xaxt = xaxt, yaxt = yaxt, ...)
}  


#' Overview plot for purity object
#' 
#' @description
#' Shows a set of plots () for .
#' 
#' @param x
#' a  (object of class \code{purity})
#' @param comp
#' for which component(s) make the plot for
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{purity}} function.
#' 
#' @export
plot.purity = function(x, comp = 1:x$ncomp, ...)
{   
  par(mfrow = c(2, 2))
  plotSpectra(x, comp = comp)
  
  if (is.null(attr(x$conc, 'width')))
    plotConcentrations(x, comp = comp)
  else
    image(x, ncomp = comp[1])
  plotPurity(x, ncomp = comp[1])
  plotVariance(x)
  par(mfrow = c(1, 1))
}

#' Print method for purity object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' unmixing results (object of class \code{purity})
#' @param ...
#' other arguments
#'
#' @export 
print.purity = function(x, ...)
{
  cat("\nUnmixing results (class purity)\n")
  cat('\nCall:\n')
  print(x$call)
  cat('\nMajor fields:\n')
  cat("$spec - resolved spectra of pure components\n")
  cat("$conc - resolved concentrations of pure components\n")
  cat("$ncomp - number of resolved components\n")
  cat("$purvars - vector with indices for pure variables\n")
  cat("$purityspec - matrix with purity values for each resolved components\n")
  cat("$stdspec - matrix with weighted standard deviation values\n")
  cat("$resvar - vector with residual variance for resolved components\n")
  cat("$purity - vector with purity values for resolved components\n")
  cat("$wavelength - vector with wavelength/wavebands/wavenumbers\n")
  cat('\nTry summary(obj) and plot(obj) to see the results.\n')   
}

#' Summary method for purity object
#' 
#' @description
#' Shows summary statistics for purity objects.
#' 
#' @param object
#' unmixing results (object of class \code{purity})
#' @param ...
#' other arguments
#' 
#' @export 
summary.purity = function(object, ...)
{
  cat('\nSummary statistics for unmixing results:\n')
  cat(sprintf('Number of resolved components: %d\n', object$ncomp))
  cat(sprintf('Total explained variance: %.1f%%\n', (1 - object$resvar[object$ncomp]) * 100))
  cat('\n')
  res = rbind(object$purevars, object$wavelength[object$purevars], 
              round(100 * diff(c(0, 1 - res$resvar)), 1),
              round(object$purity, 3))
  colnames(res) = paste('C', 1:object$ncomp, sep = '')
  rownames(res) = c('Indices', 'Wavelength', 'Expvar, %', 'Purity')
  show(res)
  
  cat('\n')
}
