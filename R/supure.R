#' Makes unmixing of set of spectra to pure components.
#'
#' @param data
#' a matrix with spectral values
#' @param  ncomp
#' number of pure components to identify (optional)
#' @param  suggest
#' a sequence with variable numbers suggested as pure variables for each component
#' @param  savgol.deriv
#' a derivative order (if 0 no derivative will be used)
#' @param  savgol.width
#' a width of filter used to smooth signal prior to take a derivative
#' @param  savgol.porder
#' a polynomial order used for smoothing
#' @param do.norm
#' logical, normalize the spectral data to unit area prior to unmixing or not
#' @param  wavelength
#' a vector with wavelength if necessary
#' 
#' @return
#' The method returns an object of \code{supure} class with following fields.
#' \code{supure } object - list with several fields, including:
#' \item{spec }{resolved spectra of pure components.}
#' \item{conc }{resolved concentrations of pure components.}
#' \item{purvars }{vector with indices for pure variables.}
#' \item{purity }{vector with values for purity function.}
#'  
#'  
#' @export 
supure = function(spectra, ncomp, suggest = NULL, savgol.deriv = 2, 
                  savgol.width = NULL, savgol.porder = 2, do.norm = TRUE,
                  wavelength = NULL)
{
  ## find pure variables 
  purity = NULL
  purevars = suggest
  
  # prepare data and take the derivative (if needed)
  D = spectra
  
  # loop to locate the variables
  # ... !!! ...

  ## resolve the pure spectra and concentrations
  
  # get the initial estimates of concentrations
  if (savgol.deriv > 0)
  {  
    deriv = -savgol(D, savgol.width, savgol.porder, savgol.deriv)
    deriv[deriv < 0] = 0
    Dr = deriv[, purevars, drop = F]  
  }
  else
  {  
    Dr = D[, purevars, drop = F]  
  }
  
  # resolve spectra and concentrations
  St = t(D) %*% Dr %*% solve(t(Dr) %*% Dr)
  Ct = D %*% St %*% solve(t(St) %*% St)
  
  # scale
  if (savgol.deriv > 0)
    f = as.matrix(rowSums(deriv))
  else
    f = as.matrix(rowSums(D))
 
  a = solve(t(Ct) %*% Ct) %*% t(Ct) %*% f
  A = diag(as.vector(a))
  
  Ct = Ct %*% A
  St = St %*% A %*% solve(t(A) %*% A)
  
  ## return the results
  res = list()
  
  res$spec = St
  res$conc = Ct
  res$purity = purity
  res$purevars = purevars
  res$wavelength = wavelength
  
  res$call = match.call()   
  class(res) = "supure"
  res
}

#' Savytzky-Golay filter
#' 
#' @description
#' Applies Savytzky-Golay filter to the rows of data matrix
#' 
#' @param data
#' a matrix with data values
#' @param width
#' width of the filter window
#' @param porder
#' order of polynomial used for smoothing
#' @param dorder
#' order of derivative to take (0 - no derivative)
#' 
#' @export
savgol = function(data, width = 3, porder = 1, dorder = 0)
{
  nobj = nrow(data)
  nvar = ncol(data)
  
  pdata = matrix(0, ncol = nvar, nrow = nobj)
  
  for (i in 1:nobj)
  {
    d = data[i, ]
    
    w = (width - 1)/2                        
    f  = pinv(outer(-w:w, 0:porder, FUN = "^"))  
    d = convolve(d, rev(f[dorder + 1, ]), type = "o")     
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

#' Overview plot for supure object
#' 
#' @description
#' Shows a set of plots () for .
#' 
#' @param x
#' a  (object of class \code{pca})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{supure}} function.
#' 
#' @export
plot.supure = function(x, ...)
{   
}

#' Print method for SUPURE object
#' 
#' @method print supure
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' unmixing results (object of class \code{supure})
#' @param ...
#' other arguments
#'
#' @export 
print.supure = function(x, ...)
{
}

#' Summary method for SUPURE object
#' 
#' @description
#' Shows some statistics () for the .
#' 
#' @param object
#' unmixing results (object of class \code{supure})
#' @param ...
#' other arguments
#' 
#' @export 
summary.supure = function(object, ...)
{
}
