#' Makes unmixing of set of spectra to pure components.
#'
#' @param data
#' a matrix with spectral values
#' @param  ncomp
#' number of pure components to identify (optional)
#' @param use.deriv
#' shall a second derivative be used (0 - no, 1 - only for unmixing, 2 - both for locating pure 
#' variables and for unmixing) 
#' @param  savgol
#' parameters for Savitzky-Goley filter: derivative order, width of filter and polynomial order
#' @param exclude
#' vector of column numbers that should be excluded when pure variables are being detected
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
supure = function(spectra, ncomp, offset = 5, use.deriv = 0, savgol = c(2, 5, 2), exclude = NULL, 
                  wavelength = NULL)
{
  ## 1. find pure variables
 
  nspec = nrow(spectra)
  nvar = ncol(spectra)
  colind = rep(TRUE, nvar)
  
  if (!is.null(exclude))
  {  
    if (max(exclude)> nvar || min(exclude) < 1)
      stop('Wrong values for the parameter "exclude"!')
    colind[exclude] = FALSE
  }
  
  # prepare data and take the derivative (if needed)
  deriv = NULL
  if (use.deriv > 1)
  {  
    deriv = -savgol(D, savgol)
    deriv[deriv < 0] = 0
    D = deriv
  }
  
  D = D[, colind]
  
  if (is.null(wavelength))
  {
    wavelength = 1:nvar
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
  }  
  
  ## 2. resolve the pure spectra and concentrations and return results
  if (use.deriv > 0)
    res = unmix(spectra, purevars, savgol)
  else
    res = unmix(spectra, purevars)
  
  res$use.deriv = use.deriv
  res$savgol = savgol
  res$purity = purity
  res$purityspec = purityspec
  res$stdspec = stdspec
  res$purevars = purevars
  res$wavelength = wavelength
  res$ncomp = length(purevars)
  res$call = match.call()   
  class(res) = "supure"
  res
  
}

#' Unmix spectral data with least squares method
#'
#' @description
#' `unmix` decomposes the spectral data to a matrix with concentrations and 
#' a matrix with spectra of pure components using set of pure variables.
#' 
#' @param D
#' matrix with the spectral data
#' @param purevars
#' a vector with pure variables
#' @param  savgol
#' parameters for Savitzky-Goley filter: derivative order, width of filter and polynomial order
#'  
#' @return 
#' Returns a list with two matrices, `conc` for concentrations and `spec` for spectra
#'  
#' @export 
unmix = function(D, purevars, savgol = c(0, 1, 1))
{
  
  # get the initial estimates of concentrations
  if (savgol[1] > 0)
  {  
    deriv = -savgol(D, savgol)
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
  if (savgol[1] > 0)
    f = as.matrix(rowSums(deriv))
  else
    f = as.matrix(rowSums(D))
  
  a = solve(t(Ct) %*% Ct) %*% t(Ct) %*% f
  
  if (length(a) == 1)
    A = a
  else
    A = diag(as.vector(a))
  
  Ct = Ct %*% A
  St = St %*% A %*% solve(t(A) %*% A)
 
  colnames(Ct) = colnames(St) = paste('C', 1:ncol(Ct), sep = '')
  rownames(Ct) = rownames(D)
  rownames(St) = colnames(D)
 
  ncomp = ncol(Ct)
  
  resvar = rep(0, ncomp)
  for (i in 1:ncomp)
  {  
    exp = Ct[, 1:i, drop = F] %*% t(St[, 1:i, drop = F])
    res = D - exp
    resvar[i] = sum(res^2) / sum(D^2)
  }
  
  res = list(
    conc = Ct,
    spec = St,
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
#' @param width
#' width of the filter window
#' @param porder
#' order of polynomial used for smoothing
#' @param dorder
#' order of derivative to take (0 - no derivative)
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

#' Normalization of spectra
#' 
#' @description 
#' normalize spectra to unit area
#' 
#' @param spectra
#' matrix of vector with spectral data
#' 
#' @export
areanorm = function(spectra)
{
  spectra = as.matrix(spectra)
  w = rowSums(abs(spectra))
  normspectra = sweep(spectra, 1, w, '/')
  normspectra
}  

#' Returns colors from colorbrewer2 palette
#' 
#' @param ncols
#' number of colors
#' 
#' @export
getColors = function(ncols = 1)
{   
  palette = c("#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", 
              "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F")
  
  if (ncols == 3)
  {
    colvals = c("#3288BD", "#91CF60", "#D53E4F")
  } 
  else
  {
    colfunc = colorRampPalette(palette, space = 'Lab')      
    colvals = colfunc(ncols)      
  }  
  
  colvals
}

#' Purity spectra plot  
#' 
#' @description 
#' Shows a plot with purity spectra for each or user selected component
#' 
#' @param obj
#' object with supure results
#' @param ncomp
#' which component to make the plot for
#' @param col
#' colors for the std, purity and original spectra
#' @param lty
#' line type for the component spectra
#' @param lwd
#' line width for the component spectra
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param main
#' main title for the plot
#' 
#' @export
plotPurity = function(obj, ncomp = 1, origspec = NULL,
                       col = c('#62B7EF', '#D84E5E', '#F0F0F0'), lty = 1, lwd = 1, 
                       xlab = 'Wavenumbers', ylab = c('Weighted std', 'Purity'),  
                       main = sprintf('Purity of component #%d', ncomp), ...)
{
  if (length(col) != 3)
    stop('Parameter "col" should have three values!')

  if (length(ncomp) != 1)
    stop('Parameter "ncomp" should have one value!')
  
  purityspec = obj$purityspec[ncomp, , drop = F]
  stdspec = obj$stdspec[ncomp, , drop = F]
  wavelength = obj$wavelength
  
  xlim = c(min(wavelength), max(wavelength))
  ylim = c(min(stdspec), max(stdspec))
  dx = xlim[2] - xlim[1]
  dy = ylim[2] - ylim[1]
  
  
  # correct margin to have two y-axis
  par(mar = c(5, 4, 4, 5) + .1)
  # first plot with std and spectra  
  plot(0,
       type = 'n',
       ylab = ylab[1],
       xlab = xlab,
       main = main,
       xlim = c(xlim[1], xlim[2] + dx/20),
       ylim = c(ylim[1], ylim[2] + dy/20)
  )   
  
  # if original data is provided scale and show it
  if (!is.null(origspec))
  {
    origspec = (origspec - min(origspec))/(max(origspec) - min(origspec))
    origspec = ylim[1] + origspec * (ylim[2] - ylim[1])
    matlines(wavelength, t(origspec), type = 'l', col = col[3], lty = lty, lwd = lwd)  
  }  
  
  # show standard deviation plot
  lines(wavelength, t(stdspec), lty = lty, lwd = lwd, col = col[1])
  
  # show legend
  l = legend(0, 0, c('Std', 'Purity'), col = col, cex = 0.80, lty = lty, lwd = lwd, plot = F)
  legend(xlim[2] - l$rect$w + dx/20, ylim[2] + dy/20,  c('Std', 'Purity'), col = col, 
         lty = lty, lwd = lwd, cex = 0.80, box.col = '#909090', bg = '#ffffff') 
  
  # second plot with purity spectrum
  par(new = TRUE)
  matplot(wavelength, t(purityspec), type = 'l', lty = lty, lwd = lwd, col = col[2],
          xlim = c(xlim[1], xlim[2] + dx/20),
          xaxt = 'n', 
          yaxt = 'n', xlab = '', ylab = '')
  axis(4)
  mtext(ylab[2], side = 4, line = 3, cex = 0.7)
  
  points(wavelength[obj$purevars[ncomp]], purityspec[obj$purevars[ncomp]], pch = 16, col = col[2])
}  
  

#' Spectral plot  
#' 
#' @description 
#' Shows a plot with original spectral data, evaluated spectra of pure components and selected 
#' pure variables
#' 
#' @param obj
#' object with supure results
#' @param comp
#' which components make the plot for
#' @param show.labels
#' logical, show or not labels for the selected pure variables on the plot
#' @param col
#' color for the component spectra
#' @param lty
#' line type for the component spectra
#' @param lwd
#' line width for the component spectra
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param main
#' main title for the plot
#' 
#' @export
plotSpectra = function(obj, comp = 1:obj$ncomp, origspec = NULL, show.labels = T,
                       col = getColors(length(comp)), lty = 1, lwd = 1, 
                       xlab = 'Wavenumbers', ylab = 'Intensity',  
                       main = 'Resolved spectra', ...)
{
  purespec = areanorm(t(obj$spec[, comp, drop = F]))
  
  purevars = obj$purevars
  wavelength = obj$wavelength

  if (!is.null(origspec))  
  {
    origspec = areanorm(obj$data)
    data = rbind(origspec, purespec)
  } 
  else
  {
    data = rbind(purespec)
  }  
  
  xlim = c(min(wavelength), max(wavelength))
  ylim = c(min(data), max(data))
  dx = xlim[2] - xlim[1]
  dy = ylim[2] - ylim[1]
  
  if (!is.null(origspec))  
    matplot(wavelength, t(origspec), type = 'l', col = '#F0F0F0', lty = lty, lwd = lwd,
          ylab = ylab,
          xlab = xlab,
          main = main,
          xlim = c(xlim[1] - dx/20, xlim[2] ),
          ylim = c(ylim[1], ylim[2] + dy/20),
          ...)  
  else
    matplot(0, 0, type = 'n', col = '#F0F0F0', lty = lty, lwd = lwd,
            ylab = ylab,
            xlab = xlab,
            main = main,
            xlim = c(xlim[1] - dx/20, xlim[2] ),
            ylim = c(ylim[1], ylim[2] + dy/20),
            ...)  
    
  matlines(wavelength, t(purespec), lty = lty, lwd = lwd, col = col)  
  legend(xlim[1] - dx/16, ylim[2] + dy/25, colnames(obj$spec)[comp], col = col, lty = lty, 
         lwd = lwd, cex = 0.85, box.col = '#909090') 
  
  if (show.labels)
  {  
    for (i in 1:length(comp))
    {
      text(wavelength[purevars[comp[i]]], purespec[i, purevars[comp[i]]], 
           sprintf('p%d', comp[i]), cex = 0.85, pos = 3)
    } 
  }  
} 

#' Residuals plot for supure object
#' 
#' @param obj
#' object with supure results
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param main
#' main title for the plot
#' 
#' @export
plotResiduals = function(obj, type = 'b', col = getColors(1), pch = 16, 
                         xlab = 'Number of components', ylab = 'Variance, %', 
                         main = 'Residual variance', ...)
{
  plot(1:obj$ncomp, 100 * obj$resvar, type = type, col = col, pch = pch, xlab = xlab, ylab = ylab,
       main = main, xaxt = 'n', ...)
  axis(1, at = 1:obj$ncomp, labels = 1:obj$ncomp)
}

#' Compares resolved and original pure component spectra
#'
#' @export
plotSpectraComparison = function(obj, ospectra, ncomp = 1, col = getColors(2), lwd = c(2, 1), 
                                 lty = c(1, 1), xlab = 'Wavenumbers', ylab = 'Intensity', 
                                 main = NULL, show.stat = F, ...)
{
  rspectrum = t(obj$spec[, ncomp, drop = F])
  ospectrum = t(ospectra[, ncomp, drop = F])
   
  ospectrum = areanorm(ospectrum)
  rspectrum = areanorm(rspectrum)
  
  res = sum((ospectrum - rspectrum)^2) / sum(ospectrum^2) 
  
  if (is.null(main))
  {
    if (show.stat)
      main = sprintf('Spectra comparison for component #%d (%.3f%%)', ncomp, 100 * res)
    else  
      main = sprintf('Spectra comparison for component #%d', ncomp)
  }  

  xlim = c(min(obj$wavelength), max(obj$wavelength))
  ylim = c(min(c(ospectrum, rspectrum)), max(c(ospectrum, rspectrum)))
  dx = xlim[2] - xlim[1]
  dy = ylim[2] - ylim[1]
  
  plot(obj$wavelength, ospectrum, type = 'l', main = main, xlab = xlab, ylab = ylab, 
       xlim = c(xlim[1] - dx/20, xlim[2]), ylim = c(ylim[1], ylim[2] + dy/20),
       col = col[1], lty = lty[1], lwd = lwd[1], ...)
  lines(obj$wavelength, rspectrum, col = col[2], lty = lty[2], lwd = lwd[2])
  
  abline(v = obj$wavelength[obj$purevars[ncomp]], lty = 2, lwd = 2, col = '#909090') 
  legend(xlim[1] - dx/15, ylim[2] + dy/25, c('Original', 'Resolved'), col = col, lty = lty, 
         lwd = lwd, cex = 0.80, box.col = '#909090')
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
