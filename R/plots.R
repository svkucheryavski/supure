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

#' Purity spectra plot  
#' 
#' @description 
#' Shows a plot with purity spectra for each or user selected component
#' 
#' @param obj
#' object with unmixing results
#' @param ncomp
#' which component to make the plot for
#' @param origspec
#' original spectral data used for unmixing
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
                      main = sprintf('Purity for C%d', ncomp))
{
  if (class(obj) != 'purity')
    stop('This plot can be maid only for "purity" objects!')
  
  if (length(col) != 3)
    stop('Parameter "col" should have three values!')
  
  if (length(ncomp) != 1)
    stop('Parameter "ncomp" should have one value!')
  
  purityspec = obj$purityspec[ncomp, , drop = F]
  stdspec = obj$stdspec[ncomp, , drop = F]
  
  if (!obj$by.spec)
    wavelength = 1:nrow(obj$conc)
  else  
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
    if (!obj$by.spec)
      origspec = t(origspec)
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
#' Shows a plot with resolved spectra
#' 
#' @param obj
#' object with unmixing results
#' @param comp
#' which components make the plot for
#' @param origspec
#' original spectral data used for unmixing
#' @param show.labels
#' logical, show or not labels for the selected pure variables on the plot
#' @param show.legend
#' logical, show or not a legend for the plot
#' @param legend.str
#' a vector with legend items
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
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
#' @param ...
#' other maplot parameters
#' 
#' @export
plotSpectra = function(obj, comp = 1:obj$ncomp, origspec = NULL, show.labels = T,
                       show.legend = T, legend.str = NULL, ylim = NULL, xlim = NULL, 
                       col = getColors(length(comp)), lty = 1, lwd = 1, 
                       xlab = 'Wavenumbers', ylab = 'Intensity',  
                       main = 'Resolved spectra', ...)
{
  purespec = areanorm(t(obj$spec[, comp, drop = F]))
  
  purevars = obj$purevars
  wavelength = obj$wavelength
  
  if (!is.null(origspec))  
  {
    origspec = areanorm(origspec)
    data = rbind(origspec, purespec)
  } 
  else
  {
    data = rbind(purespec)
  }  
  
  if (is.null(xlim))
    xlim = c(min(wavelength), max(wavelength))
  
  if (is.null(ylim))
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
  
  if (show.legend)
  {
    if (is.null(legend.str))
      legend.str = colnames(obj$spec)[comp]
    legend(xlim[1] - dx/16, ylim[2] + dy/25, legend.str, col = col, lty = lty, 
           lwd = lwd, cex = 0.85, box.col = '#909090') 
    
  }  
  
  if (show.labels & !is.null(purevars))
  {  
    for (i in 1:length(comp))
    {
      text(wavelength[purevars[comp[i]]], purespec[i, purevars[comp[i]]], 
           sprintf('p%d', comp[i]), cex = 0.85, pos = 3)
    } 
  }  
} 

#' Concentrations plot  
#' 
#' @description 
#' Shows a plot with resolved concentration profiles  
#' 
#' @param obj
#' object with unmixing results
#' @param comp
#' which components make the plot for
#' @param origspec
#' original spectral data used for unmixing
#' @param col
#' color for the plot lines
#' @param lty
#' line type for plot lines
#' @param lwd
#' line width for plot lines
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param main
#' main title for the plot
#' @param ...
#' other parameters for matplot function
#' 
#' @export
plotConcentrations = function(obj, comp = 1:obj$ncomp, origspec = NULL, 
                              col = getColors(length(comp)), lty = 1, lwd = 1, 
                              xlab = 'Mixings', ylab = '',  
                              main = 'Resolved concentrations', ...)
{
  pureconc = areanorm(t(obj$conc[, comp, drop = F]))
  x = 1:nrow(obj$conc)
  
  
  if (!is.null(origspec))  
  {
    origconc = areanorm(t(origspec))
    data = rbind(origconc, pureconc)
  } 
  else
  {
    origconc = NULL
    data = pureconc
  }  
  
  xlim = c(min(x), max(x))
  ylim = c(min(data), max(data))
  dx = xlim[2] - xlim[1]
  dy = ylim[2] - ylim[1]
  
  if (!is.null(origconc))  
    matplot(x, t(origconc), type = 'l', col = '#F0F0F0', lty = lty, lwd = lwd,
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
  
  matlines(x, t(pureconc), lty = lty, lwd = lwd, col = col)  
  legend(xlim[1] - dx/16, ylim[2] + dy/25, colnames(obj$spec)[comp], col = col, lty = lty, 
         lwd = lwd, cex = 0.85, box.col = '#909090') 
} 

#' Residual variance plot for results of unmixing
#' 
#' @param obj
#' object with unmixing results
#' @param type
#' type of the plot
#' @param col
#' color for line and markers
#' @param pch
#' marker symbol
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param main
#' main title for the plot
#' @param ...
#' other parameters for plot function
#'  
#' @export
plotVariance = function(obj, type = 'b', col = getColors(1), pch = 16, 
                        xlab = 'Number of components', ylab = 'Variance, %', 
                        main = 'Residual variance', ...)
{
  plot(1:obj$ncomp, 100 * obj$resvar, type = type, col = col, pch = pch, xlab = xlab, ylab = ylab,
       main = main, xaxt = 'n', ...)
  axis(1, at = 1:obj$ncomp, labels = 1:obj$ncomp)
}

#' Compares resolved and original pure component spectra
#'
#' @description 
#' shows two spectra - a resolved and an original for selected pure component in the same plot
#' 
#' @param obj
#' object with unmixing results
#' @param purespec
#' matrix with pure spectra for the components (spectra in rows)
#' @param ncomp
#' for which component show the comparison for
#' @param col
#' color for the spectra (two values - one for each spectrum)
#' @param lwd
#' line width for the spectra (two values)
#' @param lty
#' line type for the spectra (two values)
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param main
#' main title for the plot
#' @param show.stat
#' logical, show or not relative residual variance between the spectra
#' @param ...
#' other plot parameters

#' @export
plotSpectraComparison = function(obj, purespec, ncomp = 1, col = getColors(2), lwd = c(2, 1), 
                                 lty = c(1, 1), xlab = 'Wavenumbers', ylab = 'Intensity', 
                                 main = NULL, show.stat = F, ...)
{
  rspectrum = t(obj$spec[, ncomp, drop = F])
  ospectrum = t(purespec[, ncomp, drop = F])
  
  ospectrum = areanorm(ospectrum)
  rspectrum = areanorm(rspectrum)
  
  res = sum((ospectrum - rspectrum)^2) / sum(ospectrum^2) 
  
  if (is.null(main))
  {
    if (show.stat)
      main = sprintf('Spectra comparison for C%d (%.3f%%)', ncomp, 100 * res)
    else  
      main = sprintf('Spectra comparison for C%d', ncomp)
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

