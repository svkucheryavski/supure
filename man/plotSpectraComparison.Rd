% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plots.R
\name{plotSpectraComparison}
\alias{plotSpectraComparison}
\title{Compares resolved and original pure component spectra}
\usage{
plotSpectraComparison(obj, purespec, ncomp = 1, col = getColors(2),
  lwd = c(2, 1), lty = c(1, 1), xlab = "Wavenumbers",
  ylab = "Intensity", main = NULL, show.stat = F, ...)
}
\arguments{
\item{obj}{object with unmixing results}

\item{purespec}{matrix with pure spectra for the components (spectra in rows)}

\item{ncomp}{for which component show the comparison for}

\item{col}{color for the spectra (two values - one for each spectrum)}

\item{lwd}{line width for the spectra (two values)}

\item{lty}{line type for the spectra (two values)}

\item{xlab}{label for x axis}

\item{ylab}{label for y axis}

\item{main}{main title for the plot}

\item{show.stat}{logical, show or not relative residual variance between the spectra}

\item{...}{other plot parameters}
}
\description{
shows two spectra - a resolved and an original for selected pure component in the same plot
}

