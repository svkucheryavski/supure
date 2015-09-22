\name{supure}
\alias{supure}

\title{
Package for spectral unmixing with pure variables approach 
}

\description{
\code{supure} is an R package implementing methods for spectral unmixing (curve resolution) based on
the pure variables approach. 
}

\details{

\code{supure} is an R package implementing methods for spectral unmixing (curve resolution) based on
the pure variables approach. It includes variants of SIMPLSIMA [1,2] and a simple implementation of 
MCR-ALS with angle conctraint [3,4]. The package supplements a chapter in a book "Resolving spectral
mixtures, Data handling in Science and Technology", edited by Cyril Ruckebusch which will be published 
by Elsevier in 2016. It also contains all datasets used in the chapter and a vignette with some of 
the examples used there as well.

There are two main classes: \code{\link{purity}} for unmixing pure variable approach and its 
variations and \code{\link{mcrals}} for the angle constrained MCR-ALS. Each class has several methods 
for visualisation and tuning the results of unmixing.

There are four datasets, which can be loaded using command \code{data(puredata)}: 
\tabular{ll}{
  \code{carbs} \tab generated Raman spectra of ribose, fructose and lactose mixtures.\cr
  \code{simdata} \tab generated UV/Vis spectra of polyaromatic hydrocarbons mixtures.\cr
  \code{mir} \tab Mid infrared spectra of ethanol and glucose.\cr
  \code{hsi} \tab Raman hyperspectral image of oil-in-water emulsion.\cr
}

Each dataset is a list with following fields:
\tabular{ll}{
  \code{D} \tab data matrix with spectra of mixtures.\cr
  \code{S} \tab optional, matrix with spectra of pure components.\cr
  \code{С} \tab optional, matrix with concentrations of pure components.\cr
  \code{wavelength} \tab wavelength or wavenumbers.\cr
}

}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\references{
1. W. Windig, J. Guilment, Anal. Chem. 63 (1991) 1425–1432.

2. W. Windig, Chemom. Intell. Lab. Syst. 23 (1994) 71–86. 

3. W. Windig, M.R. Keenan, Appl. Spectrosc. 65 (2011) 349–357.
}