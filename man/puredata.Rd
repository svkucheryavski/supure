\name{puredata}
\alias{puredata}
\docType{data}
\title{
Datasets for spectral unmixing examples
}
\description{
\code{puredata} contains four different spectral datasets used in the chapter and also
in help examples of this package.
}

\usage{data(puredata)}
\format{
  There are four datasets: 
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

\details{

The first dataset, \code{carbs}, consists of Raman spectra of fructose, lactose and ribose generated 
from individual spectra with relative concentrations followed a simplex lattice design with four 
levels. The individual spectra were downloaded from publicly available SPECARB library [1]. It has 
clear pure variables for each component, e.g. 818 cm-1 for lactose, 357 cm-1 for fructose and 
542 cm-1 for ribose). 

The second dataset, \code{simdata}, consists of UV/Vis spectra of three polyaromatic carbohydrates 
and their mixtures. Values for concentration of the components in the mixtures were generated 
randomly, using uniform distribution with limits from 0 to 1 (Mole). The pure spectra can be resolved
using second derivative and following set of pure variables: 267 nm, 303 nm and 311 nm.

The third dataset, \code{mir} contains IR spectra of binary aqueous ethanol–glucose solutions. The 
spectra were acquired with a Matrix-MF FT-IR spectrophotometer with a diamond crystal probe produced 
by Art Photonics GmbH in a range 680…2900 cm-1 [2]. There are no pure spectra or concentrations 
available for this dataset.

The fourth dataset, \code{hsi}, is a 60x60 pixel Raman hyperspectral image of an oil-in-water 
emulsion [3]. The depicted area consists of an aqueous phase (with mostly glycerol and in a smaller 
concentration polawax), two droplet phases (one is polawax and another with parabens, polawax and 
some panthenol) and a structural phase with a crystalline polawax. The matrix with data, \code{D},
in this case has two additional attributes, \code{height} and \code{weight} for height and width
of the image in pixels, which allows to recognize this data by the methods implemented in the 
package as an image.

}

\references{
1. http://www.models.life.ku.dk/~specarb/specarb.html

2. A. Bogomolov, et al., Sensors Actuators B Chem. (2015). 

3. J.J. Andrew, M. a. Browne, I.E. Clark, T.M. Hancewicz, A.J. Millichope, Appl. Spectrosc. 52 (1998)
}

