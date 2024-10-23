
Morpho [![Unit Tests](https://github.com/zarquon42b/Morpho/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zarquon42b/Morpho/actions/workflows/R-CMD-check.yaml)  [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)  [![Downloads](https://cranlogs.r-pkg.org/badges/last-month/Morpho?color=brightgreen)](https://cranlogs.r-pkg.org/badges/last-month/Morpho)
======
__Morpho__ provides a rich toolset for Geometric Morphometrics and mesh processing in R. This includes (among other stuff) mesh deformations based on reference points, permutation tests, detection of outliers, processing of sliding semi-landmarks, im- and export of a variety of triangular surface mesh files.


#### Installation of the R-package Morpho from CRAN: ####

Prerequisites for Debian/Ubuntu (needed to build dependency package [*sf*](https://cran.r-project.org/package=sf))

   `sudo apt install libudunits2-dev libgeos-dev`
   

Prerequisites for Fedora (needed to build dependency package [*sf*](https://cran.r-project.org/package=sf))

   `dnf install udunits2-devel geos-devel gdal-devel proj-devel`
   

Within R:
       
	install.packages("Morpho")


#### Installation of the R-package Morpho (development snapshot) using *devtools*: ####

##### Install prerequisites #####

1. Install *devtools* from within R (Ubuntu/Debian users will have to install *libcurl4-gnutls-dev* beforehand):

		install.packages("devtools")

	**Make sure to have the latest versions of Rcpp and RcppArmadillo installed!!**
	   

2. Install build environment
    * **Windows:** Install latest version of *[Rtools](http://cran.r-project.org/bin/windows/Rtools)*
During installation of *Rtools* make sure to install the *toolchain*, and to select *"Edit the system path"* (and confirming the installers suggestions).
    * **OSX:** Install *[XCODE](https://developer.apple.com/xcode/)*
	
3. Install prerequisites for Debian/Ubuntu (needed to build dependency package  [*sf*](https://cran.r-project.org/package=sf))

  `sudo apt install libudunits2-dev libgeos-dev`


##### Install Morpho #####

Run the following command in R:
        
	require(devtools)
	install_github("zarquon42b/Morpho")


