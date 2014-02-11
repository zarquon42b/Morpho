Morpho
======
__Morpho__ provides a rich toolset for Geometric Morphometrics and mesh processing in R. This includes (among other stuff) mesh deformations based on reference points, permutation tests, detection of outliers, processing of sliding semi-landmarks, im- and export of a variety of triangular surface mesh files.


#### Installation of the R-package Morpho from CRAN: ####

Within R:
       
       install.packages("Morpho")


#### Installation of the R-package Morpho from sourceforge/github (latest release): ####
1. Make sure to work with the latest version of R and install dependencies (type the following commands into your R terminal): 
     
            
        install.packages(c("rgl", "MASS","doParallel","colorRamps","yaImpute","RcppArmadillo"))


  Also required are the packages 'Matrix' and 'parallel' which usually are already installed as R's recommended packages.


2. Download the version suitable for your OS from [sourceforge](https://sourceforge.net/projects/morpho-rpackage/) or [github](https://github.com/zarquon42b/Morpho/releases). Either the compiled package (for Windows and OS X) or the source tarball (Linux).

3. Installation command from within R: 
   
        install.packages("Path_to_downloaded_package_Morpho[Version_OS]",repos=NULL)

4. check if the package can be loaded:
        
        load package: library(Morpho)

#### Installation of the R-package Morpho (development snapshot) using *devtools*: ####
##### install prerequisites #####

1. install *devtools* from within R (Ubuntu/Debian users will have to install *libcurl4-gnutls-dev* beforehand):

        install.packages("devtools")

    * **Make sure to have the latest versions of Rcpp and RcppArmadillo installed!!**

2. Install build environment
    * **Windows:** Install latest version of *[Rtools](http://cran.r-project.org/bin/windows/Rtools)*
During installation of *Rtools* make sure to install the *toolchain*, and to select *"Edit the system path"* (and confirming the installers suggestions).
    * **OSX:** Install *[XCODE](https://developer.apple.com/xcode/)*

##### install Morpho #####

Run the following command in R:
        
        require(devtools)
        install_github("zarquon42b/Morpho", local=FALSE")


