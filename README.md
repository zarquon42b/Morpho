#### Installation of the R-package Morpho (latest release): ####
   0. Make sure to work with the latest version of R and install dependencies (type the following commands into your R terminal): 
     
            
        install.packages(c("rgl", "MASS","doParallel","colorRamps","yaImpute"))


Also required are the packages 'Matrix' and 'parallel' which usually are already installed as R's recommended packages.


* Download the version suitable for your OS from [sourceforge](https://sourceforge.net/projects/morpho-rpackage/) or [github](https://github.com/zarquon42b/Morpho/releases). Either the compiled package (for Windows and OS X) or the source tarball (Linux).

* Installation command from within R: 
   
        install.packages("Path_to_downloaded_package_Morpho[Version_OS]",repos=NULL)

* check if the package can be loaded:
        
        load package: library(Morpho)

#### Installation of the R-package Morpho (development snapshot) using *devtools*: ####
##### install prerequisites ####

1. install *devtools* from within R (Ubuntu/Debian users will have to install *libcurl4-gnutls-dev* beforehand):

        install.packages("devtools")

2. Install build environment
    * **Windows:** Install latest version of *[Rtools](http://cran.r-project.org/bin/windows/Rtools)*
During installation of *Rtools* make sure to install the *toolchain*, and to select *"Edit the system path"* (and confirming the installers suggestions).
    * **OSX:** Install *[XCODE](https://developer.apple.com/xcode/)*


3. Run the following command in R:
        
        require(devtools)
        install_url("https://github.com/zarquon42b/Rvcg/archive/master.zip")



#### Installation of the command line tools (strongly recommended):###
   To  gain full functionality of the Morpho R-package (using sliding landmarks, importing meshfiles other than ascii ply files, etc), it is required to download and install the latest version of the command line programms.

1. Download the binaries appropriate for your OS from [https://github.com/zarquon42b/trimesh-cxx/releases](https://github.com/zarquon42b/trimesh-cxx/releases)

* Install files:

    * **Windows:** simply double click trimesh-tools.msi and follow installer instructions.
    * **OSX:**
        1. extract the .tgz archive and copy these files e.g. to /usr/bin/ 
        * you also need the QT libraries QtGui and QtCore installed. If you don't know how to do that, you can download the complete QT libraries from http://qt.nokia.com/downloads and install them on your machine.
    * **Linux:** 
        * Debian/Ubuntu: Please use my [PPA](https://launchpad.net/~zarquon42/+archive/ppa). The package is called *trimesh-tools*.
        * All other systems: Compile the binaries yourself ([see below](\#compilation-of-command-line-tools)).

* Test if the system finds the files (necessary for Morpho).
    1. open a command line terminal
    * type *ply2ascii*
    * if everything is alright, you see a help message.
	
   
#### Compilation of command line tools:   

  1. Install QT-SDK (http://qt-project.org/downloads)
  * make sure to have C++ compilers installed
  * Get the latest sources
       1. Download the [tarball](https://github.com/zarquon42b/trimesh-cxx/archive/0.2.5.tar.gz)
       * or use git to obtain the latest snapshot
                    
            git clone git://git.code.sf.net/p/morpho-rpackage/trimesh-cxx trimesh-cxx
                
	

######Compilation

        cd trimesh-cxx
        qmake 
        make

The compiled binaries are found in the folder trimesh-cxx/bin which is created during the compilation process.

Copy these binaries to any convenient directory and make sure that the directory is added to the system's PATH variable.