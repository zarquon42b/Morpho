\name{NEWS}
\title{News for Package 'Morpho'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in Morpho version 0.25 (2013-07-24)}{
  \itemize{
      \item updated COPYRIGHT file
      \item simplified code in CVA and CVA.crova and finally mentioned Marta Rufino for providing all examples in docu for CVA (sorry Marta, I simply forgot)
      \item fixed issue when n-groups > n variables in CVA and CVA.crova
      \item    replaced nn2 (RANN) with ann (yaImpute) because RANN got (temporarily) banned from CRAN because of violation of policy
      \item    added readLandmarks.csv and docu to read landmarks stored in spreadsheet text files (e.g. comma separated)
      \item    added check for empty groups in typprobClass and covW
      \item    fixed some docu and added examples in pcaplot3d and rotonmat
      \item    changed cube creation in deformGrid3d (much faster now).
      \item    added options weights and centerweight in rotonmat
      \item    added confirmation for plotting spheres in deformGrid when n > 1000 (can be slow)
      \item    pcaplot3d now works if data has only 1 PC
      \item    removed useless code remains 
      \item    added docu for anonymize, r2morphoj, r2morphologika, NNshapeReg, kendalldist
      \item    added bindArr to concatenate matrices/arrays
      \item    switched from depends to import for colorRamps, MASS, RANN
      \item    changed default core value in closemeshKD from detectCores() to 1.
      \item    added option 'keep' in rmVertex 
      \item    changed default method in meshDist.mesh3d from vcglib to morpho
      \item    renamed regdist.raw to regdist and removed original regdist (too prone to user error)
    }
  }
    