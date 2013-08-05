\name{NEWS}
\title{News for Package 'Morpho'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in Morpho version 0.26-00 (2013-NN-NN)}{
  \itemize{
    \item white space formatting in code
    \item added option "mshape" to function 'orp' to specify a custom pole for projection
    \item made function 'fx', 'PCAcrova' and 'calcTang_U' interntal by
    prefixing fucntion name with '.'
    \item removed unnecessary sinus calculation in 'angle.calc'
    \item added docu for function 'barycenter'
    \item added a lot of examples to docu and documented data(nose)
    \item added documentation for 'angle.calc'
    \item added example in 'relaxLM'
    \item minor change in DESCRIPTION
  }
}


\section{Changes in Morpho version 0.25-1 (2013-08-02)}{
  \itemize{
    \item added docu for new data set "boneData"
    \item added examples in 'procSym' and replaced abind with 'bindArr' in 'fixLMtps' and 'procSym'
    \item added additional example data (landmarks and a bone surface mesh)
    \item added example in tps3d.Rd
    \item updated DESCRIPTION and README
    \item updated README added sourceforge and github download links
    \item moved Distance calculation in 'CVA' to '.CVAdists' and radically cleaned code; removal of redundant variables in 'groupPCA'
    \item added README.Rd
    \item added option reflect to 'procSym' (and underlying functions) and added copyright note for code copied from shapes package
    \item added NEWS.Rd
  }
}

\section{Changes in Morpho version 0.25 (2013-07-24)}{
  \itemize{
    \item updated COPYRIGHT file
    \item simplified code in 'CVA' and 'CVA.crova' and finally mentioned Marta Rufino for providing all examples in docu for 'CVA' (sorry Marta, I simply forgot)
    \item fixed issue when n-groups > n variables in 'CVA' and 'CVA.crova'
    \item replaced 'nn2' (RANN) with 'ann' (yaImpute) because RANN got (temporarily) banned from CRAN because of violation of policy
    \item added 'readLandmarks.csv' and docu to read landmarks stored in spreadsheet text files (e.g. comma separated)
    \item added check for empty groups in 'typprobClass' and 'covW'
    \item fixed some docu and added examples in 'pcaplot3d' and 'rotonmat'
    \item changed cube creation in 'deformGrid3d' (much faster now).
    \item added options weights and centerweight in 'rotonmat'
    \item added confirmation for plotting spheres in 'deformGrid3d' when n > 1000 (can be slow)
    \item 'pcaplot3d' now works if data has only 1 PC
    \item removed useless code remains 
    \item added docu for 'anonymize', 'r2morphoj', 'r2morphologika', 'NNshapeReg', 'kendalldist'
    \item added 'bindArr' to concatenate matrices/arrays
    \item switched from depends to import for colorRamps, MASS, RANN
    \item changed default core value in 'closemeshKD' from detectCores() to 1.
    \item added option 'keep' in 'rmVertex' 
    \item changed default method in 'meshDist.mesh3d' from vcglib to morpho
    \item renamed 'regdist.raw' to 'regdist' and removed original 'regdist' (too prone to user error)
  }
}
