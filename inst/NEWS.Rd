\name{NEWS}
\title{News for Package 'Morpho'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in Morpho version 0.26-00 (2013-NN-NN)}{
  \subsection{R-code}{
    \itemize{
      \item added function 'checkLM'. Useful to visualize if projection by
  'placePatch' was successful by browsing through sample. 
      \item added functions 'placePatch', 'createAtlas' and 'plotAtlas'
  for defining templates and functions to project (semi-) landmarks from
  an atlas onto all specimen in a given sample.
      \item fixed missing '-1' subtraction in 'covW'
      \item added function 'mergeMeshes' to concatenate multiple triangular meshes.
      \item made function 'procAOVsym' (procrustes Anova for object
      symmetry) public.
      \item made functions 'cSize', 'mcNNindex', 'meshcube', 'meshres', 'ProcGPA' public.
      \item added option 'info' in 'read.mpp'
      \item merged 'CreateL' and 'CreateL2D'
      \item rewrote large parts of 'CVA', 'groupPCA' and 'permudist' to
  make the code more efficient.
      \item 'covPCA' now allows permutation tests to check if pairwise
      distances between covariance matrices are significant and reports rownames for two-group case.
      \item removed unnecessary complicated check for empty groups by simply using function 'factor' in 'covW' and 'typprobClass'
      \item white space formatting in code
      \item added option "mshape" to function 'orp' to specify a custom pole for projection
      \item made function 'fx', 'PCAcrova' and 'calcTang_U'
  interntal-only by prefixing function name with '.'
      \item removed unnecessary sinus calculation in 'angle.calc'
    }
  }
  \subsection{Package documentation}{
    \itemize{
      \item added docu and examples for all newly exposed functions (see above)
      \item added examples for 'find.outliers'
      \item added more examples for covPCA, performing permutation tests
      \item added docu for function 'barycenter'
      \item added a lot of examples to docu and documented data(nose)
      \item added documentation for 'angle.calc'
      \item added example in 'relaxLM'
      \item added section "details" in procSym, added example in slider3d
      \item added example for procAOVsym
    }
  }
  \subsection{Misc}{
    \itemize{
      \item changed foreach, parallel and rgl (needs rgl version >=
      0.93.963) from 'Depends' to 'Imports' (affecting DESCRIPTION and NAMESPACE)
      \item massive NAMESPACE editing by switching to manual export instead of regular expression export
    }
  }
}


\section{Changes in Morpho version 0.25-1 (2013-08-02)}{
  \subsection{R-code}{
    \itemize{
      \item replaced abind with 'bindArr' in 'fixLMtps' and 'procSym'
      \item moved Distance calculation in 'CVA' to '.CVAdists' and radically cleaned code;
      \item removal of redundant variables in 'groupPCA'
      \item added option 'reflect' to 'procSym' (and underlying functions) and added copyright note for code copied from shapes package
    }
  }
  \subsection{Package documentation}{
    \itemize{
      \item added docu for new data set "boneData"
      \item added examples in 'procSym'
      \item added example for 'tps3d'
    }
  }
  \subsection{Misc}{
    \itemize{
      \item added additional example data (landmarks and a bone surface mesh)
      \item updated DESCRIPTION and README
      \item updated README added sourceforge and github download links
      \item added README.Rd
      \item added NEWS.Rd

    }
  }
}


\section{Changes in Morpho version 0.25 (2013-07-24)}{
  \subsection{R-code}{
    \itemize{
      \item simplified code in 'CVA' and 'CVA.crova'
      \item fixed issue when n-groups > n variables in 'CVA' and 'CVA.crova'
      \item replaced 'nn2' (RANN) with 'ann' (yaImpute) because RANN got (temporarily) banned from CRAN because of violation of policy
      \item added 'readLandmarks.csv' to read landmarks stored in spreadsheet text files (e.g. comma separated)
      \item added check for empty groups in 'typprobClass' and 'covW'
      \item changed cube creation in 'deformGrid3d' (much faster now).
      \item added options weights and centerweight in 'rotonmat'
      \item added confirmation for plotting spheres in 'deformGrid3d' when n > 1000 (can be slow)
      \item 'pcaplot3d' now works if data has only 1 PC
      \item removed useless code remains
      \item added 'bindArr' to concatenate matrices/arrays
      \item changed default core value in 'closemeshKD' from detectCores() to 1.
      \item added option 'keep' in 'rmVertex' 
      \item changed default method in 'meshDist.mesh3d' from vcglib to morpho
      \item renamed 'regdist.raw' to 'regdist' and removed original 'regdist' (too prone to user error)
    }
  }
  \subsection{Package documentation}{
    \itemize{
      \item finally mentioned Marta Rufino for providing all examples in docu for 'CVA' (sorry Marta, I simply forgot)
      \item fixed some docu and added examples in 'pcaplot3d' and
      'rotonmat'
      \item added docu for 'anonymize', 'r2morphoj', 'r2morphologika', 'NNshapeReg', 'kendalldist'
    }
  }
  \subsection{Misc}{
    \itemize{
      \item updated COPYRIGHT file
      \item switched from depends to import for colorRamps, MASS, RANN
    }
  }
}

