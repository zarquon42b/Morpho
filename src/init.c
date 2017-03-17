#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP addoCpp(SEXP);
extern SEXP ang_calcC(SEXP, SEXP);
extern SEXP ang_calcM(SEXP, SEXP);
extern SEXP armaGinvCpp(SEXP, SEXP);
extern SEXP arrMean3Cpp(SEXP);
extern SEXP asymPermuteCpp(SEXP, SEXP, SEXP);
extern SEXP barycenterCpp(SEXP, SEXP);
extern SEXP covPCAwrap(SEXP, SEXP, SEXP, SEXP);
extern SEXP createL(SEXP, SEXP);
extern SEXP CVAdists(SEXP, SEXP, SEXP, SEXP);
extern SEXP edgePlane(SEXP, SEXP, SEXP);
extern SEXP face_zero(SEXP);
extern SEXP fastSubsetMeans(SEXP, SEXP, SEXP, SEXP);
extern SEXP meshresCpp(SEXP, SEXP);
extern SEXP permudistArma(SEXP, SEXP, SEXP);
extern SEXP points2mesh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scaleprocCpp(SEXP);
extern SEXP tpsfx(SEXP, SEXP, SEXP, SEXP);
extern SEXP tweakU(SEXP, SEXP, SEXP, SEXP);
extern SEXP updateFaceNormals(SEXP, SEXP);
extern SEXP updateVertexNormals(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"addoCpp",                (DL_FUNC) &addoCpp,                1},
    {"ang_calcC",           (DL_FUNC) &ang_calcC,           2},
    {"ang_calcM",           (DL_FUNC) &ang_calcM,           2},
    {"armaGinvCpp",            (DL_FUNC) &armaGinvCpp,            2},
    {"arrMean3Cpp",            (DL_FUNC) &arrMean3Cpp,            1},
    {"asymPermuteCpp",         (DL_FUNC) &asymPermuteCpp,         3},
    {"barycenterCpp",          (DL_FUNC) &barycenterCpp,          2},
    {"covPCAwrap",          (DL_FUNC) &covPCAwrap,          4},
    {"createL",             (DL_FUNC) &createL,             2},
    {"CVAdists",            (DL_FUNC) &CVAdists,            4},
    {"edgePlane",           (DL_FUNC) &edgePlane,           3},
    {"face_zero",           (DL_FUNC) &face_zero,           1},
    {"fastSubsetMeans",     (DL_FUNC) &fastSubsetMeans,     4},
    {"meshresCpp",             (DL_FUNC) &meshresCpp,             2},
    {"permudistArma",       (DL_FUNC) &permudistArma,       3},
    {"points2mesh",         (DL_FUNC) &points2mesh,         8},
    {"scaleprocCpp",           (DL_FUNC) &scaleprocCpp,           1},
    {"tpsfx",               (DL_FUNC) &tpsfx,               4},
    {"tweakU",              (DL_FUNC) &tweakU,              4},
    {"updateFaceNormals",   (DL_FUNC) &updateFaceNormals,   2},
    {"updateVertexNormals", (DL_FUNC) &updateVertexNormals, 3},
    {NULL, NULL, 0}
};

void R_init_Morpho(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
