#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _LinCDE_LinCDECdf(SEXP, SEXP, SEXP);
extern SEXP _LinCDE_LinCDEQuantiles(SEXP, SEXP, SEXP);
extern SEXP _LinCDE_LinCDESplit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _LinCDE_rcpp_hello_world();

static const R_CallMethodDef CallEntries[] = {
    {"_LinCDE_LinCDECdf",        (DL_FUNC) &_LinCDE_LinCDECdf,        3},
    {"_LinCDE_LinCDEQuantiles",  (DL_FUNC) &_LinCDE_LinCDEQuantiles,  3},
    {"_LinCDE_LinCDESplit",      (DL_FUNC) &_LinCDE_LinCDESplit,      7},
    {"_LinCDE_rcpp_hello_world", (DL_FUNC) &_LinCDE_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

void R_init_LinCDE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
