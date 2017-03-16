#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
    Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BLPestimatoR_getDelta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); extern SEXP BLPestimatoR_getExpMu(SEXP, SEXP, SEXP, SEXP, SEXP); extern SEXP BLPestimatoR_getSijMod(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
     {"BLPestimatoR_getDelta",  (DL_FUNC) &BLPestimatoR_getDelta,  13},
     {"BLPestimatoR_getExpMu",  (DL_FUNC) &BLPestimatoR_getExpMu,   5},
     {"BLPestimatoR_getSijMod", (DL_FUNC) &BLPestimatoR_getSijMod,  3},
     {NULL, NULL, 0}
};

void R_init_BLPestimatoR(DllInfo *dll)
{
     R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
     R_useDynamicSymbols(dll, FALSE);
}
