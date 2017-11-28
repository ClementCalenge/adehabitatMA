#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bufligr(void *, void *, void *, void *, void *, void *, void *, void *);
extern void distxyr(void *, void *, void *, void *, void *);
extern void erodil(void *, void *, void *, void *, void *);
extern void getcontourc(void *, void *, void *, void *, void *, void *);
extern void lcontour(void *, void *, void *, void *);
extern void regrouascnumr(void *, void *, void *, void *, void *, void *);
extern void regroufacascr(void *, void *, void *, void *, void *, void *, void *, void *);
extern void seqeticorr(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bufligr",       (DL_FUNC) &bufligr,       8},
    {"distxyr",       (DL_FUNC) &distxyr,       5},
    {"erodil",        (DL_FUNC) &erodil,        5},
    {"getcontourc",    (DL_FUNC) &getcontourc,    6},
    {"lcontour",      (DL_FUNC) &lcontour,      4},
    {"regrouascnumr", (DL_FUNC) &regrouascnumr, 6},
    {"regroufacascr", (DL_FUNC) &regroufacascr, 8},
    {"seqeticorr",    (DL_FUNC) &seqeticorr,    3},
    {NULL, NULL, 0}
};

void R_init_adehabitatMA(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
