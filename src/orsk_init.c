#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .Fortran calls */
extern void F77_NAME(oddsratio)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
	    {"oddsratio", (DL_FUNC) &F77_NAME(oddsratio), 10},
	        {NULL, NULL, 0}
};

void R_init_orsk(DllInfo *dll)
{
	    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
	        R_useDynamicSymbols(dll, FALSE);
}

