#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spOccupancy.h"

static const R_CallMethodDef CallEntries[] = {
    {"PGOcc", (DL_FUNC) &PGOcc, 19},
    {"PGOccREDet", (DL_FUNC) &PGOccREDet, 32},
    {"PGOccREOcc", (DL_FUNC) &PGOccREOcc, 32},
    {"PGOccREBoth", (DL_FUNC) &PGOccREBoth, 42},
    {"spPGOcc", (DL_FUNC) &spPGOcc, 34}, 
    {"spPGOccNNGP", (DL_FUNC) &spPGOccNNGP, 39},
    {"msPGOcc", (DL_FUNC) &msPGOcc, 28},
    {"spMsPGOcc", (DL_FUNC) &spMsPGOcc, 43},
    {"spMsPGOccNNGP", (DL_FUNC) &spMsPGOccNNGP, 49},
    {"intPGOcc", (DL_FUNC) &intPGOcc, 26},
    {"spIntPGOcc", (DL_FUNC) &spIntPGOcc, 41},
    {"spIntPGOccNNGP", (DL_FUNC) &spIntPGOccNNGP, 47},
    {NULL, NULL, 0}
};

void R_init_spOccupancy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

