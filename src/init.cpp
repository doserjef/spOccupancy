#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spOccupancy.h"

static const R_CallMethodDef CallEntries[] = {
    {"PGOcc", (DL_FUNC) &PGOcc, 23},
    {"PGOccREDet", (DL_FUNC) &PGOccREDet, 33},
    {"PGOccREOcc", (DL_FUNC) &PGOccREOcc, 33},
    {"PGOccREBoth", (DL_FUNC) &PGOccREBoth, 43},
    {"spPGOcc", (DL_FUNC) &spPGOcc, 38}, 
    {"spPGOccRE", (DL_FUNC) &spPGOccRE, 48},
    {"spPGOccNNGP", (DL_FUNC) &spPGOccNNGP, 43},
    {"spPGOccNNGPRE", (DL_FUNC) &spPGOccNNGPRE, 53},
    {"spPGOccPredict", (DL_FUNC) &spPGOccPredict, 14},
    {"spPGOccNNGPPredict", (DL_FUNC) &spPGOccNNGPPredict, 16},
    {"msPGOcc", (DL_FUNC) &msPGOcc, 32},
    {"msPGOccREBoth", (DL_FUNC) &msPGOccREBoth, 52},
    {"msPGOccREDet", (DL_FUNC) &msPGOccREDet, 42},
    {"msPGOccREOcc", (DL_FUNC) &msPGOccREOcc, 42},
    {"spMsPGOcc", (DL_FUNC) &spMsPGOcc, 47},
    {"spMsPGOccRE", (DL_FUNC) &spMsPGOccRE, 57},
    {"spMsPGOccNNGP", (DL_FUNC) &spMsPGOccNNGP, 53},
    {"spMsPGOccNNGPRE", (DL_FUNC) &spMsPGOccNNGPRE, 63},
    {"spMsPGOccPredict", (DL_FUNC) &spMsPGOccPredict, 15},
    {"spMsPGOccNNGPPredict", (DL_FUNC) &spMsPGOccNNGPPredict, 17},
    {"intPGOcc", (DL_FUNC) &intPGOcc, 29},
    {"spIntPGOcc", (DL_FUNC) &spIntPGOcc, 44},
    {"spIntPGOccNNGP", (DL_FUNC) &spIntPGOccNNGP, 50},
    {NULL, NULL, 0}
};

void R_init_spOccupancy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

