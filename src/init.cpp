#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spOccupancy.h"

static const R_CallMethodDef CallEntries[] = {
    {"PGOcc", (DL_FUNC) &PGOcc, 25},
    {"PGOccREDet", (DL_FUNC) &PGOccREDet, 35},
    {"PGOccREOcc", (DL_FUNC) &PGOccREOcc, 35},
    {"PGOccREBoth", (DL_FUNC) &PGOccREBoth, 45},
    {"spPGOcc", (DL_FUNC) &spPGOcc, 40}, 
    {"spPGOccRE", (DL_FUNC) &spPGOccRE, 50},
    {"spPGOccNNGP", (DL_FUNC) &spPGOccNNGP, 45},
    {"spPGOccNNGPRE", (DL_FUNC) &spPGOccNNGPRE, 55},
    {"spPGOccPredict", (DL_FUNC) &spPGOccPredict, 14},
    {"spPGOccNNGPPredict", (DL_FUNC) &spPGOccNNGPPredict, 16},
    {"msPGOcc", (DL_FUNC) &msPGOcc, 34},
    {"msPGOccREBoth", (DL_FUNC) &msPGOccREBoth, 54},
    {"msPGOccREDet", (DL_FUNC) &msPGOccREDet, 44},
    {"msPGOccREOcc", (DL_FUNC) &msPGOccREOcc, 44},
    {"spMsPGOcc", (DL_FUNC) &spMsPGOcc, 49},
    {"spMsPGOccRE", (DL_FUNC) &spMsPGOccRE, 59},
    {"spMsPGOccNNGP", (DL_FUNC) &spMsPGOccNNGP, 55},
    {"spMsPGOccNNGPRE", (DL_FUNC) &spMsPGOccNNGPRE, 65},
    {"spMsPGOccPredict", (DL_FUNC) &spMsPGOccPredict, 15},
    {"spMsPGOccNNGPPredict", (DL_FUNC) &spMsPGOccNNGPPredict, 17},
    {"intPGOcc", (DL_FUNC) &intPGOcc, 31},
    {"spIntPGOcc", (DL_FUNC) &spIntPGOcc, 46},
    {"spIntPGOccNNGP", (DL_FUNC) &spIntPGOccNNGP, 52},
    {NULL, NULL, 0}
};

void R_init_spOccupancy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

