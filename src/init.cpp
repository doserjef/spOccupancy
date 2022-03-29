#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spOccupancy.h"

static const R_CallMethodDef CallEntries[] = {
    {"PGOcc", (DL_FUNC) &PGOcc, 35},
    {"spPGOcc", (DL_FUNC) &spPGOcc, 50}, 
    {"spPGOccNNGP", (DL_FUNC) &spPGOccNNGP, 56},
    {"spPGOccPredict", (DL_FUNC) &spPGOccPredict, 15},
    {"spPGOccNNGPPredict", (DL_FUNC) &spPGOccNNGPPredict, 17},
    {"msPGOcc", (DL_FUNC) &msPGOcc, 43},
    {"spMsPGOcc", (DL_FUNC) &spMsPGOcc, 58},
    {"spMsPGOccNNGP", (DL_FUNC) &spMsPGOccNNGP, 64},
    {"spMsPGOccPredict", (DL_FUNC) &spMsPGOccPredict, 16},
    {"spMsPGOccNNGPPredict", (DL_FUNC) &spMsPGOccNNGPPredict, 18},
    {"intPGOcc", (DL_FUNC) &intPGOcc, 31},
    {"spIntPGOcc", (DL_FUNC) &spIntPGOcc, 46},
    {"spIntPGOccNNGP", (DL_FUNC) &spIntPGOccNNGP, 52},
    {"lfMsPGOcc", (DL_FUNC) &lfMsPGOcc, 44},
    {"sfMsPGOccNNGP", (DL_FUNC) &sfMsPGOccNNGP, 61},
    {"sfMsPGOccNNGPPredict", (DL_FUNC) &sfMsPGOccNNGPPredict, 20},
    {"lfJSDM", (DL_FUNC) &lfJSDM, 25},
    {"sfJSDMNNGP", (DL_FUNC) &sfJSDMNNGP, 42},
    {NULL, NULL, 0}
};

void R_init_spOccupancy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

