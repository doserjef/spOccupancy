#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spOccupancy.h"

static const R_CallMethodDef CallEntries[] = {
    {"PGOcc", (DL_FUNC) &PGOcc, 35},
    {"spPGOcc", (DL_FUNC) &spPGOcc, 52}, 
    {"spPGOccNNGP", (DL_FUNC) &spPGOccNNGP, 58},
    {"spPGOccPredict", (DL_FUNC) &spPGOccPredict, 15},
    {"spPGOccNNGPPredict", (DL_FUNC) &spPGOccNNGPPredict, 17},
    {"msPGOcc", (DL_FUNC) &msPGOcc, 43},
    {"spMsPGOcc", (DL_FUNC) &spMsPGOcc, 59},
    {"spMsPGOccNNGP", (DL_FUNC) &spMsPGOccNNGP, 65},
    {"spMsPGOccPredict", (DL_FUNC) &spMsPGOccPredict, 16},
    {"spMsPGOccNNGPPredict", (DL_FUNC) &spMsPGOccNNGPPredict, 18},
    {"intPGOcc", (DL_FUNC) &intPGOcc, 31},
    {"spIntPGOcc", (DL_FUNC) &spIntPGOcc, 48},
    {"spIntPGOccNNGP", (DL_FUNC) &spIntPGOccNNGP, 54},
    {"lfMsPGOcc", (DL_FUNC) &lfMsPGOcc, 44},
    {"sfMsPGOccNNGP", (DL_FUNC) &sfMsPGOccNNGP, 62},
    {"sfMsPGOccNNGPPredict", (DL_FUNC) &sfMsPGOccNNGPPredict, 20},
    {"lfJSDM", (DL_FUNC) &lfJSDM, 25},
    {"sfJSDMNNGP", (DL_FUNC) &sfJSDMNNGP, 44},
    {"tPGOcc", (DL_FUNC) &tPGOcc, 44},
    {"stPGOccNNGP", (DL_FUNC) &stPGOccNNGP, 64},
    {"stPGOccNNGPPredict", (DL_FUNC) &stPGOccNNGPPredict, 21},
    {"svcPGBinomNNGP", (DL_FUNC) &svcPGBinomNNGP, 45},
    {"svcPGOccNNGPPredict", (DL_FUNC) &svcPGOccNNGPPredict, 20},
    {"svcPGOccNNGP", (DL_FUNC) &svcPGOccNNGP, 59},
    {"svcTPGBinomNNGP", (DL_FUNC) &svcTPGBinomNNGP, 51},
    {"svcTPGOccNNGPPredict", (DL_FUNC) &svcTPGOccNNGPPredict, 24},
    {"svcTPGOccNNGP", (DL_FUNC) &svcTPGOccNNGP, 65},
    {"intMsPGOcc", (DL_FUNC) &intMsPGOcc, 46},
    {"postHocLM", (DL_FUNC) &postHocLM, 21},
    {"svcMsPGOccNNGP", (DL_FUNC) &svcMsPGOccNNGP, 64},
    {"svcMsPGOccNNGPPredict", (DL_FUNC) &svcMsPGOccNNGPPredict, 22},
    {"svcTMsPGOccNNGP", (DL_FUNC) &svcTMsPGOccNNGP, 61},
    {"svcTMsPGOccNNGPPredict", (DL_FUNC) &svcTMsPGOccNNGPPredict, 26},
    {"stMsPGOccNNGP", (DL_FUNC) &svcTMsPGOccNNGP, 58},
    {NULL, NULL, 0}
};

void R_init_spOccupancy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

