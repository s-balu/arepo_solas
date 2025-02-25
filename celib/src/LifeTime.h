#pragma once

#define CELIB_LIFETIME_Z_P98  5
#define CELIB_LIFETIME_M_P98  30

#define CELIB_LIFETIME_Z_S02  1
#define CELIB_LIFETIME_M_S02  4

extern double *CELibLifeTimeZ;
extern double *CELibLifeTimeLogZ;
extern double CELibLifeTimeMass[CELIB_LIFETIME_M_P98];
extern double CELibLifeTimeZMF[CELIB_LIFETIME_Z_P98][CELIB_LIFETIME_M_P98];
