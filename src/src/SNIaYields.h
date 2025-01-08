#pragma once

struct CELibStructSNIaYields{
    char Name[MaxCharactersInLine];
    double Mej;
    double Energy;
    double Elements[CELibYield_Number];
};

void CELibInitSNIaYields(void);
void CELibGetSNIaModelNames(char *ModelNames);
void CELibGetSNIaCurrentModelName(char *ModelName);
void CELibGetSNIaCurrentYieldModelName(char *ModelName);
double CELibGetSNIaIntegratedRateGreggioRenzini(const double Age, const double Metallicity);
double CELibGetSNIaIntegratedRatePowerLaw(const double Age);
double CELibGetSNIaIntegratedRateVogelsberger(const double Age);

double CELibGetSNIaEarliestExplosionTimePortinari(const double Metallicity);
double CELibGetSNIaEarliestExplosionTimePowerLaw(void);
double CELibGetSNIaEarliestExplosionTimeVogelsberger(void);

