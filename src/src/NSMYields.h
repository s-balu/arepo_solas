#pragma once 

struct CELibStructNSMYields{
    char Name[MaxCharactersInLine];
    double Mej;
    double Energy;
    double Elements[CELibYield_Number];
};

void CELibInitNSMYields(void);
void CELibGetNSMModelNames(char *ModelNames);
void CELibGetNSMCurrentModelName(char *ModelName);
void CELibGetNSMCurrentYieldModelName(char *ModelName);
double CELibGetNSMIntegratedRatePowerLaw(const double Age);
double CELibGetNSMEarliestExplosionTimePowerLaw(void);

