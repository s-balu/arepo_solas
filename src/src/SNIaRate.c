#include "config.h"
#include "SNIaYields.h"
#include "SNIaRate.h"

/*! \file SNIaRate.c
 * \brief This file holds functions regarding SNIa rate.
 */

#define CELibSNIaRateTableSize  (1000)
#define CELibSNIaRateAgeMax     (1.38e10)

struct CELibStructSNIaRate{
    double AgeMin;
    double AgeMax; // = 1.38e10; // Hubble time, following Plank's first date.
    double Rate[CELibSNIaRateTableSize]; // The normalized cumulative SNIa rate.
    double Age[CELibSNIaRateTableSize];  // The Age bin for the cumulative SNIa rate.
};

/*!
 * CELibSNIaRate holds the number of type Ia supernovae as functions of age and
 * metallicity. The suffix is used to distinguish the amount of the metallicity.
 * Since the power-law type models do not consider a metallicity dependence,
 * only CELibSNIaRate[0] is used.
 */
struct CELibStructSNIaRate *CELibSNIaRate;
static int LifeTimeTableSizeMetallicity;

/*!
 * This function initializes the SNIa rate calculation routine. 
 */
void CELibInitSNIaRate(void){

    LifeTimeTableSizeMetallicity = CELibGetLifeTimeTableSizeMetallicity();
    CELibSNIaRate = realloc(CELibSNIaRate,
            sizeof(struct CELibStructSNIaRate)*LifeTimeTableSizeMetallicity);

    if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_GR83){ // Greggio & Renzini model.
        for(int i=0;i<LifeTimeTableSizeMetallicity;i++){
            double Metallicity = CELibLifeTimeZ[i];
            CELibSNIaRate[i].AgeMin = CELibGetSNIaEarliestExplosionTimePortinari(Metallicity);
            CELibSNIaRate[i].AgeMax = CELibSNIaRateAgeMax; // Hubble time, following Plank's first year date.
            double LogTimeBin[] = {log10(CELibSNIaRate[i].AgeMin),log10(CELibSNIaRate[i].AgeMax)};
            double dLogTime = (LogTimeBin[1]-LogTimeBin[0])/CELibSNIaRateTableSize;
            for(int k=0;k<CELibSNIaRateTableSize;k++){
                CELibSNIaRate[i].Age[k] = pow(10.0,dLogTime*k+LogTimeBin[0]);
                CELibSNIaRate[i].Rate[k] = CELibGetSNIaIntegratedRateGreggioRenzini(CELibSNIaRate[i].Age[k],Metallicity);
            }
        }
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLaw){ // Power-law model based on Maoz & Mannucci (2012).
        CELibSNIaRate[0].AgeMin = CELibGetSNIaEarliestExplosionTimePowerLaw();
        CELibSNIaRate[0].AgeMax = CELibSNIaRateAgeMax; // Hubble time, following Plank's first year date.
        double LogTimeBin[] = {log10(CELibSNIaRate[0].AgeMin),log10(CELibSNIaRate[0].AgeMax)};
        double dLogTime = (LogTimeBin[1]-LogTimeBin[0])/CELibSNIaRateTableSize;
        for(int k=0;k<CELibSNIaRateTableSize;k++){
            CELibSNIaRate[0].Age[k] = pow(10.0,dLogTime*k+LogTimeBin[0]);
            CELibSNIaRate[0].Rate[k] = CELibGetSNIaIntegratedRatePowerLaw(CELibSNIaRate[0].Age[k]);
        }
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLawV13){ // Power-law model based on Vogelsberger et al (2013).
        CELibSNIaRate[0].AgeMin = CELibGetSNIaEarliestExplosionTimeVogelsberger();
        CELibSNIaRate[0].AgeMax = CELibSNIaRateAgeMax; // Hubble time, following Plank's first year date.
        double LogTimeBin[] = {log10(CELibSNIaRate[0].AgeMin),log10(CELibSNIaRate[0].AgeMax)};
        double dLogTime = (LogTimeBin[1]-LogTimeBin[0])/CELibSNIaRateTableSize;
        for(int k=0;k<CELibSNIaRateTableSize;k++){
            CELibSNIaRate[0].Age[k] = pow(10.0,dLogTime*k+LogTimeBin[0]);
            CELibSNIaRate[0].Rate[k] = CELibGetSNIaIntegratedRateVogelsberger(CELibSNIaRate[0].Age[k]);
        }
    } else {
        fprintf(stderr,"Incorrect SNIaType = %d is used.\n",CELibRunParameters.SNIaType);
        fflush(stderr);
    }

    if(CELibRunParameters.TestMode){
        CELibWriteSNIaRate("./CELib/CELibSNIaRate");
    }

    return ;
}

/*!
 * This function writes SNIa rate as a function of age.
 */
void CELibWriteSNIaRate(char OutDir[]){

    MakeDir(OutDir);

    if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_GR83){
        for(int i=0;i<LifeTimeTableSizeMetallicity;i++){
            FILE *fp;
            char fname[MaxCharactersInLine];
            Snprintf(fname,"%s/CELibSNIaRate.%02d.dat",OutDir,i);
            FileOpen(fp,fname,"w");
            double Metallicity = CELibLifeTimeZ[i];
            fprintf(fp,"#CELib SNIa rate based on the Portinari model.\n");
            fprintf(fp,"#Metallicity:%g",Metallicity);
            fprintf(fp,"#Defined range = %g year - %g year\n",CELibSNIaRate[i].AgeMin,CELibSNIaRate[i].AgeMax);
            for(int k=0;k<CELibSNIaRateTableSize;k++){
                fprintf(fp,"%g %g\n",CELibSNIaRate[i].Age[k],CELibSNIaRate[i].Rate[k]);
            }
            fclose(fp);
        }
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLaw){
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/CELibSNIaRate.PowerLaw.dat",OutDir);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#CELib SNIa rate based on the power-law.\n");
        fprintf(fp,"#No metallicity dependence.\n");
        fprintf(fp,"#Defined range = %g year - %g year\n",CELibSNIaRate[0].AgeMin,CELibSNIaRate[0].AgeMax);
        fprintf(fp,"#Age #SNIaRate(Number/1Msun)\n");
        for(int k=0;k<CELibSNIaRateTableSize;k++){
            fprintf(fp,"%g %g\n",CELibSNIaRate[0].Age[k],CELibSNIaRate[0].Rate[k]);
        }
        fclose(fp);
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLawV13){
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/CELibSNIaRate.Vogel.dat",OutDir);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#CELib SNIa rate based on the Vogelsberger model.\n");
        fprintf(fp,"#No metallicity dependence.\n");
        fprintf(fp,"#Defined range = %g year - %g year\n",CELibSNIaRate[0].AgeMin,CELibSNIaRate[0].AgeMax);
        fprintf(fp,"#Age #SNIaRate(Number/1Msun)\n");
        for(int k=0;k<CELibSNIaRateTableSize;k++){
            fprintf(fp,"%g %g\n",CELibSNIaRate[0].Age[k],CELibSNIaRate[0].Rate[k]);
        }
        fclose(fp);
    }

    return;
}


/*!
 * This function returns the explosion time of type Ia SNe for an SSP particle
 * with a given initial mass and a given rate. The time is based on the
 * Portinari et al. (1998)'s model.
 */
static double pCELibGetSNIaExplosionTimePortinariFromTable(const int TableID, const double InitialMass, const double Rate){

    int LValue = 0;
    int RValue = CELibSNIaRateTableSize-1;
    if(Rate <  InitialMass*CELibSNIaRate[TableID].Rate[LValue]) return CELibSNIaRate[TableID].Age[LValue];
    if(Rate >  InitialMass*CELibSNIaRate[TableID].Rate[RValue]) return 10*CELibSNIaRateAgeMax;

    int Pivot = (LValue+RValue)>>1;
    do{
        if(Rate > InitialMass*CELibSNIaRate[TableID].Rate[Pivot]){
            LValue = Pivot;
        } else {
            RValue = Pivot;
        }

        Pivot = (RValue+LValue)>>1;
        if(RValue - LValue == 1){
            double dagedrate = (CELibSNIaRate[TableID].Age[RValue] - CELibSNIaRate[TableID].Age[LValue])
                      /(InitialMass*(CELibSNIaRate[TableID].Rate[RValue] - CELibSNIaRate[TableID].Rate[LValue]));
            return dagedrate*(Rate - InitialMass*CELibSNIaRate[TableID].Rate[LValue]) + CELibSNIaRate[TableID].Age[LValue];
        }
    } while((RValue - LValue) != 0);
}


/*!
 * This function returns the explosion time of type Ia SNe for an SSP particle
 * with a given initial mass and a given rate. The time is based on power-law
 * models.
 */
static double pCELibGetSNIaExplosionTimePowerLaw(const double InitialMass, const double Rate){

    int LValue = 0;
    int RValue = CELibSNIaRateTableSize-1;
    if(Rate <  InitialMass*CELibSNIaRate[0].Rate[LValue]) return CELibSNIaRate[0].Age[LValue];
    if(Rate >  InitialMass*CELibSNIaRate[0].Rate[RValue]) return 10*CELibSNIaRateAgeMax;

    int Pivot = (LValue+RValue)>>1;
    do{
        if(Rate > InitialMass*CELibSNIaRate[0].Rate[Pivot]){
            LValue = Pivot;
        } else {
            RValue = Pivot;
        }

        Pivot = (RValue+LValue)>>1;
        if(RValue - LValue == 1){
            double dagedrate = (CELibSNIaRate[0].Age[RValue] - CELibSNIaRate[0].Age[LValue])
                      /(InitialMass*(CELibSNIaRate[0].Rate[RValue] - CELibSNIaRate[0].Rate[LValue]));
            return dagedrate*(Rate - InitialMass*CELibSNIaRate[0].Rate[LValue]) + CELibSNIaRate[0].Age[LValue];
        }
    } while((RValue - LValue) != 0);
}


/*!
 * This function returns the next explosion time of type Ia SNe for an SSP
 * particle. To estimate the explosion time, the DTD, which is selected in
 * CELibRunParameter.SNIaType flags, is used. If the expected explosion time
 * becomes out of range, this function returns -1. To evaluate the explosion
 * time, this function requires one random number, Rate, in the range of 0-1,
 * the metallicity of the SSP particle, Metallicity, the initial mass of the SSP
 * particle, InitialMass, and the count of Ia explosion in the SSP particle,
 * Count. If the metallicity dependent model (the Portinari et al. 1998's model)
 * is adopted, an interpolation is used in order to obtain the explosion time
 * for an intermediate metallicity.  This is done in this funtion automatically.
 */
double CELibGetSNIaExplosionTime(const double Rate, const double Metallicity, const double InitialMass, const int Count){

    double CurrentCount = CELibRunParameters.SNIaNassociation*Count;
    double CurrentRate  = CELibRunParameters.SNIaNassociation*Rate;

    if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_GR83){ // Greggio & Renzini 1983 model
        if(Metallicity <= CELibLifeTimeZ[0]){
            return pCELibGetSNIaExplosionTimePortinariFromTable(0,InitialMass,CurrentCount+CurrentRate);
        } else if(Metallicity >= CELibLifeTimeZ[LifeTimeTableSizeMetallicity-1]){
            return pCELibGetSNIaExplosionTimePortinariFromTable(LifeTimeTableSizeMetallicity-1,InitialMass,CurrentCount+CurrentRate);
        } else {
            int IndexMetallicity = 0;
            for(int i=1;i<LifeTimeTableSizeMetallicity;i++){
                if(CELibLifeTimeZ[i] > Metallicity){
                    IndexMetallicity = i-1;
                    break;
                }
            }
            double LogMetallicity = log10(Metallicity);
            double ExplosionTime[] = {pCELibGetSNIaExplosionTimePortinariFromTable(IndexMetallicity+1,InitialMass,CurrentCount+CurrentRate),
                                      pCELibGetSNIaExplosionTimePortinariFromTable(IndexMetallicity,InitialMass,CurrentCount+CurrentRate)};
            if(ExplosionTime[0] < 0.e0) return NONE;
            if(ExplosionTime[1] < 0.e0) return NONE;
            double dagedmetal = (ExplosionTime[0]-ExplosionTime[1])
                                /(CELibLifeTimeLogZ[IndexMetallicity+1]-CELibLifeTimeLogZ[IndexMetallicity]);
            return dagedmetal*(LogMetallicity-CELibLifeTimeLogZ[IndexMetallicity])+ExplosionTime[1];
        }
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLaw){ // Power-law model
        return pCELibGetSNIaExplosionTimePowerLaw(InitialMass,CurrentCount+CurrentRate);
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLawV13){ // Power-law model (Vogelsberger et al. 2013)
        return pCELibGetSNIaExplosionTimePowerLaw(InitialMass,CurrentCount+CurrentRate);
    } else {
        return 0.0;
    }
}

