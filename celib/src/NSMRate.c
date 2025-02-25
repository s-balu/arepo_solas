#include "config.h"
#include "NSMYields.h"
#include "NSMRate.h"

/*! \file NSMRate.c
 * \brief This file has functions regarding NSM rate.
 */

#define CELibNSMRateTableSize  (1000)
#define CELibNSMRateAgeMax     (1.38e10)

struct CELibStructNSMRate{
    double AgeMin;
    double AgeMax; // = 1.38e10; // Hubble time, following Plank's first date.
    double Rate[CELibNSMRateTableSize]; // The normalized cumulative NSM rate.
    double Age[CELibNSMRateTableSize];  // The Age bin for the cumulative NSM rate.
};

/*
 * CELibNSMRate holds the number of type NSMs as functions of age and
 * metallicity. The suffix is used to distinguish the amount of the metallicity.
 * For the case of the simple power law type rate, the only CELibNSMRate[0] is
 * used.
 */
struct CELibStructNSMRate *CELibNSMRate;
static int LifeTimeTableSizeMetallicity;

/*!
 * This fuction initializes NSM rate calculation routines. Before call this
 * function, an appropriate NSM mode should be selected.
 */
void CELibInitNSMRate(void){

    LifeTimeTableSizeMetallicity = CELibGetLifeTimeTableSizeMetallicity();
    CELibNSMRate = realloc(CELibNSMRate,
            sizeof(struct CELibStructNSMRate)*LifeTimeTableSizeMetallicity);

    if(CELibRunParameters.NSMType == CELibNSMRateModelID_PowerLaw){
        CELibNSMRate[0].AgeMin = CELibGetNSMEarliestExplosionTimePowerLaw();
        CELibNSMRate[0].AgeMax = 5*CELibNSMRateAgeMax; // Hubble time, following Plank's first date.
        double LogTimeBin[] = {log10(CELibNSMRate[0].AgeMin),log10(CELibNSMRate[0].AgeMax)};
        double dLogTime = (LogTimeBin[1]-LogTimeBin[0])/CELibNSMRateTableSize;
        for(int k=0;k<CELibNSMRateTableSize;k++){
            CELibNSMRate[0].Age[k] = pow(10.0,dLogTime*k+LogTimeBin[0]);
            CELibNSMRate[0].Rate[k] = CELibGetNSMIntegratedRatePowerLaw(CELibNSMRate[0].Age[k]);
        }
    } else {
        fprintf(stderr,"Incorrect NSMType = %d is selected.\n",CELibRunParameters.NSMType);
        fflush(stderr);
    }

    if(CELibRunParameters.TestMode){
        CELibWriteNSMRate("./CELib/CELibNSM");
    }

    return ;
}

/*!
 * This function writes a NSM rate as a function of age.
 */
void CELibWriteNSMRate(char OutDir[]){

    MakeDir(OutDir);

    FILE *fp;
    char fname[MaxCharactersInLine];
    if(CELibRunParameters.NSMType == CELibNSMRateModelID_PowerLaw){
        Snprintf(fname,"%s/CELibNSMRate.PowerLaw.dat",OutDir);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#CELib NSM rate based on the power law model.\n");
        fprintf(fp,"#No metallicity dependence.\n");
        fprintf(fp,"#Defined range = %g year - %g year\n",CELibNSMRate[0].AgeMin,CELibNSMRate[0].AgeMax);
        fprintf(fp,"#Age #NSMRate(Number/1Msun)\n");
        for(int k=0;k<CELibNSMRateTableSize;k++){
            fprintf(fp,"%g %g\n",CELibNSMRate[0].Age[k],CELibNSMRate[0].Rate[k]);
        }
        fclose(fp);
    } 

    return;
}


/*!
 * This function returns the explosion time of NSMs for an SSP particle
 * with a given initial mass and a given rate.
 */
static double pCELibGetNSMExplosionTimePowerLaw(const double InitialMass, const double Rate){

    int LValue = 0;
    int RValue = CELibNSMRateTableSize-1;
    if(Rate <  InitialMass*CELibNSMRate[0].Rate[LValue]) return CELibNSMRate[0].Age[LValue];
    if(Rate >  InitialMass*CELibNSMRate[0].Rate[RValue]) return 10*CELibNSMRateAgeMax;

    int Pivot = (LValue+RValue)>>1;
    do{
        if(Rate > InitialMass*CELibNSMRate[0].Rate[Pivot]){
            LValue = Pivot;
        } else {
            RValue = Pivot;
        }

        Pivot = (RValue+LValue)>>1;
        if(RValue - LValue == 1){
            double dagedrate = (CELibNSMRate[0].Age[RValue] - CELibNSMRate[0].Age[LValue])
                      /(InitialMass*(CELibNSMRate[0].Rate[RValue] - CELibNSMRate[0].Rate[LValue]));
            return dagedrate*(Rate - InitialMass*CELibNSMRate[0].Rate[LValue]) + CELibNSMRate[0].Age[LValue];
        }
    } while((RValue - LValue) != 0);
}


/*!
 * This function returns the next explosion time of NSMs for an SSP particle.
 * The next explosion time is evaluated according to the adopted DTD function.
 * If the expected explosion time becomes out of range, this function returns
 * -1. To evaluate the explosion time, this function requires a random number,
 * Rate, in the range of 0-1, the metallicity of the SSP particle Metallicity,
 * the initial mass of the SSP particle InitialMass and the count of NSM
 * explosion in the SSP particle Count. So far, the rate does not depend on the
 * metallicity and thus any value is available. This is reserved for future
 * extension. 
 */
double CELibGetNSMFeedbackTime(const double Rate, const double Metallicity, const double InitialMass, const int Count){

    double CurrentCount = CELibRunParameters.NSMNassociation*Count;
    double CurrentRate  = CELibRunParameters.NSMNassociation*Rate;

    if(CELibRunParameters.NSMType == CELibNSMRateModelID_PowerLaw){ // Power-law model
        return pCELibGetNSMExplosionTimePowerLaw(InitialMass,CurrentCount+CurrentRate);
    } else {
        return 0.0;
    }
}

