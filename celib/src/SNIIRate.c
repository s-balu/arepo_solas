#include "config.h"
#include "SNIIYields.h"
#include "SNIIRate.h"
#include "LifeTime.h"

/*! \file SNIIRate.c
 * \brief This file has functions regarding SNII rate.
 */

#define CELibSNIIRateTableSize  (100)

struct CELibStructSNIIRate{
    double Metallicity; // Initial metallicity.
    double AgeMin;
    double AgeMax;
    double Rate[CELibSNIIRateTableSize]; // The normalized cumulative SNII rate.
    double Age[CELibSNIIRateTableSize];  // The Age bin for the cumulative SNII rate.
};

/*!
 * CELibSNIIRate holds the number of type II supernovae as functions of age and
 * metallicity. The suffix is used to distinguish the amount of the metallicity.
 */
struct CELibStructSNIIRate *CELibSNIIRate;

/*!
 * This function generates the SNII rate tables.
 */
void CELibInitSNIIRate(void){

    CELibSNIIRate = realloc(CELibSNIIRate,
            sizeof(struct CELibStructSNIIRate)*CELibSNIIYields_Metallicity);

    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        double Metallicity = CELibSNIIYieldsZ[i];
        double MassMax, MassMin;
    
        if(CELibRunParameters.SNIIExplosionTimeFollowYields == 1){
            MassMax = CELibGetSNIIYieldsMaxExplosionMass(i);
            MassMin = CELibGetSNIIYieldsMinExplosionMass(i);
        } else {
            MassMax = CELibRunParameters.SNIIUpperMass;
            MassMin = CELibRunParameters.SNIILowerMass;
        }

        CELibSNIIRate[i].AgeMin = CELibGetLifeTimeofStar(MassMax,Metallicity);
        CELibSNIIRate[i].AgeMax = CELibGetLifeTimeofStar(MassMin,Metallicity);

        if(CELibRunParameters.TestMode){
            fprintf(stderr,"//\t Show SNII Explosion time, yields = %d\n",CELibRunParameters.SNIIYieldsTableID);
            fprintf(stderr,"//\t\t SNII Explosion mass range for Z = %g = %g %g\n",
                    Metallicity,MassMin,MassMax);
            fprintf(stderr,"//\t\t  AgeMin(%g) = %g, AgeMax(%g) = %g\n",
                    MassMax,CELibSNIIRate[i].AgeMin,MassMin,CELibSNIIRate[i].AgeMax);
        }

        int CurrentIMF = CELibRunParameters.IMFType;
        if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
            if((CELibRunParameters.PopIIIIMF == 1)&&(Metallicity < CELibRunParameters.PopIIIMetallicity)){
                CurrentIMF = CELibIMF_Susa;
            }
        }

        double Normalization = IntegralSimpson(MassMin,MassMax,
                        CELibRunParameters.IntegrationSteps,CELibIMF[CurrentIMF].IMFFunctionPerMass);

        double dAge = (CELibSNIIRate[i].AgeMax-CELibSNIIRate[i].AgeMin)/(CELibSNIIRateTableSize-1);
        for(int k=0;k<CELibSNIIRateTableSize;k++){
            CELibSNIIRate[i].Age[k] = CELibSNIIRate[i].AgeMin + k*dAge;
            double Mass = CELibGetDyingStellarMass(CELibSNIIRate[i].Age[k],Metallicity);
            if(Mass == __CELib_NoValue__){
                CELibSNIIRate[i].Rate[k] = 0.e0;
            } else {
                CELibSNIIRate[i].Rate[k] = IntegralSimpson(Mass,MassMax,
                        CELibRunParameters.IntegrationSteps,CELibIMF[CurrentIMF].IMFFunctionPerMass);
            }
            CELibSNIIRate[i].Rate[k] /= Normalization;
        }
        CELibSNIIRate[i].Rate[CELibSNIIRateTableSize-1] = 1.0;
    }

    if(CELibRunParameters.TestMode){
        CELibWriteSNIIRate("./CELib/CELibSNIIRate");
        CELibTestSNIIRate("./CELib/CELibSNIIRate");
        CELibTestSNIIRateTable("./CELib/CELibSNIIRate");
    }

    return ;
}

/*!
 * This function writes all SNII rate tables in files.
 */
void CELibWriteSNIIRate(char OutDir[]){
    
    MakeDir(OutDir);
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/SNIIRate.%02d.dat",OutDir,i);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#CELib SNII Rate\n");
        fprintf(fp,"#Metallicity:%g\n",CELibSNIIYieldsZ[i]);
        for(int k=0;k<CELibSNIIRateTableSize;k++){
            fprintf(fp,"%g %g\n",CELibSNIIRate[i].Age[k],CELibSNIIRate[i].Rate[k]);
        }
        fclose(fp);
    }

    return ;
}

/*!
 * This function picks up the explosion time from the lifetime table. The
 * parameter Rate should be normalized and the range should be 0-1.0.
 */
static double pCELibGetSNIIExplosionTimeFromTable(const int TableID, const double Rate){

    int LValue = 0;
    int RValue = CELibSNIIRateTableSize-1;
    if(Rate >  CELibSNIIRate[TableID].Rate[RValue]) return NONE;
    if(Rate <  CELibSNIIRate[TableID].Rate[LValue]) return NONE;

    int Pivot = (LValue+RValue)>>1;
    do{
        if(Rate > CELibSNIIRate[TableID].Rate[Pivot]){
            LValue = Pivot;
        } else {
            RValue = Pivot;
        }
        Pivot = (RValue+LValue)>>1;
        if(RValue - LValue == 1){
            double dagedrate = (CELibSNIIRate[TableID].Age[RValue] - CELibSNIIRate[TableID].Age[LValue])
                            /(CELibSNIIRate[TableID].Rate[RValue] - CELibSNIIRate[TableID].Rate[LValue]);
            return dagedrate*(Rate - CELibSNIIRate[TableID].Rate[LValue]) + CELibSNIIRate[TableID].Age[LValue];
        }

    } while((RValue - LValue) != 0);
}

/*!
 * This function evaluates the explosion time of a SSP particle and then returns
 * it (in the unit of year).  The explosion time computed in this function is
 * IMF weighted one. If the expected explosion time becomes out of range, this
 * function returns NONE(=-1). To evaluate the explosion time, this function
 * requires (1) a random number, Rate, in the range of 0-1, and (2) Metallicity
 * of the SSP particle. Interpolation in metallicity is carried out in this
 * funtion automatically.
 */
double CELibGetSNIIExplosionTime(const double Rate, const double Metallicity){

    if(Metallicity <= CELibSNIIYieldsZ[0]){ // Use the lowest metallicity star's lifetime.
        return pCELibGetSNIIExplosionTimeFromTable(0,Rate);
    } else if(Metallicity >= CELibSNIIYieldsZ[CELibSNIIYields_Metallicity-1]){ // Use the highest metallicity star's lifetime. 
        return pCELibGetSNIIExplosionTimeFromTable(CELibSNIIYields_Metallicity-1,Rate);
    } else { // Use metallicity interpolated stellar lifetime.
        int IndexMetallicity = 0;
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(CELibSNIIYieldsZ[i] > Metallicity){
                IndexMetallicity = i-1;
                break;
            }
        }
        double ExplosionTime[] = {pCELibGetSNIIExplosionTimeFromTable(IndexMetallicity+1,Rate),
                                  pCELibGetSNIIExplosionTimeFromTable(IndexMetallicity,Rate)};
        if(ExplosionTime[0] < 0.0) return NONE;
        if(ExplosionTime[1] < 0.0) return NONE;
        double dagedmetal = (ExplosionTime[0]-ExplosionTime[1])
                            /(CELibSNIIYieldsZ[IndexMetallicity+1]-CELibSNIIYieldsZ[IndexMetallicity]);
        return dagedmetal*(Metallicity-CELibSNIIYieldsZ[IndexMetallicity])+ExplosionTime[1];
    }
}

/*!
 * This function writes the explosion time of type II SN as a function of rate.
 * Interpolation in metallicity is carried out.
 */
void CELibTestSNIIRate(char OutDir[]){

    const int NMetallicity = 10;
    const int NRate = 100;

    double MetallicityMin = 1.e-5;
    double MetallicityMax = 0.1;
    double MetallicityLogMin = log10(MetallicityMin);
    double MetallicityLogMax = log10(MetallicityMax);
    double dMetallicity = (MetallicityLogMax-MetallicityLogMin)/NMetallicity;
 
    MakeDir(OutDir);
    for(int i=0;i<NMetallicity;i++){
        double Metallicity = pow(10.0,dMetallicity*i+MetallicityLogMin);

        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/SNIIRateTest.%02d.dat",OutDir,i);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#CELib SNII Rate Test\n");
        fprintf(fp,"#Metallicity:%g\n",Metallicity);
        double dRate = 1.0/NRate;
        for(int k=0;k<NRate;k++){
            double Rate = dRate*k;
            fprintf(fp,"%g %g\n",CELibGetSNIIExplosionTime(Rate,Metallicity),Rate);
        }
        fclose(fp);
    }
    return ;
}

/*!
 * This function writes the explosion time of type II SN as a function of rate.
 * The table data are directly written.
 */
void CELibTestSNIIRateTable(char OutDir[]){

    const int NRate = 200;

    MakeDir(OutDir);
    for(int i=0;i<CELibSNIIYields_Metallicity;i++){
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/SNIIRateTableTest.%02d.dat",OutDir,i);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#CELib SNII Rate Test\n");
        fprintf(fp,"#Metallicity:%g\n",CELibSNIIYieldsZ[i]);
        double dRate = 1.0/NRate;
        for(int k=0;k<NRate;k++){
            double Rate = dRate*k;
            fprintf(fp,"%g %g\n",pCELibGetSNIIExplosionTimeFromTable(i,Rate),Rate);
        }
        fclose(fp);
    }
    return ;
}
