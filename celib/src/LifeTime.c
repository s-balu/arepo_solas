#include "config.h"
#include "LifeTimeTable.h"

double (*CELibGetLifeTimeofStar)(const double, const double);
double (*CELibGetDyingStellarMass)(const double, const double);

/*!
 * Sizes of the lifetime table. 
 */
int CELibLifeTime_Metallicity; 
int CELibLifeTime_Mass; 

double *CELibLifeTimeZ;
double *CELibLifeTimeLogZ;

static void pCELibFittingLifeTime(void);


int CELibGetLifeTimeTableSizeMetallicity(void){
    return CELibLifeTime_Metallicity;
}

/*!
 * This function initializes all data concerning the stellar lifetimes.
 */
void CELibInitLifeTime(void){


    CELibLifeTime_Metallicity = CELIB_LIFETIME_Z_P98; 
    free(CELibLifeTimeZ);
    free(CELibLifeTimeLogZ);
    CELibLifeTimeZ = malloc(sizeof(double)*CELibLifeTime_Metallicity);
    CELibLifeTimeLogZ = malloc(sizeof(double)*CELibLifeTime_Metallicity);

    CELibLifeTime_Mass = CELIB_LIFETIME_M_P98;

    for(int i=0;i<CELibLifeTime_Metallicity;i++){
        CELibLifeTimeZ[i] = CELibLifeTimeZ_P98[i];
        CELibLifeTimeLogZ[i] = log10(CELibLifeTimeZ_P98[i]);
    }

    // Set function pointers.
    if(CELibRunParameters.LifeTimeType == CELibLifeTime_Original){ 
        // Adopt table data.
        CELibGetLifeTimeofStar = CELibGetLifeTimeofStarTable;
        CELibGetDyingStellarMass = CELibGetDyingStellarMassTable;
    } else if(CELibRunParameters.LifeTimeType == CELibLifeTime_LSF){
        // Adopt fitted data.
        CELibGetLifeTimeofStar = CELibGetLifeTimeofStarLSF;
        CELibGetDyingStellarMass = CELibGetDyingStellarMassLSF;
    }

    // Evaluate fitting coefficients.
    pCELibFittingLifeTime();

    return ;
}


/*!
 * This function returns an interpolated value of the mass whose lifetime is
 * just finished at the given age of a SSP particle. The parameter "LifeTime", which
 * represents the age of the SSP particle, should be in the unit of "year".
 */
static inline double pCELibGetDyingStellarMassFromTable(const int TableID, const double LifeTime){

    int LValue = 0;
    int RValue = CELIB_LIFETIME_M_P98-1;
    if(LifeTime <= CELibLifeTimeZMF_P98[TableID][RValue]) return __CELib_NoValue__;
    if(LifeTime > CELibLifeTimeZMF_P98[TableID][LValue]) return __CELib_NoValue__;

    int Pivot = (LValue+RValue)>>1;
    do{
        if(LifeTime < CELibLifeTimeZMF_P98[TableID][Pivot]){
            LValue = Pivot;
        } else {
            RValue = Pivot;
        }

        Pivot = (RValue+LValue)>>1;
        if(RValue - LValue == 1){
            double dmdage = (CELibLifeTimeMass_P98[RValue] - CELibLifeTimeMass_P98[LValue])
                            /(CELibLifeTimeZMF_P98[TableID][RValue] - CELibLifeTimeZMF_P98[TableID][LValue]);
            return dmdage*(LifeTime - CELibLifeTimeZMF_P98[TableID][LValue]) + CELibLifeTimeMass_P98[LValue];
        }

    } while((RValue - LValue) != 0);

    return CELibLifeTimeMass_P98[RValue];
}

/*!
 * This function returns the mass of the star whose life is just finished at the
 * given age. The matallicity dependence is taken into account. The parameter
 * "LifeTime" should be given in the unit of "year", while the parameter
 * "metallicity" should be [0,1].
 */
double CELibGetDyingStellarMassTable(const double LifeTime, const double Metallicity){

    if(Metallicity <= CELibLifeTimeZ[0]){
        return pCELibGetDyingStellarMassFromTable(0,LifeTime);
    } else if(Metallicity >= CELibLifeTimeZ[CELibLifeTime_Metallicity-1]){
        return pCELibGetDyingStellarMassFromTable(CELibLifeTime_Metallicity-1,LifeTime);
    } else {
        int TableID = 0;
        for(int i=1;i<CELibLifeTime_Metallicity;i++){
            if(CELibLifeTimeZ_P98[i] > Metallicity){
                TableID = i-1;
                break;
            }
        }
        double MassLower = pCELibGetDyingStellarMassFromTable(TableID,LifeTime);
        double MassUpper = pCELibGetDyingStellarMassFromTable(TableID+1,LifeTime);
        double grad = ((MassUpper - MassLower)/(CELibLifeTimeZ_P98[TableID+1] - CELibLifeTimeZ_P98[TableID]));
        return grad*(Metallicity - CELibLifeTimeZ_P98[TableID]) + MassLower;
    }
}


/*!
 * This function returns the lifetime of a star of which mass is "Mass".  "Mass"
 * should be in the unit of the solar mass. "TableID" is the metallicity ID of
 * the look up table.
 */
static inline double pCELibGetLifeTimeofStarFromTable(const int TableID, const double Mass){

    int LValue = 0;
    int RValue = CELibLifeTime_Mass-1;
    if(Mass < CELibLifeTimeMass_P98[LValue]) return CELibLifeTimeZMF_P98[TableID][LValue];
    if(Mass > CELibLifeTimeMass_P98[RValue]) return CELibLifeTimeZMF_P98[TableID][RValue];

    int Pivot = (LValue+RValue)>>1;
    do{
        if(Mass < CELibLifeTimeMass_P98[Pivot]){
            RValue = Pivot;
        } else { 
            LValue = Pivot;
        }

        Pivot = (RValue+LValue)>>1;
        if(RValue - LValue == 1){
            double dtdm = (CELibLifeTimeZMF_P98[TableID][RValue] - CELibLifeTimeZMF_P98[TableID][LValue])
                            /(CELibLifeTimeMass_P98[RValue] - CELibLifeTimeMass_P98[LValue]);
            return dtdm*(Mass - CELibLifeTimeMass_P98[LValue]) + CELibLifeTimeZMF_P98[TableID][LValue];
        }

    } while((RValue - LValue) != 0);

    return CELibLifeTimeMass_P98[RValue];
}

/*!
 * This function returns the lifetime of a star (year) whose mass and
 * metallicity are "Mass" and "Metallicity", respectively. "Mass" should be in
 * the unit of the solar mass.
 */
double CELibGetLifeTimeofStarTable(const double Mass, const double Metallicity){

    if(Metallicity <= CELibLifeTimeZ[0]){
        return pCELibGetLifeTimeofStarFromTable(0,Mass);
    } else if(Metallicity >= CELibLifeTimeZ[CELibLifeTime_Metallicity-1]){
        return pCELibGetLifeTimeofStarFromTable(CELibLifeTime_Metallicity-1,Mass);
    } else {
        int TableID = 0;
        for(int i=1;i<CELibLifeTime_Metallicity;i++){
            if(CELibLifeTimeZ_P98[i] > Metallicity){
                TableID = i-1;
                break;
            }
        }
        double AgeLower = pCELibGetLifeTimeofStarFromTable(TableID,Mass);
        double AgeUpper = pCELibGetLifeTimeofStarFromTable(TableID+1,Mass);
        double grad = ((AgeUpper-AgeLower)/(CELibLifeTimeZ_P98[TableID+1]-CELibLifeTimeZ_P98[TableID]));
        return grad*(Metallicity-CELibLifeTimeZ_P98[TableID])+AgeLower;
    }
}

/*
 * @def
 * The order of least square fitting.
 */
#define FitDim (6)

/*
 * Valiables below are the pointer and its sizes.
 */
static double **FittingCoef;   // LifeTime
static double **FittingCoefMd; // Mass 

/*!
 * This fucntion evaluates coefficients for lifetime tables using the least
 * square fitting. 
 */
static void pCELibFittingLifeTimei(const double Mass[restrict], const double LT[restrict], const int Length, 
        double FittingCoef[restrict], double FittingCoefMd[restrict]){

    double a[FitDim+1][FitDim+2],s[2*FitDim+1],t[FitDim+1];
    double am[FitDim+1][FitDim+2],sm[2*FitDim+1],tm[FitDim+1];

    for(int k=0;k<2*FitDim+1;k++){
        s[k] = 0.e0;
        sm[k] = 0.e0;
    }
    for(int k=0;k<FitDim+1;k++){
        t[k] = 0.e0;
        tm[k] = 0.e0;
    }
    
    for(int k=0;k<Length;k++){
        for(int j=0;j<2*FitDim+1;j++){
            s[j] += pow(log10(Mass[k]),j);
            sm[j] += pow(log10(LT[k]),j);
        }
        for(int j=0;j<FitDim+1;j++){
            t[j] += pow(log10(Mass[k]),j)*log10(LT[k]);
            tm[j] += pow(log10(LT[k]),j)*log10(Mass[k]);
        }
    }

    for(int k=0;k<FitDim+1;k++){
        for(int j=0;j<FitDim+1;j++){
            a[k][j] = s[k+j];
            am[k][j] = sm[k+j];
        }
        a[k][FitDim+1] = t[k];
        am[k][FitDim+1] = tm[k];
    }

    for(int k=0;k<FitDim+1;k++){
        double p = a[k][k];
        double pm = am[k][k];
        for(int j=k;j<FitDim+2;j++){
            a[k][j] = a[k][j]/p;
            am[k][j] = am[k][j]/pm;
        }
        for(int j=0;j<FitDim+1;j++){
            if(j!=k){
                double d = a[j][k];
                double dm = am[j][k];
                for(int l=k;l<FitDim+2;l++){
                    a[j][l] = a[j][l]-d*a[k][l];
                    am[j][l] = am[j][l]-dm*am[k][l];
                }
            }
        }

    }

    for(int k=0;k<FitDim+1;k++){
        FittingCoef[k] = a[k][FitDim+1];
        FittingCoefMd[k] = am[k][FitDim+1];
    }

    return ;
}

/*!
 * This function evaluates coefficients of stellar lifetime data for the whole
 * metallicity range.
 */
static void pCELibFittingLifeTime(void){

    // Allocate the array for LSF coefficients.
    FittingCoef = malloc(sizeof(double *)*CELibLifeTime_Metallicity);
    FittingCoefMd = malloc(sizeof(double *)*CELibLifeTime_Metallicity);
    for(int i=0;i<CELibLifeTime_Metallicity;i++){
        FittingCoef[i] = malloc(sizeof(double)*(FitDim+1));
        FittingCoefMd[i] = malloc(sizeof(double)*(FitDim+1));
    }

    // Set mass.
    double tmpMass[CELIB_LIFETIME_M_P98+CELIB_LIFETIME_M_S02];
    for(int i=0;i<CELIB_LIFETIME_M_P98;i++){
        tmpMass[i] = CELibLifeTimeMass_P98[i];
    }
    for(int i=0;i<CELIB_LIFETIME_M_S02;i++){
        tmpMass[CELIB_LIFETIME_M_P98+i] = CELibLifeTimeMass_S02[i];
    }
    
    double tmpLT[CELIB_LIFETIME_M_P98+CELIB_LIFETIME_M_S02];

    for(int i=0;i<CELibLifeTime_Metallicity;i++){
        int Length;
    
        for(int k=0;k<CELibLifeTime_Mass;k++){
            tmpLT[k] = CELibLifeTimeZMF_P98[i][k];
        }
        Length = CELibLifeTime_Mass;
    
        pCELibFittingLifeTimei(tmpMass,tmpLT,Length,FittingCoef[i],FittingCoefMd[i]);
    }

    return ;
}

/*!
 * This function returns the lifetime of a star whose mass is "Mass" based on
 * the fitted lifetime data. 
 */
static inline double pCELibGetLifeTimeofStarUsingLSF(const int TableID, const double Mass){


    double LogMass = log10(Mass);
    double LogAge=0.e0;
    for(int k=0;k<FitDim+1;k++){
        LogAge += FittingCoef[TableID][k]*pow(LogMass,k);
    }
    return pow(10.0,LogAge);
}

/*!
 * This function returns the stellar lifetime (year) of a star whose mass and
 * metallicity are "Mass" and "Metallicity", respectively, based on the least
 * square fitting results.  "Mass" is in the unit of the solar mass. 
 */
double CELibGetLifeTimeofStarLSF(const double Mass, const double Metallicity){

    if(Metallicity <= CELibLifeTimeZ[0]){
        return pCELibGetLifeTimeofStarUsingLSF(0,Mass);
    } else if(Metallicity >= CELibLifeTimeZ[CELibLifeTime_Metallicity-1]){
        return pCELibGetLifeTimeofStarUsingLSF(CELibLifeTime_Metallicity-1,Mass);
    } else {
        int TableID = 0;
        for(int i=1;i<CELibLifeTime_Metallicity;i++){
            if(CELibLifeTimeZ[i] > Metallicity){
                TableID = i-1;
                break;
            }
        }
        double AgeLower = pCELibGetLifeTimeofStarUsingLSF(TableID,Mass);
        double AgeUpper = pCELibGetLifeTimeofStarUsingLSF(TableID+1,Mass);
        double grad = ((AgeUpper-AgeLower)/(CELibLifeTimeZ[TableID+1]-CELibLifeTimeZ[TableID]));
        return grad*(Metallicity-CELibLifeTimeZ[TableID])+AgeLower;
    }
}

/*!
 * This function writes lifetimes of stars for given masses based on data
 * obtained by the least square fitting. 
 */
static void pCELibStellarLifeTimeDumpInterpolatedValuesLSF(const char OutDir[]){

    const int NMetal = 10;
    const int NMass = 1000;
    MakeDir(OutDir);
    
    double Zmin = CELibLifeTimeLogZ[0];  
    double Zmax = CELibLifeTimeLogZ[CELibLifeTime_Metallicity-1];  
    double dMetal = (Zmax-Zmin)/NMetal;

    double Mmax = 1000;
    double Mmin = 0.1;
    double dMass = (Mmax-Mmin)/NMass;

    FILE *fp;
    char fname[MaxCharactersInLine];

    for(int i=0;i<NMetal;i++){
        double Metallicity = pow(10.0,dMetal*i+CELibLifeTimeLogZ[0]);

        Snprintf(fname,"%s/CELibLifeTimeLSF.%02d",OutDir,i);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#%g\n",Metallicity);
        fprintf(fp,"#Mass #Age\n");
        for(int k=0;k<NMass;k++){
            double Mass = dMass*k+Mmin;
            fprintf(fp,"%g %g\n",Mass,CELibGetLifeTimeofStarLSF(Mass,Metallicity));
        }
        fclose(fp);
    }

    return ;
}

/*!
 * This function returns an interpolated value of the mass whose lifetime is
 * just finished at a given age of a SSP particle. The parameter "LifeTime", which
 * represents the age of the SSP particle, should be in the unit of "year".
 */
static inline double pCELibGetDyingStellarMassUsingLSF(const int TableID, const double LifeTime){

    double LogLifeTime = log10(LifeTime);
    double LogMass=0.e0;
    for(int k=0;k<FitDim+1;k++){
        LogMass += FittingCoefMd[TableID][k]*pow(LogLifeTime,k);
    }
    return pow(10.0,LogMass);
}

/*!
 * This function returns the mass of the star of which life is just finished.
 * Unlike CELibGetDyingStellarMass, this function uses the functions obtained by
 * the least square fitting.
 * This function needs two input parameters. One is the age, and the other is
 * the metallicity. The parameter "LifeTime" should be given in the unit of "year",
 * while "metallicity" should be [0,1].
 */
double CELibGetDyingStellarMassLSF(const double LifeTime, const double Metallicity){

    if(Metallicity <= CELibLifeTimeZ[0]){
        return pCELibGetDyingStellarMassUsingLSF(0,LifeTime);
    } else if(Metallicity >= CELibLifeTimeZ[CELibLifeTime_Metallicity-1]){
        return pCELibGetDyingStellarMassUsingLSF(CELibLifeTime_Metallicity-1,LifeTime);
    } else {
        int TableID = 0;
        for(int i=1;i<CELibLifeTime_Metallicity;i++){
            if(CELibLifeTimeZ[i] > Metallicity){
                TableID = i-1;
                break;
            }
        }
        double MassLower = pCELibGetDyingStellarMassUsingLSF(TableID,LifeTime);
        double MassUpper = pCELibGetDyingStellarMassUsingLSF(TableID+1,LifeTime);
        double grad = ((MassUpper - MassLower)/(CELibLifeTimeZ[TableID+1]-CELibLifeTimeZ[TableID]));
        return grad*(Metallicity-CELibLifeTimeZ[TableID])+MassLower;
    }
}

/*!
 * This function returns next event time. 
 */
double CELibGetNextEventTimeStarbyStar(struct CELibStructNextEventTimeStarbyStarInput Input, const int Type){

    switch (Type){
        case CELibFeedbackType_SNII:
            {
                /*
                double Time;
                if(Input.noPopIII == 1){
                    Time = CELibGetLifeTimeofStar(Input.InitialMass_in_Msun,CELibLifeTimeZ[1]);
                } else {
                    Time = CELibGetLifeTimeofStar(Input.InitialMass_in_Msun,Input.Metallicity);
                }
                */
                double Time = CELibGetLifeTimeofStar(Input.InitialMass_in_Msun,Input.Metallicity);
                if(Time < 0.0){
                    fprintf(stderr,"Input rate is incorrect.\n");
                    return 0.e0;
                } else { 
                    return Time;
                }
            }
        default:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            return 0.e0;
    }
}