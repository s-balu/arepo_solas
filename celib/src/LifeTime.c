#include "config.h"
#include "LifeTimeTable.h"


// Define the function pointer
double (*CELibGetDyingStellarMass)(const double, const double) = NULL;
double (*CELibGetLifeTimeofStar)(const double, const double) = NULL;

/*!
 * Sizes of the lifetime table. 
 */
int CELibLifeTime_Metallicity; 
int CELibLifeTime_Mass; 

double *CELibLifeTimeZ;
double *CELibLifeTimeLogZ;

static void pCELibStellarLifeTimeDumpInterpolatedValues(const char OutDir[]);
static void pCELibDyingStarMassDump(const char OutDir[]);
static void pCELibDyingStarMassDumpInterpolatedValues(const char OutDir[]);
static void pCELibFittingLifeTime(void);
static void pCELibStellarLifeTimeDumpInterpolatedValuesLSF(const char OutDir[]);
static void pCELibDyingStarMassDumpInterpolatedValuesLSF(const char OutDir[]);
static void pCELibDyingStarMassDumpLSF(const char OutDir[]);
static void pCELibShowZeroZStarLifeTime(void);


int CELibGetLifeTimeTableSizeMetallicity(void){
    return CELibLifeTime_Metallicity;
}

/*!
 * This function initializes all data concerning the stellar lifetimes.
 */
void CELibInitLifeTime(void){

    if(CELibRunParameters.PopIIILifeTime == 1){
        // Adopt Schaerer's popIII lifetime data 
        CELibLifeTime_Metallicity = CELIB_LIFETIME_Z_P98+1; 
        free(CELibLifeTimeZ);
        free(CELibLifeTimeLogZ);
        CELibLifeTimeZ = malloc(sizeof(double)*CELibLifeTime_Metallicity);
        CELibLifeTimeLogZ = malloc(sizeof(double)*CELibLifeTime_Metallicity);

        CELibLifeTime_Mass = CELIB_LIFETIME_M_P98;

        CELibLifeTimeZ[0] = CELibRunParameters.PopIIIMetallicity;
        CELibLifeTimeLogZ[0] = log10(CELibRunParameters.PopIIIMetallicity);
        for(int i=0;i<CELIB_LIFETIME_Z_P98;i++){
            CELibLifeTimeZ[i+1] = CELibLifeTimeZ_P98[i];
            CELibLifeTimeLogZ[i+1] = log10(CELibLifeTimeZ_P98[i]);
        }
    } else {
        // Use only Portinari et al. (1998)'s lifetime data.
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

    // Write log.
    if(CELibRunParameters.TestMode){
        pCELibDyingStarMassDump("./CELib/CELibLife");
        pCELibDyingStarMassDumpInterpolatedValues("./CELib/CELibLife");
        pCELibStellarLifeTimeDumpInterpolatedValues("./CELib/CELibLife");

        pCELibDyingStarMassDumpLSF("./CELib/CELibLife");
        pCELibDyingStarMassDumpInterpolatedValuesLSF("./CELib/CELibLife");
        pCELibStellarLifeTimeDumpInterpolatedValuesLSF("./CELib/CELibLife");
    }

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

/*!
 * This function writes the lifetime table of Portinari et al. (1998) in files.
 */
static void pCELibStellarLifeTimeDumpInterpolatedValues(const char OutDir[]){

    const int NMetal = 10;
    const int NMass = 100;
    MakeDir(OutDir);
    
    double dMetal = (log10(CELibLifeTimeZ_P98[CELIB_LIFETIME_Z_P98-1]) - log10(CELibLifeTimeZ_P98[0]))/NMetal;
    double dMass = (CELibLifeTimeMass_P98[CELIB_LIFETIME_M_P98-1] - CELibLifeTimeMass_P98[0])/NMass;

    FILE *fp;
    char fname[MaxCharactersInLine];

    for(int i=0;i<NMetal;i++){
        double Metallicity = pow(10.0,dMetal*i+log10(CELibLifeTimeZ_P98[0]));

        Snprintf(fname,"%s/CELibLifeTime.%02d",OutDir,i);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#%g\n",Metallicity);
        fprintf(fp,"#Mass #Age\n");
        for(int k=0;k<NMass;k++){
            double Mass = dMass*k+CELibLifeTimeMass_P98[0];
            fprintf(fp,"%g %g\n",Mass,CELibGetLifeTimeofStar(Mass,Metallicity));
        }
        fclose(fp);
    }

    return ;
}

/*!
 * This function writes ages and stellar masses who finish their lives.
 * Data are on the grid points of the lifetime table.
 */
static void pCELibDyingStarMassDump(const char OutDir[]){

    MakeDir(OutDir);
    
    FILE *fp;
    char fname[MaxCharactersInLine];

    const int NAge = 100;
    double dAge = (log10(1.5e10)-log10(1.e6))/NAge;

    for(int i=0;i<CELIB_LIFETIME_Z_P98;i++){

        Snprintf(fname,"%s/CELibDyingMass.%02d",OutDir,i);
        FileOpen(fp,fname,"w");

        for(int k=0;k<NAge;k++){
            double Age = pow(10.0,dAge*k+log10(1.e6));
            fprintf(fp,"%g %g\n",Age,pCELibGetDyingStellarMassFromTable(i,Age));
        }
        fclose(fp);
    }

    return ;
}

/*!
 * This function writes ages and stellar masses who finish their lives.  Unlike
 * pCELibDyingStarMassDump(const char OutDir[]), the interpolation is applied in
 * the direction of metallicity.
 */
static void pCELibDyingStarMassDumpInterpolatedValues(const char OutDir[]){

    MakeDir(OutDir);

    const int NMetal = 10;
    double dMetal = (log10(CELibLifeTimeZ[CELibLifeTime_Metallicity-1]) - log10(CELibLifeTimeZ[0]))/NMetal;
    
    FILE *fp;
    char fname[MaxCharactersInLine];

    const int NAge = 10000;
    double dAge = (log10(1.5e10)-log10(1.e6))/NAge;

    for(int i=0;i<NMetal;i++){

        double Metallicity = pow(10.0,dMetal*i+log10(CELibLifeTimeZ_P98[0]));

        Snprintf(fname,"%s/CELibDyingMassInterpolated.%02d",OutDir,i);
        FileOpen(fp,fname,"w");

        for(int k=0;k<NAge;k++){
            double Age = pow(10.0,dAge*k+log10(1.e6));
            fprintf(fp,"%g %g\n",Age,CELibGetDyingStellarMassTable(Age,Metallicity));
        }
        fclose(fp);
    }

    return ;
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
 * This function shows the coefficients of the fitted stellar lifetimes.
 */
static void pCELibLifeTimeShowFittingCoefficients(void){

    fprintf(stderr,"//\t Show coefficients of lifetime fittings\n");
    for(int i=0;i<CELibLifeTime_Metallicity;i++){
        fprintf(stderr,"//\t\t Z = %g\n",CELibLifeTimeZ[i]);
        for(int k=0;k<FitDim+1;k++){
            fprintf(stderr,"//\t\t LT Coef[%d][%d] = %g\n",i,k,FittingCoef[i][k]);
        }
        for(int k=0;k<FitDim;k++){
            fprintf(stderr,"//\t\t MD Coef[%d][%d] = %g\n",i,k,FittingCoefMd[i][k]);
        }
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
        if(CELibRunParameters.PopIIILifeTime == 1){
            if(i == 0){
                for(int k=0;k<CELIB_LIFETIME_M_P98;k++){
                    tmpLT[k] = CELibLifeTimeZMF_P98[i][k];
                }
                for(int k=0;k<CELIB_LIFETIME_M_S02;k++){
                    tmpLT[CELIB_LIFETIME_M_P98+k] = CELibLifeTimeZMF_S02[k];
                }
                Length = CELIB_LIFETIME_M_P98+CELIB_LIFETIME_M_S02;
            } else {
                for(int k=0;k<CELibLifeTime_Mass;k++){
                    tmpLT[k] = CELibLifeTimeZMF_P98[i-1][k];
                }
                Length = CELibLifeTime_Mass;
            }
        } else {
            for(int k=0;k<CELibLifeTime_Mass;k++){
                tmpLT[k] = CELibLifeTimeZMF_P98[i][k];
            }
            Length = CELibLifeTime_Mass;
        }
        pCELibFittingLifeTimei(tmpMass,tmpLT,Length,FittingCoef[i],FittingCoefMd[i]);
    }

    if(CELibRunParameters.TestMode == true)
        pCELibLifeTimeShowFittingCoefficients();

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
 * This function writes the mass of stars whose lifetime is just finished.  The
 * data with the least-square fitting are used. The interpolation is used for
 * the age.
 */
static void pCELibDyingStarMassDumpLSF(const char OutDir[]){

    MakeDir(OutDir);
    
    FILE *fp;
    char fname[MaxCharactersInLine];

    const int NLifeTime = 100;
    double dLifeTime = (log10(1.5e10)-log10(1.e6))/NLifeTime;

    for(int i=0;i<CELibLifeTime_Metallicity;i++){

        Snprintf(fname,"%s/CELibDyingMassLSF.%02d",OutDir,i);
        FileOpen(fp,fname,"w");

        for(int k=0;k<NLifeTime;k++){
            double LifeTime = pow(10.0,dLifeTime*k+log10(1.e6));
            fprintf(fp,"%g %g\n",LifeTime,pCELibGetDyingStellarMassUsingLSF(i,LifeTime));
        }
        fclose(fp);
    }

    return ;
}

/*!
 * This function writes the mass of stars whose lifetime is just finished.  The
 * data with the least-square fitting are used. The interpolation is used for
 * both age and metallicity.
 */
static void pCELibDyingStarMassDumpInterpolatedValuesLSF(const char OutDir[]){

    MakeDir(OutDir);

    const int NMetal = 10;
    double Zmin = CELibLifeTimeLogZ[0];  
    double Zmax = CELibLifeTimeLogZ[CELibLifeTime_Metallicity-1];  
    double dMetal = (Zmax-Zmin)/NMetal;

    FILE *fp;
    char fname[MaxCharactersInLine];

    const int NAge = 1000;
    double dAge = (log10(1.5e10)-log10(1.e6))/NAge;

    for(int i=0;i<NMetal;i++){
        double Metallicity = pow(10.0,dMetal*i+CELibLifeTimeLogZ[0]);

        Snprintf(fname,"%s/CELibDyingMassInterpolatedLSF.%02d",OutDir,i);
        FileOpen(fp,fname,"w");

        for(int k=0;k<NAge;k++){
            double Age = pow(10.0,dAge*k+log10(1.e6));
            fprintf(fp,"%g %g\n",Age,CELibGetDyingStellarMassLSF(Age,Metallicity));
        }
        fclose(fp);
    }

    return ;
}


/*!
 * This function shows the Schaerer (2002)'s zero-metal stars lifetime.
 */
static void pCELibShowZeroZStarLifeTime(void){

    double Mass[] = {100,150,200,220,300,500,1000};

    fprintf(stderr,"//\t Show lifetime of zero metal stars.\n");
    for(int i=0;i<7;i++){
        double LogMass = log10(Mass[i]);
        double LogAge=0.e0;
        for(int k=0;k<FitDim+1;k++){
            LogAge += FittingCoef[0][k]*pow(LogMass,k);
        }
        fprintf(stderr,"//\t For %g Msun stars, lifetime = %g yr.\n",Mass[i],pow(10.0,LogAge));
    }

    return ;
}
