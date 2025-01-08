#include "config.h"
#include "SNIaYields.h"
#include "../data/Iwamoto+1999/Iwamoto+1999_Struct.h"
#include "../data/Maeda+2010/Maeda+2010_Struct.h"
#include "../data/Seitenzahl+2013/Seitenzahl+2013_Struct.h"
#include "../data/Travaglio+2004/Travaglio+2004_Struct.h"

/*! \file SNIaYields.c
 * \brief This file has functions regarding SNIa yields.
 */

/*
 * Variables for the active yield table size.
 */
static int CELibSNIaYields_Metallicity; 

struct CELibStructSNIaYields *CELibSNIaYields;

static char SNIaYieldsElementName[CELibYield_Number][MaxCharactersInLine];
static void pCELibSNIaYieldsSetElementName(void){

    strcpy(SNIaYieldsElementName[CELibYield_H],    "H"); 
    strcpy(SNIaYieldsElementName[CELibYield_He],   "He"); 
    strcpy(SNIaYieldsElementName[CELibYield_C],    "C"); 
    strcpy(SNIaYieldsElementName[CELibYield_N],    "N"); 
    strcpy(SNIaYieldsElementName[CELibYield_O],    "O"); 
    strcpy(SNIaYieldsElementName[CELibYield_Ne],   "Ne"); 
    strcpy(SNIaYieldsElementName[CELibYield_Mg],   "Mg"); 
    strcpy(SNIaYieldsElementName[CELibYield_Si],   "Si"); 
    strcpy(SNIaYieldsElementName[CELibYield_S],    "S"); 
    strcpy(SNIaYieldsElementName[CELibYield_Ca],   "Ca"); 
    strcpy(SNIaYieldsElementName[CELibYield_Fe],   "Fe"); 
    strcpy(SNIaYieldsElementName[CELibYield_Ni],   "Ni"); 
    strcpy(SNIaYieldsElementName[CELibYield_Eu],   "Eu"); 

    return ;
}


static void pCELibDumpSNIaYieldsRate(const char OutDir[]);
static void pCELibDumpSNIaYieldsTables(const char OutDir[]);
static void pCELibDumpSNIaYieldsMainTable(const char OutDir[]);

/*
 * This function returns the total return mass of an element Index.
 */
static double pCELibGetSNIaYieldsMetalNoDependYield(const int Index, const double Metallicity){
    return CELibSNIaYields->Elements[Index]; 
}

/*
 * This function returns the metallicity dependent total return mass of an
 * element Index. This function uses the yield tables of N100, N100_Z0.5,
 * N100_Z0.1, N100_Z0.01. We assume that the metallicity ($^{22}$Ne) is the same
 * as the solar abundance ratio of Asplund et al. (2009).
 */
#define Z_solar (0.013)
static double pCELibGetSNIaYieldsMetalDependYield(const int Index, const double Metallicity){

    if(Metallicity <= 0.01*Z_solar){
        return CELibSNIaYieldsS13[16].Elements[Index]; 
    } else  if(Metallicity < 0.1*Z_solar){
        // Use 16 and 15
        double Grad = (CELibSNIaYieldsS13[15].Elements[Index]-CELibSNIaYieldsS13[16].Elements[Index])/(0.1*Z_solar-0.01*Z_solar);
        return Grad*(Metallicity-0.01*Z_solar)+CELibSNIaYieldsS13[16].Elements[Index];
    } else  if(Metallicity < 0.5*Z_solar){
        // Use 15 and 14
        double Grad = (CELibSNIaYieldsS13[14].Elements[Index]-CELibSNIaYieldsS13[15].Elements[Index])/(0.5*Z_solar-0.1*Z_solar);
        return Grad*(Metallicity-0.1*Z_solar)+CELibSNIaYieldsS13[15].Elements[Index];
    } else  if(Metallicity < Z_solar){
        // Use 14 and 7
        double Grad = (CELibSNIaYieldsS13[7].Elements[Index]-CELibSNIaYieldsS13[14].Elements[Index])/(Z_solar-0.5*Z_solar);
        return Grad*(Metallicity-0.5*Z_solar)+CELibSNIaYieldsS13[14].Elements[Index];
    } else {
    //} else if(Metallicity >= Z_solar){
        return CELibSNIaYieldsS13[7].Elements[Index]; 
    }
}

static double (*pCELibSNIaYieldsGetYield)(const int, const double);

void CELibInitSNIaYields(void){

    pCELibSNIaYieldsSetElementName();
    if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_I99){
        CELibSNIaYields = CELibSNIaYieldsI99+CELibRunParameters.SNIaYieldsTableModelID;
        if(CELibRunParameters.SNIaYieldsTableModelID >= CELibSNIaYieldsTableModelID_I99_Number){
            fprintf(stderr,"SNIa table model ID range over. %d %d\n",CELibRunParameters.SNIaYieldsTableID,CELibRunParameters.SNIaYieldsTableModelID);
            fprintf(stderr,"SNIa table model ID range over. 1\n");
            exit(1);
        }
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_M10){
        CELibSNIaYields = CELibSNIaYieldsM10+CELibRunParameters.SNIaYieldsTableModelID;
        if(CELibRunParameters.SNIaYieldsTableModelID >= CELibSNIaYieldsTableModelID_M10_Number){
            fprintf(stderr,"SNIa table model ID range over. 2\n");
            exit(1);
        }
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_S13){
        if(CELibRunParameters.SNIaYieldsTableModelID == 17){
            CELibSNIaYields = NULL;
        } else {
            CELibSNIaYields = CELibSNIaYieldsS13+CELibRunParameters.SNIaYieldsTableModelID;
            if(CELibRunParameters.SNIaYieldsTableModelID >= CELibSNIaYieldsTableModelID_S13_Number){
                fprintf(stderr,"SNIa table model ID range over. 3\n");
                exit(1);
            }
        }
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_T04){
        CELibSNIaYields = CELibSNIaYieldsT04+CELibRunParameters.SNIaYieldsTableModelID;
        if(CELibRunParameters.SNIaYieldsTableModelID >= CELibSNIaYieldsTableModelID_T04_Number){
            fprintf(stderr,"SNIa table model ID range over. 4\n");
            exit(1);
        }
    }

    if(CELibSNIaYields==NULL){
        pCELibSNIaYieldsGetYield = pCELibGetSNIaYieldsMetalDependYield;
    } else {
        pCELibSNIaYieldsGetYield = pCELibGetSNIaYieldsMetalNoDependYield;
    }

    if(CELibRunParameters.TestMode){
        MakeDir("./CELib");
        pCELibDumpSNIaYieldsRate("./CELib/CELibSNIa");
        pCELibDumpSNIaYieldsTables("./CELib/CELibSNIa");
        pCELibDumpSNIaYieldsMainTable("./CELib/CELibSNIa");
    }

    return ;
}


void CELibGetSNIaModelNames(char *ModelNames){

    snprintf(ModelNames,MaxCharactersInLine,
            "Type 0: Greggio & Renzini 1983\nType 1: Power law\nType 2: Power law (Vogelsberger+2013)\n");

    return ;
}

void CELibGetSNIaCurrentModelName(char *ModelName){

    if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_GR83){
        Snprintf(ModelName,"Type %d: Greggio & Renzini 1983",CELibRunParameters.SNIaType);
    } else if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_PowerLaw){
        Snprintf(ModelName,"Type %d: Power law",CELibRunParameters.SNIaType);
    } else if(CELibRunParameters.SNIaType == 2){
        Snprintf(ModelName,"Type %d: Power law (Vogelsberger+2013)",CELibRunParameters.SNIaType);
    } else {
        Snprintf(ModelName,"Type %d: NULL",CELibRunParameters.SNIaType);
    }
}

void CELibGetSNIaCurrentYieldModelName(char *ModelName){
    if(CELibSNIaYields == NULL){
        sprintf(ModelName,"Seitenzahl+2013, metal depend version.\n");
    } else {
        Snprintf(ModelName,"%s\n",CELibSNIaYields->Name);
    }
    return ;
}


static int pCELibDTDStep = 1000;
/*
 * DTD used in Portinari+1998 and it is originaly derived from Greggio & Renzini
 * 1983.
 */
static double f(const double mu){
    return 24*SQ(mu);
}

static double pCELibGetDTDGreggioRenzini83(const double MassB, const double IntLower){
    return IntegralSimpson(IntLower,0.5,100,&f);
}

static double pCELibGetSNIaRateIntegrandGreggioRenzini83(const double MassB, const double DyingStellarMass){
    const int CurrentIMFType = CELibRunParameters.IMFType;
    double IntLower = fmin(0.5,fmax(DyingStellarMass/MassB,(MassB-0.5*CELibRunParameters.SNIaUpperMassBinary)/MassB));
    return CELibIMF[CurrentIMFType].IMFFunctionPerMass(MassB)*pCELibGetDTDGreggioRenzini83(MassB,IntLower);
}

double CELibGetSNIaEarliestExplosionTimePortinari(const double Metallicity){
    return CELibGetLifeTimeofStar(CELibRunParameters.SNIaUpperMass,Metallicity);
}

/*!
 * The normalization for SNIa in the Portinari+ model is 0.05-0.08 in the case
 * M_{B,low} = 3 Msun. If one adopts M_{B,low} = 2 Msun, the quantity should be
 * 0.03-0.05. See Portinari et al. A&A, 334, 505, 1998.
 */
#define pCELib_SNIa_N0_FOR_Portinari (0.07)
double CELibGetSNIaIntegratedRateGreggioRenzini(const double Age, const double Metallicity){

    double DyingStellarMass = CELibGetDyingStellarMass(Age,Metallicity);

    if(DyingStellarMass == __CELib_NoValue__) 
        return 0.e0;

    return pCELib_SNIa_N0_FOR_Portinari*IntegralSimpsonOneExtraParameter(CELibRunParameters.SNIaLowerMass,
            CELibRunParameters.SNIaUpperMass,
            pCELibDTDStep,&pCELibGetSNIaRateIntegrandGreggioRenzini83,DyingStellarMass);
}

/*!
 * DTD used in Vogelsberger et al. (2013). Parameters are adopted from Maoz et
 * al. MNRAS, 426, 3282, 2012.
 */
#define pCELib_SNIa_N0_FOR_DTDVogelsberger (1.3e-3)
#define pCELib_SNIa_S_FOR_DTDVogelsberger (1.12)
#define pCELib_SNIa_Tau_FOR_DTDVogelsberger (4.e7)
/*!
 * This function returns the DTD of Vogelsberger et al. (2013).
 */
static double pCELibGetDTDVogelsberger(const double Age){
    return pCELib_SNIa_N0_FOR_DTDVogelsberger*pow(Age/pCELib_SNIa_Tau_FOR_DTDVogelsberger,-pCELib_SNIa_S_FOR_DTDVogelsberger)
        *((pCELib_SNIa_S_FOR_DTDVogelsberger-1.0)/(pCELib_SNIa_Tau_FOR_DTDVogelsberger));
}

/*!
 * This function returns the beginning time of SNIa for the Vogelsberger et al.
 * (2013) model. 
 */
double CELibGetSNIaEarliestExplosionTimeVogelsberger(void){
    return pCELib_SNIa_Tau_FOR_DTDVogelsberger;
}

/*!
 * This function returns the cumulative SNIa rate based on parameters of Vogelsberger et al. (2013).
 */
double CELibGetSNIaIntegratedRateVogelsberger(const double Age){
    if(Age < pCELib_SNIa_Tau_FOR_DTDVogelsberger){
        return 0.e0;
    } else {
        return IntegralSimpson(pCELib_SNIa_Tau_FOR_DTDVogelsberger,Age,pCELibDTDStep,&pCELibGetDTDVogelsberger);
    }
}

/*
 * Power law DTD with parameters of Maoz & Mannucci 2012. The integrated SNIa
 * number per Msun = 2.2e-3.  Maoz+ MNRAS, 426, 3283 (2012) provides almost
 * similar DTD with somewhat different parameters (S=1.12 and NSNIa/Msun =
 * 1.3e-3).
 */
#define pCELib_SNIa_N0_FOR_DTDMaozMannucci (4e-13)
#define pCELib_SNIa_S_FOR_DTDMaozMannucci (1.0)
#define pCELib_SNIa_Tau_FOR_DTDMaozMannucci (4.e7)
/*!
 * This function returns the DTD of Maoz and Mannucci (2012).
 */
static double pCELibGetDTDPowerLaw(const double Age){
    return pCELib_SNIa_N0_FOR_DTDMaozMannucci*pow(Age*1.e-9,-pCELib_SNIa_S_FOR_DTDMaozMannucci);
}

/*!
 * This function returns the beginning time of SNIa for the Maoz and Mannucci
 * (2102) model
 */
double CELibGetSNIaEarliestExplosionTimePowerLaw(void){
    return pCELib_SNIa_Tau_FOR_DTDMaozMannucci;
}

/*!
 * This function returns the cumulative SNIa rate based on parameters of Maoz &
 * Mannucci (2012).
 */
double CELibGetSNIaIntegratedRatePowerLaw(const double Age){
    if(Age < pCELib_SNIa_Tau_FOR_DTDMaozMannucci){
        return 0.e0;
    } else {
        return IntegralSimpson(pCELib_SNIa_Tau_FOR_DTDMaozMannucci,Age,pCELibDTDStep,&pCELibGetDTDPowerLaw);
    }
}

/*!
 * This function writes integrated SNIa rates for various models and
 * metallicities.
 *
 * One can find that the first time derivatives of these integrated rates have
 * saw-toothed shape. This is because the age table of CELib is constructed by the
 * linear interpolation. To obtain the smoothed rate in the first derivative,
 * one might adopt a smoothed interpolation, say spline.
 */
static void pCELibDumpSNIaYieldsRate(const char OutDir[]){

    MakeDir(OutDir);

    int Steps = 400;
    double HubbleTime = 1.38e10; // in year

    double dAge = log10(HubbleTime)/Steps;

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/SNIaRate.dat",OutDir);

    FileOpen(fp,fname,"w");
    for(int i=0;i<Steps;i++){
        double Age = pow(10.0,dAge*(i+0.5));

        double Mass = CELibGetDyingStellarMass(Age,0.e0);
        double RateMM = CELibGetSNIaIntegratedRatePowerLaw(Age);
        double RateVog = CELibGetSNIaIntegratedRateVogelsberger(Age);
        double RatePortinari0 = CELibGetSNIaIntegratedRateGreggioRenzini(Age,0.0);
        double RatePortinari1 = CELibGetSNIaIntegratedRateGreggioRenzini(Age,0.02);
        double RatePortinari2 = CELibGetSNIaIntegratedRateGreggioRenzini(Age,0.05);
        fprintf(fp,"%g %g %g %g %g %g %g %g %g\n",Age,Mass,RateMM,RateVog,
                RatePortinari0,RatePortinari1,RatePortinari2,
                CELibGetSNIaIntegratedRateGreggioRenzini(pow(10.0,dAge*(i+1)),0.0)-
                CELibGetSNIaIntegratedRateGreggioRenzini(pow(10.0,dAge*(i)),0.0),
                CELibGetSNIaIntegratedRateGreggioRenzini(pow(10.0,dAge*(i+1)),0.02)-
                CELibGetSNIaIntegratedRateGreggioRenzini(pow(10.0,dAge*(i)),0.02));
    }

    fclose(fp);

    return ;
}

/*!
 * This function writes yields table used for SNIa. Under construction.
 */
static void pCELibDumpSNIaYieldsTables(const char OutDir[]){

    MakeDir(OutDir);

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/SNIaYieldsTable.dat",OutDir);

    FileOpen(fp,fname,"w");

    for(int k=0;k<CELibSNIaYieldsTableModelID_I99_Number;k++){
        fprintf(fp,"%s \n",CELibSNIaYieldsI99[k].Name);
        for(int i=0;i<CELibYield_Number;i++){
            fprintf(fp,"%g #%s\n",CELibSNIaYieldsI99[k].Elements[i],SNIaYieldsElementName[i]);
        }
        fprintf(fp,"\n\n");
    }
    fclose(fp);

    return ;
}

/*!
 * This function writes the adopted yields table of SNIa in a file.
 */
static void pCELibDumpSNIaYieldsMainTable(const char OutDir[]){

    MakeDir(OutDir);

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/SNIaYieldMainTable.dat",OutDir);

    FileOpen(fp,fname,"w");
    char ModelName[MaxCharactersInLine];
    CELibGetSNIaCurrentYieldModelName(ModelName);

    fprintf(fp,"#%s\n",ModelName);

    if(CELibSNIaYields==NULL){
        double Z[] = {0.0,0.0001,0.0005,0.001,0.005,0.01,0.02};
        int NMetal = sizeof(Z)/sizeof(double);
        for(int k=0;k<NMetal;k++){
            fprintf(fp,"#Metallicity = %g\n",Z[k]);
            for(int i=0;i<CELibYield_Number;i++){
                fprintf(fp,"%g %s\n",pCELibSNIaYieldsGetYield(i,Z[k]),SNIaYieldsElementName[i]);
            }
            fprintf(fp,"\n");
        }
    } else {
        for(int i=0;i<CELibYield_Number;i++){
            fprintf(fp,"%g %s\n",CELibSNIaYields->Elements[i],SNIaYieldsElementName[i]);
            fprintf(fp,"%g %s\n",pCELibSNIaYieldsGetYield(i,0.0),SNIaYieldsElementName[i]);
        }
    }

    fclose(fp);

    return ;
}

/*!
 * This function returns energy, masses of individual elements and the sum of them
 * (ejecta mass) released from SNeIa.  It also returns the remnant mass.
 */
struct CELibStructFeedbackOutput CELibGetSNIaFeedback(struct CELibStructFeedbackInput Input){

    struct CELibStructFeedbackOutput SNIaFeedback;

    SNIaFeedback.Energy = CELibRunParameters.SNIaNassociation*CELibRunParameters.SNIaEnergy;
    for(int i=0;i<CELibYield_Number;i++){
        SNIaFeedback.Elements[i] = CELibRunParameters.SNIaNassociation*pCELibSNIaYieldsGetYield(i,Input.Metallicity);
    }

    SNIaFeedback.EjectaMass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        SNIaFeedback.EjectaMass += SNIaFeedback.Elements[i];
    }
    SNIaFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor - SNIaFeedback.EjectaMass;

    return SNIaFeedback;
}



