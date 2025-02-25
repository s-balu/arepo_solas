#include "config.h"
#include "NSMYields.h"

/*! \file NSMYields.c
 * \brief This file has functions regarding NSM yields.
 */

#define CELibNSMYieldsTableModelNumber 1

/*!
 * The yields table for NSM.
 */
struct CELibStructNSMYields *CELibNSMYields;
struct CELibStructNSMYields CELibNSMYields_W14[CELibNSMYieldsTableModelNumber] = {
	{
	.Name = "Wanajo+2014 model.",
    .Mej = 2e-5,
    .Energy = 0.e0,
	.Elements[CELibYield_H]  = 0.0,
	.Elements[CELibYield_He] = 0.0,
	.Elements[CELibYield_C]  = 0.0,
	.Elements[CELibYield_N]  = 0.0,
	.Elements[CELibYield_O]  = 0.0,
	.Elements[CELibYield_Ne] = 0.0,
	.Elements[CELibYield_Mg] = 0.0,
	.Elements[CELibYield_Si] = 0.0,
	.Elements[CELibYield_S]  = 0.0,
	.Elements[CELibYield_Ca] = 0.0,
	.Elements[CELibYield_Fe] = 0.0,
	.Elements[CELibYield_Ni] = 0.0,
	.Elements[CELibYield_Eu] = 2e-5, // 2.e-5 Msun per 1 NSM.
	},
};

/*!
 * Set yield names.
 */
static char NSMYieldsElementName[CELibYield_Number][MaxCharactersInLine];
static void pCELibNSMYieldsSetElementName(void){

    strcpy(NSMYieldsElementName[CELibYield_H],     "H"); 
    strcpy(NSMYieldsElementName[CELibYield_He],    "He"); 
    strcpy(NSMYieldsElementName[CELibYield_C],     "C"); 
    strcpy(NSMYieldsElementName[CELibYield_N],     "N"); 
    strcpy(NSMYieldsElementName[CELibYield_O],     "O"); 
    strcpy(NSMYieldsElementName[CELibYield_Ne],    "Ne"); 
    strcpy(NSMYieldsElementName[CELibYield_Mg],    "Mg"); 
    strcpy(NSMYieldsElementName[CELibYield_Si],    "Si"); 
    strcpy(NSMYieldsElementName[CELibYield_S],     "S"); 
    strcpy(NSMYieldsElementName[CELibYield_Ca],    "Ca"); 
    strcpy(NSMYieldsElementName[CELibYield_Fe],    "Fe"); 
    strcpy(NSMYieldsElementName[CELibYield_Ni],    "Ni"); 
    strcpy(NSMYieldsElementName[CELibYield_Eu],    "Eu"); 

    return ;
}


static void pCELibNSMYieldsDumpRate(const char OutDir[]);
static void pCELibNSMYieldsDumpYieldsTable(const char OutDir[]);


/*!
 * This function initializes the NSM yields table.
 */
void CELibInitNSMYields(void){

    CELibNSMYields = CELibNSMYields_W14;

    // Integrate the mass^-1 weighted IMF from NSMLowerMass to NSMUpperMass.
    double Normalization = IntegralSimpson(CELibRunParameters.NSMLowerMass,
                CELibRunParameters.NSMUpperMass,CELibRunParameters.IntegrationSteps,
                CELibIMF[CELibRunParameters.IMFType].IMFFunctionPerMass);
    CELibRunParameters.NSMNumberPerMass = Normalization*CELibRunParameters.NSMFraction;

    CELibRunParameters.NSMDTDNormalization = 1.0;
    CELibRunParameters.NSMDTDNormalization = CELibRunParameters.NSMNumberPerMass
                                             /CELibGetNSMIntegratedRatePowerLaw(1.e+10);

    if(CELibRunParameters.TestMode){
        fprintf(stderr,"//\t NSM parameters.\n");
        fprintf(stderr,"//\t\t NSM per mass = %g\n",CELibRunParameters.NSMNumberPerMass);
        fprintf(stderr,"//\t\t DTD normalization = %g\n",CELibRunParameters.NSMDTDNormalization);
        fprintf(stderr,"\n");
    }

    pCELibNSMYieldsSetElementName();

    if(CELibRunParameters.TestMode){
        pCELibNSMYieldsDumpRate("./CELib/CELibNSM");
        pCELibNSMYieldsDumpYieldsTable("./CELib/CELibNSM");
    }

    return ;
}

/*!
 * This function gives the adopted model data.
 */
void CELibGetNSMModelNames(char *ModelNames){

    snprintf(ModelNames,MaxCharactersInLine,
            "Type 0: Power law.\n");
            //"Type 0: Power law.\nType 1: Single event.\n");

    return ;
}

/*!
 * This function sets the current NSM model type into ModelName.
 */
void CELibGetNSMCurrentModelName(char *ModelName){
    if(CELibRunParameters.NSMType == CELibNSMRateModelID_PowerLaw){
        Snprintf(ModelName,"Type %d: Power law.",CELibRunParameters.NSMType);
    // } else if(CELibRunParameters.NSMType == 1){ // Not implemented
        // Snprintf(ModelName,"Type %d: Single event.",CELibRunParameters.NSMType);
    } else {
        Snprintf(ModelName,"Type %d: NULL",CELibRunParameters.NSMType);
    }
}

/*!
 * This function sets the yields table name into ModelName.
 */
void CELibGetNSMCurrentYieldModelName(char *ModelName){
    int TableID = CELibRunParameters.NSMYieldsTableID;
    Snprintf(ModelName,"%s\n",CELibNSMYields[TableID].Name);
    return ;
}


/*
 * Power law type DTD. The power law index is set to -1 and the time scale of
 * the power law function is 100 Myr (0.1Gyr) for the fiducial model. But should
 * be adjustable.
 */
#define pCELib_NSM_N0_FOR_DTDPowerLaw (4e-13)  // Normalization

/*!
 * This function returns the value of DTD at a given age.
 */
static double pCELibGetNSMDTDPowerLaw(const double Age_in_Year){
    return CELibRunParameters.NSMDTDNormalization*
        pow(Age_in_Year,CELibRunParameters.NSMDTDPowerLawIndex);
}

/*
 * This function returns the shortest time of NSM event in year.
 */
double CELibGetNSMEarliestExplosionTimePowerLaw(void){
    return CELibRunParameters.NSMDTDOffsetForPower;
}

static int pCELibDTDStep = 1000;

/*!
 * This function returns the integrated (cumulative) NSM DTD until the given
 * age.
 */
double CELibGetNSMIntegratedRatePowerLaw(const double Age_in_Year){
    if(Age_in_Year < CELibRunParameters.NSMDTDOffsetForPower){
        return 0.e0;
    } else {
        return IntegralSimpson(CELibRunParameters.NSMDTDOffsetForPower,Age_in_Year,
                pCELibDTDStep,&pCELibGetNSMDTDPowerLaw);
    }
}


/*!
 * This function writes the integrated (cumulative) NSM rates for various models and
 * metallicities.
 */
static void pCELibNSMYieldsDumpRate(const char OutDir[]){

    MakeDir(OutDir);

    int Steps = 400;
    double HubbleTime = 1.38e10; // in year

    double dAge = log10(HubbleTime)/Steps;

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/NSMrate.dat",OutDir);

    FileOpen(fp,fname,"w");
    for(int i=0;i<Steps;i++){
        double Age = pow(10.0,dAge*(i+0.5));

        double Mass = CELibGetDyingStellarMass(Age,0.e0);
        double RatePowerLaw = CELibGetNSMIntegratedRatePowerLaw(Age);
        fprintf(fp,"%g %g %g\n",Age,Mass,RatePowerLaw);
    }

    fclose(fp);

    return ;
}

/*!
 * Dump all yields tables.
 */
static void pCELibNSMYieldsDumpYieldsTable(const char OutDir[]){

    MakeDir(OutDir);

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/NSMYieldsTable.dat",OutDir);

    FileOpen(fp,fname,"w");

    for(int k=0;k<CELibNSMYieldsTableModelNumber;k++){
        fprintf(fp,"%s \n",CELibNSMYields->Name);
        for(int i=0;i<CELibYield_Number;i++){
            fprintf(fp,"%g #%s\n",CELibNSMYields->Elements[i],NSMYieldsElementName[i]);
        }
        fprintf(fp,"\n\n");
    }
    fclose(fp);

    return ;
}

/*!
 * This function returns masses of individual elements released from a NSM event.
 */
struct CELibStructFeedbackOutput CELibGetNSMFeedback(struct CELibStructFeedbackInput Input){

    struct CELibStructFeedbackOutput NSMFeedback;

    NSMFeedback.Energy = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        NSMFeedback.Elements[i] = CELibRunParameters.NSMNassociation*
            CELibNSMYields[CELibRunParameters.NSMYieldsTableID].Elements[i];
    }

    NSMFeedback.EjectaMass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        NSMFeedback.EjectaMass += NSMFeedback.Elements[i];
    }
    NSMFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor - NSMFeedback.EjectaMass;

    return NSMFeedback;
}

