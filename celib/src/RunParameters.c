#include "config.h"
#include "Init.h"
#include "RunParameters.h"
#include "IMFFunctions.h"
#include "SNIIYields.h"
#include "SNIaYields.h"

/*! \file 
 * \brief The CELib control structure and functions which change run parameters
 * are defined in this file.
 */


/*!
 * This sturcture holds all control parameters of CELib. Initialization of the
 * yield tables is based on parameters given in this structure. When parameters in
 * this structure are changed, a user code should call the initialization
 * routine of CELib.
 */
struct CELibStructRunParameters CELibRunParameters = {
    .TestMode = false,
    .IntegrationSteps = 10000,
    .IMFType = CELibIMF_Salpeter,
    .SolarAbundancePatternType = CELibSolarAbundancePattern_A09,

    .LifeTimeType = 1,

    /* Parameters for SNeII */
    .SNIIYieldsTableID = CELibSNIIYieldsTableID_N13,
    .SNIIYieldsModificationP98 = 1,
    .SNIIHyperNovaFraction = 0.0,
    .SNIIEnergy = 1.e+51,
    .SNIIEnergyEfficiency = 1.0,
    .SNIIUpperMass = 120.0,
    .SNIILowerMass = 8.0,
    .SNIIExplosionTimeFollowYields = 1,

    /* Parameters for SNeIa */
    .SNIaType = CELibSNIaRateModelID_PowerLaw, 
    .SNIaYieldsTableID = CELibSNIaYieldsTableID_S13, 
    .SNIaYieldsTableModelID = CELibSNIaYieldsTableModelID_S13_N100, 
    .SNIaEnergy = 1.e+51,
    .SNIaEnergyEfficiency = 1.0,
    .SNIaUpperMassBinary = 12.0,
    .SNIaLowerMassBinary = 3.0,
    .SNIaNassociation = 3,
    .SNIaUpperMass = 6.0,
    .SNIaLowerMass = 3.0,

    /* Parameters for AGB mass loss */
    .AGBYieldsTableID = CELibAGBYieldsTableID_K10D14,
    .AGBBinType = CELibAGBRateModelID_LinearBin,
    .AGBBinNumber = 100,
    .AGBBinTimeInterval = 1.e8,
    .AGBBinUpperAge = 1.38e10,
    .AGBBinLowerAge = 3.0e7,
    .AGBUpperMass = 8.0,
    .AGBLowerMass = 1.0,

    /* Parameters for NSM */
    .NSMYieldsTableID = CELibNSMYieldsTableID_W14, 
    .NSMType = CELibNSMRateModelID_PowerLaw,
    .NSMNumberPerMass = 0.0,
    .NSMUpperMass = 20.0,
    .NSMLowerMass = 8.0,
    .NSMFraction = 0.002,
    .NSMDTDNormalization = 1.0,
    .NSMDTDPowerLawIndex = -1.0,
    .NSMDTDOffsetForPower = 1.e+8,
    .NSMNassociation = 1,  // The size of NSM association.
    .NSMEnergy = 0.e0,
    .NSMEnergyPerMass = 0.e0,
    .NSMEnergyEfficiency = 0.e0,

    /* Parameters for PopIII */
    .PopIIIIMF = 0,
    .PopIIISNe = 0,
    .PopIIIAGB = 0, 
    .PopIIIAGBYieldsTableID = CELibAGBZ0YieldsTableID_CL08G13, 
    .PopIIIMetallicity = 1.34e-7, //1.e-5 x Z_sun,
    .PopIIILifeTime = 0,
};


static struct CELibStructRunParameters CELibRunParametersBackup; //!<
/*
 * This function saves the current paramter set of CELib to
 * ``CELibRunParametersBackup''.
 */
void CELibSaveRunParameter(void){
    CELibRunParametersBackup = CELibRunParameters;
    return ;
}

/*
 * This function loads the parameters saved in ``CELibRunParametersBackup''. 
 */
void CELibLoadRunParameter(void){
    CELibRunParameters = CELibRunParametersBackup;
    CELibInit();
    return ;
}

/*
 * This function set the test mode flag.
 * If ``true'', CELib writes logs frequently. 
 */
void CELibSetRunParameterTestMode(const bool TestMode){
    CELibRunParameters.TestMode = TestMode;
    return ;
}

/*
 * This function set IntegrationSteps.
 */
void CELibSetRunParameterIntegrationSteps(const int IntegrationSteps){
    CELibRunParameters.IntegrationSteps = IntegrationSteps;
    return ;
}

/*
 * This function set the IMF type.
 */
void CELibSetRunParameterIMFType(const int IMFType){
    CELibRunParameters.IMFType = IMFType;
    return ;
}

/*
 * This function set the type of solar abundance pattern.
 */
void CELibSetRunParameterSolarAbundancePatternType(const int SolarAbundancePatternType){
    CELibRunParameters.SolarAbundancePatternType = SolarAbundancePatternType;
    return ;
}

/*
 * This function set the lifetime type.
 */
void CELibSetRunParameterLifeTimeType(const int LifeTimeType){
    CELibRunParameters.LifeTimeType = LifeTimeType;
    return ;
}

void CELibSetRunParameterSNIIYieldsTableID(const int SNIIYieldsTableID){
    if(SNIIYieldsTableID > 1){
        fprintf(stderr,"SNIIYieldsTableID should be < 2\n");
        abort();
    }

    CELibRunParameters.SNIIYieldsTableID = SNIIYieldsTableID;
    return ;
}

void CELibSetRunParameterSNIIYieldsModificationP98(const int SNIIYieldsModificationP98){
    CELibRunParameters.SNIIYieldsModificationP98 = SNIIYieldsModificationP98;
    return ;
}

void CELibSetRunParameterSNIIRange(const double SNIIUpperMass_in_SolarMass, const double SNIILowerMass_in_SolarMass){
    if(SNIIUpperMass_in_SolarMass <= SNIILowerMass_in_SolarMass){
        fprintf(stderr,"Range error! %s:%s:%d\n",__FUNCTION__,__FILE__,__LINE__);
    }
    CELibRunParameters.SNIIUpperMass = SNIIUpperMass_in_SolarMass;
    CELibRunParameters.SNIILowerMass = SNIILowerMass_in_SolarMass;
    return ;
}

void CELibSetRunParameterSNIIHyperNovaFraction(const double SNIIHyperNovaFraction){
    CELibRunParameters.SNIIHyperNovaFraction = SNIIHyperNovaFraction;
    return ;
}

void CELibSetRunParameterSNIIExplosionTimeFollowYields(const int SNIIExplosionTimeFollowYields){
    CELibRunParameters.SNIIExplosionTimeFollowYields = SNIIExplosionTimeFollowYields;
    return ;
}

void CELibSetRunParameterSNIaType(const int SNIaType){
    CELibRunParameters.SNIaType = SNIaType;
    return ;
}

void CELibSetRunParameterSNIaYieldsTableID(const int SNIaYieldsTableID){
    CELibRunParameters.SNIaYieldsTableID = SNIaYieldsTableID;
    return ;
}

void CELibSetRunParameterSNIaYieldsTableModelID(const int SNIaYieldsTableModelID){
    CELibRunParameters.SNIaYieldsTableModelID = SNIaYieldsTableModelID;
    return ;
}

void CELibSetRunParameterSNIaRange(const double SNIaUpperMass_in_SolarMass, const double SNIaLowerMass_in_SolarMass){
    CELibRunParameters.SNIaUpperMass = SNIaUpperMass_in_SolarMass;
    CELibRunParameters.SNIaLowerMass = SNIaLowerMass_in_SolarMass;
    return ;
}

void CELibSetRunParameterSNIaUpperMassBinary(const double SNIaUpperMassBinary){
    CELibRunParameters.SNIaUpperMassBinary = SNIaUpperMassBinary;
    return ;
}

void CELibSetRunParameterSNIaLowerMassBinary(const double SNIaLowerMassBinary){
    CELibRunParameters.SNIaLowerMassBinary = SNIaLowerMassBinary;
    return ;
}

void CELibSetRunParameterSNIaNassociation(const int SNIaNassociation){
    CELibRunParameters.SNIaNassociation = SNIaNassociation;
    return ;
}

void CELibSetRunParameterAGBYieldsTableID(const int AGBYieldsTableID){
    CELibRunParameters.AGBYieldsTableID = AGBYieldsTableID;
    return ;
}

void CELibSetRunParameterAGBRange(const double AGBUpperMass_in_SolarMass, const double AGBLowerMass_in_SolarMass){
    if(AGBUpperMass_in_SolarMass <= AGBLowerMass_in_SolarMass){
        fprintf(stderr,"Range error! %s:%s:%d\n",__FUNCTION__,__FILE__,__LINE__);
    }
    CELibRunParameters.AGBUpperMass = AGBUpperMass_in_SolarMass;
    CELibRunParameters.AGBLowerMass = AGBLowerMass_in_SolarMass;
    return ;
}

void CELibSetRunParameterAGBBinUpperAge(const int AGBBinUpperAge){
    CELibRunParameters.AGBBinUpperAge = AGBBinUpperAge;
    return ;
}

void CELibSetRunParameterAGBBinLowerAge(const int AGBBinLowerAge){
    CELibRunParameters.AGBBinLowerAge = AGBBinLowerAge;
    return ;
}

void CELibSetRunParameterAGBBinNumber(const int AGBBinNumber){
    CELibRunParameters.AGBBinNumber = AGBBinNumber;
    return ;
}

void CELibSetRunParameterAGBBinType(const int AGBBinType){
    CELibRunParameters.AGBBinType = AGBBinType;
    return ;
}

void CELibSetRunParameterAGBBinTimeInterval(const double AGBBinTimeInterval){
    CELibRunParameters.AGBBinTimeInterval = AGBBinTimeInterval;
    return ;
}

void CELibSetRunParameterNSMYieldsTableID(const int NSMYieldsTableID){
    CELibRunParameters.NSMYieldsTableID = NSMYieldsTableID;
    return ;
}

void CELibSetRunParameterNSMType(const int NSMType){
    CELibRunParameters.NSMType = NSMType;
    return ;
}

void CELibSetRunParameterNSMDTDPowerLawIndex(const double NSMDTDPowerLawIndex){
    CELibRunParameters.NSMDTDPowerLawIndex = NSMDTDPowerLawIndex;
    return ;
}

void CELibSetRunParameterNSMDTDOffsetForPower(const double NSMDTDOffsetForPower){
    CELibRunParameters.NSMDTDOffsetForPower = NSMDTDOffsetForPower;
    return ;
}

void CELibSetRunParameterNSMNassociation(const int NSMNassociation){
    CELibRunParameters.NSMNassociation = NSMNassociation;
    return ;
}

void CELibSetRunParameterPopIIIIMF(const int PopIIIIMF){
    CELibRunParameters.PopIIIIMF = PopIIIIMF;
    return ;
}

void CELibSetRunParameterPopIIISNe(const int PopIIISNe){
    CELibRunParameters.PopIIISNe = PopIIISNe;
    return ;
}

void CELibSetRunParameterPopIIIAGB(const int PopIIIAGB){
    CELibRunParameters.PopIIIAGB = PopIIIAGB;
    return ;
}

void CELibSetRunParameterPopIIIAGBYieldsTableID(const int PopIIIAGBYieldsTableID){
    CELibRunParameters.PopIIIAGBYieldsTableID = PopIIIAGBYieldsTableID;
    return ;
}

void CELibSetRunParameterPopIIIMetallicity(const double PopIIIMetallicity){
    CELibRunParameters.PopIIIMetallicity = PopIIIMetallicity;
    return ;
}

void CELibSetRunParameterPopIIILifeTime(const int PopIIILifeTime){
    CELibRunParameters.PopIIILifeTime = PopIIILifeTime;
    return ;
}

void CELibShowCurrentRunParameters(void){

    fprintf(stderr,"//\t CELib Run Parameters.\n");
    fprintf(stderr,"//\t\t Test Mode = %d\n",CELibRunParameters.TestMode);
    fprintf(stderr,"//\t\t IMF type = %d\n",CELibRunParameters.IMFType);
    char IMFName[MaxCharactersInLine];
    CELibGetIMFName(IMFName);
    fprintf(stderr,"//\t\t %s\n",IMFName);

    fprintf(stderr,"//\n");
    fprintf(stderr,"//\t\t SNII mass range = %g -- %g Msun\n",CELibRunParameters.SNIILowerMass,CELibRunParameters.SNIIUpperMass);
    //fprintf(stderr," # of Type II SNe per Msun = %g\n",CELibRunParameters.SNIIEnergyEfficiency);
    fprintf(stderr,"//\t\t # of SNe II  per Msun = %g\n",CELibRunParameters.SNIINumberPerMass);
    fprintf(stderr,"//\t\t SN II yields table ID = %d\n",CELibRunParameters.SNIIYieldsTableID);
    fprintf(stderr,"//\t\t Hypernovae blend fraction = %g\n",CELibRunParameters.SNIIHyperNovaFraction);

    fprintf(stderr,"//\n");
    fprintf(stderr,"//\t\t SNIa type = %d\n",CELibRunParameters.SNIaType);
    char ModelName[MaxCharactersInLine];
    CELibGetSNIaCurrentModelName(ModelName);
    fprintf(stderr,"//\t\t %s\n",ModelName);
    if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_GR83)
        fprintf(stderr,"//\t\t SNIa mass range = %g -- %g Msun\n",CELibRunParameters.SNIaLowerMass,CELibRunParameters.SNIaUpperMass);
    fprintf(stderr,"//\t\t SNIa event number is %d\n",CELibRunParameters.SNIaNassociation);
    fprintf(stderr,"//\t\t SN Ia yields table ID = %d\n",CELibRunParameters.SNIaYieldsTableID);
    fprintf(stderr,"//\t\t SN Ia yields table model ID = %d\n",CELibRunParameters.SNIaYieldsTableModelID);


    fprintf(stderr,"//\n");
    fprintf(stderr,"//\t\t AGB yields table ID = %d\n",CELibRunParameters.AGBYieldsTableID);
    fprintf(stderr,"//\t\t AGB mass = %g -- %g Msun\n",CELibRunParameters.AGBLowerMass,CELibRunParameters.AGBUpperMass);
    fprintf(stderr,"//\t\t AGB bin number = %d\n",CELibRunParameters.AGBBinNumber);
    fprintf(stderr,"//\t\t AGB bin type is ");
    if(CELibRunParameters.AGBBinType){
        fprintf(stderr," log\n");
    } else {
        fprintf(stderr," linear\n");
    }

    fprintf(stderr,"//\n");
    fprintf(stderr,"//\t\t NSM type = %d\n",CELibRunParameters.NSMType);
    fprintf(stderr,"//\t\t NSM yields table ID = %d\n",CELibRunParameters.NSMYieldsTableID);
    fprintf(stderr,"//\t\t NSMs mass range = %g -- %g Msun\n",CELibRunParameters.NSMLowerMass,CELibRunParameters.NSMUpperMass);
    fprintf(stderr,"//\t\t NSMs fraction = %g\n",CELibRunParameters.NSMFraction);
    fprintf(stderr,"//\t\t NSMs DTD power-law index = %g\n",CELibRunParameters.NSMDTDPowerLawIndex);
    fprintf(stderr,"//\t\t NSMs DTD offset = %g [yr]\n",CELibRunParameters.NSMDTDOffsetForPower);

    fprintf(stderr,"//\n");
    fprintf(stderr,"//\t\t Pop III IMF = %d\n",CELibRunParameters.PopIIIIMF);
    fprintf(stderr,"//\t\t Pop III Metallicity = %g\n",CELibRunParameters.PopIIIMetallicity);
    fprintf(stderr,"//\t\t Pop III Lifetime = %d\n",CELibRunParameters.PopIIILifeTime);
    fprintf(stderr,"//\t\t Pop III SNe = %d\n",CELibRunParameters.PopIIISNe);
    fprintf(stderr,"//\t\t Pop III AGB = %d\n",CELibRunParameters.PopIIIAGB);
    fprintf(stderr,"//\t\t Pop III AGB yields table ID = %d\n",CELibRunParameters.PopIIIAGBYieldsTableID);


    return ;
}

bool CELibGetRunParameterTestMode(void){
    return CELibRunParameters.TestMode;
}

int CELibGetRunParameterIntegrationSteps(void){
    return CELibRunParameters.IntegrationSteps;
}

int CELibGetRunParameterIMFType(void){
    return CELibRunParameters.IMFType;
}

void CELibGetRunParameterIMFName(char *IMFName){
    Snprintf(IMFName,"%s",CELibIMF[CELibRunParameters.IMFType].Name);
    return;
}

int CELibGetRunParameterSolarAbundancePatternType(void){
    return CELibRunParameters.SolarAbundancePatternType;
}

int CELibGetRunParameterLifeTimeType(void){
    return CELibRunParameters.LifeTimeType;
}


int CELibGetRunParameterSNIIYieldsModificationP98(void){
    return CELibRunParameters.SNIIYieldsModificationP98;
}

double CELibGetRunParameterSNIIUpperMass(void){
    return CELibRunParameters.SNIIUpperMass;
}

double CELibSetRunParameterSNIILowerMass(void){
    return CELibRunParameters.SNIILowerMass;
}

int CELibGetRunParameterSNIIYieldsTableID(void){
    return CELibRunParameters.SNIIYieldsTableID;
}

double CELibGetRunParameterSNIIHyperNovaFraction(void){
    return CELibRunParameters.SNIIHyperNovaFraction;
}

int CELibSetRunParameterSNIIExplosionTime(void){
    return CELibRunParameters.SNIIExplosionTimeFollowYields;
}

int CELibGetRunParameterSNIaYieldsTableID(void){
    return CELibRunParameters.SNIaYieldsTableID;
}

int CELibGetRunParameterSNIaYieldsTableModelID(void){
    return CELibRunParameters.SNIaYieldsTableModelID;
}

int CELibGetRunParameterSNIaType(void){
    return CELibRunParameters.SNIaType;
}

double CELibGetRunParameterSNIaUpperMass(void){
    return CELibRunParameters.SNIaUpperMass;
}

double CELibGetRunParameterSNIaLowerMass(void){
    return CELibRunParameters.SNIaLowerMass;
}

double CELibGetRunParameterSNIaUpperMassBinary(void){
    return CELibRunParameters.SNIaUpperMassBinary;
}

double CELibGetRunParameterSNIaLowerMassBinary(void){
    return CELibRunParameters.SNIaLowerMassBinary;
}

int CELibGetRunParameterSNIaNassociation(void){
    return CELibRunParameters.SNIaNassociation;
}

int CELibGetRunParameterAGBYieldsTableID(void){ 
    return CELibRunParameters.AGBYieldsTableID;
}

double CELibGetRunParameterAGBUpperMass(void){
    return CELibRunParameters.AGBUpperMass;
}

double CELibGetRunParameterAGBLowerMass(void){
    return CELibRunParameters.AGBLowerMass;
}

int CELibGetRunParameterAGBBinNumber(void){
    return CELibRunParameters.AGBBinNumber;
}

int CELibGetRunParameterAGBBinType(void){
    return CELibRunParameters.AGBBinType;
}

double CELibGetRunParameterAGBBinTimeInterval(void){
    return CELibRunParameters.AGBBinTimeInterval;
}

int CELibGetRunParameterNSMYieldsTableID(void){
    return CELibRunParameters.NSMYieldsTableID;
}

int CELibGetRunParameterNSMType(void){
    return CELibRunParameters.NSMType;
}

double CELibGetRunParameterNSMDTDPowerLawIndex(void){
    return CELibRunParameters.NSMDTDPowerLawIndex;
}

double CELibGetRunParameterNSMDTDOffsetForPower(void){
    return CELibRunParameters.NSMDTDOffsetForPower;
}

int CELibGetRunParameterNSMNassociation(void){
    return CELibRunParameters.NSMNassociation;
}

int CELibGetRunParameterPopIIIIMF(void){
    return CELibRunParameters.PopIIIIMF;
}

int CELibGetRunParameterPopIIISNe(void){
    return CELibRunParameters.PopIIISNe;
}

int CELibGetRunParameterPopIIIAGB(void){
    return CELibRunParameters.PopIIIAGB;
}

int CELibGetRunParameterPopIIIAGBYieldsTableID(void){
    return CELibRunParameters.PopIIIAGBYieldsTableID;
}

double CELibGetRunParameterPopIIIMetallicity(void){
    return CELibRunParameters.PopIIIMetallicity;
}

int CELibGetRunParameterPopIIILifeTime(void){
    return CELibRunParameters.PopIIILifeTime;
}

