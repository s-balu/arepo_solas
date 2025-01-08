#include "config.h"
#include "SNIaYields.h"
#include "../data/Iwamoto+1999/Iwamoto+1999_Struct.h"
#include "../data/Maeda+2010/Maeda+2010_Struct.h"
#include "../data/Seitenzahl+2013/Seitenzahl+2013_Struct.h"
#include "../data/Travaglio+2004/Travaglio+2004_Struct.h"


/*! \file Info.h
 * \brief This file has functions which show version information and control
 * flags.
 */

/*!
 *  This function shows information of this library.
 */
void CELibShowVersion(void){
    
    fprintf(stderr,"// CELib, chemical evolution library for galaxy formation simulations.\n");
    fprintf(stderr,"//\t This library is developed by Takayuki Saitoh.\n");
    fprintf(stderr,"//\t Version %d.%d.%d\n",
            __CELib_MajorVersion__,__CELib_MinorVersion__,__CELib_MicroVersion__);
    fprintf(stderr,"//\t Git repository hash: %s\n",GITVERSION);
    fflush(stderr);

    return ;
}


/*!
 *  This function shows current control flags.
 */
void CELibShowCurrentStatus(void){

    /// Mode
    fprintf(stderr,"\n");
    fprintf(stderr,"//__________For_Mode__\n");
    if(CELibRunParameters.TestMode == true){
        fprintf(stderr,"\tTest mode: on\n");
    }else {
        fprintf(stderr,"\tTest mode: on\n");
    }
 
    /// IMF
    fprintf(stderr,"\n");
    fprintf(stderr,"//__________For_IMF__\n");
    fprintf(stderr,"\tIMF: %s\n",CELibIMF[CELibRunParameters.IMFType].Name);
    fflush(stderr);

    /// Lifetime table
    fprintf(stderr,"\n");
    fprintf(stderr,"//__________For_Lifetime__\n");
    if(CELibRunParameters.LifeTimeType == CELibLifeTime_Original){
        fprintf(stderr,"\tOriginal stellar lifetime table data is used.\n");
    } else if(CELibRunParameters.LifeTimeType == CELibLifeTime_LSF){
        fprintf(stderr,"\tSmoothed stellar lifetime is used.\n");
    }
    fprintf(stderr,"\n");

    /// SNe II
    fprintf(stderr,"//__________For_SNe_II__\n");
    if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
        fprintf(stderr,"\tSNe II yields table: Portinari+1998\n");
    } else if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_N13){
        fprintf(stderr,"\tSNe II yields table: Nomot+2013\n");
    }
    if(CELibRunParameters.SNIIYieldsTableID==CELibSNIIYieldsTableID_P98){
        if(CELibRunParameters.SNIIYieldsModificationP98){
            fprintf(stderr,"\tModifications of yields are used\n");
        }
        fprintf(stderr,"\tSNe II yield table ID is %d\n",CELibRunParameters.SNIIYieldsTableID);
    }
    if(CELibRunParameters.SNIIYieldsTableID==CELibSNIIYieldsTableID_N13){
        fprintf(stderr,"\tHNe fraction is %g\n",CELibRunParameters.SNIIHyperNovaFraction);
    }
    if(CELibRunParameters.SNIIExplosionTimeFollowYields==0){
        fprintf(stderr,"\tSNII explosion time range is defined using stellar ages of SNIILowerMass and SNIIUpperMass\n");
    } else {
        fprintf(stderr,"\tSNII explosion time range follows the defined range of the adopted yields table\n");
    }
    fprintf(stderr,"\n");

    /// SNe Ia
    fprintf(stderr,"//__________For_SNe_Ia__\n");
    if(CELibRunParameters.SNIaType == CELibSNIaRateModelID_GR83){
        fprintf(stderr,"\tSNe Ia model is Gregio & Renzine (1983)\n");
        fprintf(stderr,"\tMinimum total mass of a binary system is %g [Msun]\n",CELibRunParameters.SNIaLowerMassBinary);
        fprintf(stderr,"\tMaximum total mass of a binary system is %g [Msun]\n",CELibRunParameters.SNIaUpperMassBinary);
    } else {
        fprintf(stderr,"\tSNe Ia model is a Power law type DTD\n");
    }
    if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_I99){
        fprintf(stderr,"\tSNe Ia yields table based on Iwamoto+1999\n");
        fprintf(stderr,"\t\tModel name is %s\n",CELibSNIaYieldsI99[CELibRunParameters.SNIaYieldsTableModelID].Name);
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_M10){
        fprintf(stderr,"\tSNe Ia yields table based on Maeda+2010\n");
        fprintf(stderr,"\t\tModel name is %s\n",CELibSNIaYieldsM10[CELibRunParameters.SNIaYieldsTableModelID].Name);
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_S13){
        fprintf(stderr,"\tSNe Ia yields table based on Seitenzahl+2013\n");
        if(CELibRunParameters.SNIaYieldsTableModelID == CELibSNIaYieldsTableModelID_S13_N100ZDepend){
            fprintf(stderr,"\t\tModel name is Seitenzahl+2013, metal depend version\n");
        } else {
            fprintf(stderr,"\t\tModel name is %s\n",CELibSNIaYieldsS13[CELibRunParameters.SNIaYieldsTableModelID].Name);
        }
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableID_T04){
        fprintf(stderr,"\tSNe Ia yields table based on Travaglio+2004\n");
        fprintf(stderr,"\t\tModel name is %s\n",CELibSNIaYieldsT04[CELibRunParameters.SNIaYieldsTableModelID].Name);
    }
    fprintf(stderr,"\tSNe Ia yields table ID is %d\n",CELibRunParameters.SNIIYieldsTableID);
    fprintf(stderr,"\tSNe Ia type is %d\n",CELibRunParameters.SNIIYieldsTableID);
    fprintf(stderr,"\tSNe Ia number of association for a cluster mode: %d\n",CELibRunParameters.SNIaNassociation);

    fprintf(stderr,"\tReleased energy per SN Ia is %g [erg]\n",CELibRunParameters.SNIaEnergy);
    fprintf(stderr,"\n");

    /// AGB
    fprintf(stderr,"//__________For_AGB__\n");
    if(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_K10){
        fprintf(stderr,"\tAGB yeilds talbe: Karakas 2010\n");
    } else if(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_K10D14){
        fprintf(stderr,"\tAGB yeilds talbe: Karakas 2010 & Doherty+2013\n");
    }
    if(CELibRunParameters.AGBBinType == CELibAGBRateModelID_LogBin){
        fprintf(stderr,"\tLog bin is used for AGB feedback\n");
    } else {
        fprintf(stderr,"\tLinear bin is used for AGB feedback\n");
    }
    fprintf(stderr,"\tTime interval for AGB feedback is %d\n",CELibRunParameters.AGBYieldsTableID);
    fprintf(stderr,"\tLower age for AGB feedback is %g\n",CELibRunParameters.AGBBinLowerAge);
    fprintf(stderr,"\tUpper age for AGB feedback is %g\n",CELibRunParameters.AGBBinUpperAge);
    fprintf(stderr,"\tLower mass for AGBs: %g [Msun]\n",CELibRunParameters.AGBLowerMass);
    fprintf(stderr,"\tUpper mass for AGBS: %g [Msun]\n",CELibRunParameters.AGBUpperMass);
    fprintf(stderr,"\n");

    /// NSM
    fprintf(stderr,"//__________For_NSM__\n");
    if(CELibRunParameters.NSMType == CELibNSMRateModelID_PowerLaw){
        fprintf(stderr,"\tPower law type DTD is used.\n");
        fprintf(stderr,"\tDTD power law index: %g\n",CELibRunParameters.NSMDTDPowerLawIndex);
        fprintf(stderr,"\tOffset for DTD: %g\n",CELibRunParameters.NSMDTDOffsetForPower);
    }
    fprintf(stderr,"\tLower mass for NSMs: %g [Msun]\n",CELibRunParameters.NSMLowerMass);
    fprintf(stderr,"\tUpper mass for NSMs: %g [Msun]\n",CELibRunParameters.NSMUpperMass);
    fprintf(stderr,"\tNSMs fraction: %g\n",CELibRunParameters.NSMFraction);
    fprintf(stderr,"\tNSMs number of association for a cluster mode: %d\n",CELibRunParameters.NSMNassociation);
    fprintf(stderr,"\tReleased energy per NSM: %g [erg]\n",CELibRunParameters.NSMEnergy);
    if(CELibRunParameters.NSMYieldsTableID == CELibNSMYieldsTableID_W14){
        fprintf(stderr,"\tNSMs yield table: Wanajo+2014\n");
    }
    fprintf(stderr,"\n");

    /// PopIII
    fprintf(stderr,"//__________For_PopIII__\n");
    if(CELibRunParameters.PopIIIIMF == 1){
        fprintf(stderr,"\tPop III mode: on\n");
    }else {
        fprintf(stderr,"\tPop III mode: off\n");
    }
    if(CELibRunParameters.PopIIISNe == 1){
        fprintf(stderr,"\tPop III SNe mode: on\n");
    }else {
        fprintf(stderr,"\tPop III SNe mode: off\n");
    }
    if(CELibRunParameters.PopIIIAGB > 0){
        fprintf(stderr,"\tPop III AGB mode: on\n");
        if(CELibRunParameters.PopIIIAGBYieldsTableID == CELibAGBZ0YieldsTableID_CL08){
            fprintf(stderr,"\tAGB yields table of Campbell & Lattanzio 2008 is used.\n");
        } else if(CELibRunParameters.PopIIIAGBYieldsTableID == CELibAGBZ0YieldsTableID_CL08G13){
            fprintf(stderr,"\tAGB yields tables of Campbell & Lattanzio 2008 and Gil-Pons+2013 are used.\n");
        }
    }else {
        fprintf(stderr,"\tPop III AGB mode: off\n");
    }
    if(CELibRunParameters.PopIIIIMF == 1){
        fprintf(stderr,"\tPop III Z = %g Zsun\n",CELibRunParameters.PopIIIMetallicity/
                    CELibGetMetalFractionForSolarChemicalComposision());
    }
    fprintf(stderr,"\n");

    return ;
}
