#include "config.h"
#include "SNIaYields.h"
#include "../data/Iwamoto+1999/Iwamoto+1999_Struct.h"
#include "../data/Maeda+2010/Maeda+2010_Struct.h"
#include "../data/Seitenzahl+2013/Seitenzahl+2013_Struct.h"
#include "../data/Travaglio+2004/Travaglio+2004_Struct.h"


/*! \file 
 * \brief CELib run parameters check function is defined in this file. 
 */

static int ErrorCount = 0;

#define CELibCheckRangeInteger(f,a,b)                      \
    if(f < a){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %d, which should be >= %d\n",f,a); \
        ErrorCount ++;                                     \
    }                                                      \
    if(f > b){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %d, which should be <= %d\n",f,b); \
        ErrorCount ++;                                     \
    }                                                      \

#define CELibCheckRangeIntegerMin(f,a)                     \
    if(f < a){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %d, which should be >= %d\n",f,a); \
        ErrorCount ++;                                     \
    }                                                      \

#define CELibCheckRangeIntegerMax(f,a)                      \
    if(f > b){                                              \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %d, which should be <= %d\n",f,b); \
        ErrorCount ++;                                     \
    }                                                      \

#define CELibCheckRangeDouble(f,a,b)                       \
    if(f < a){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %g, which should be >= %g\n",f,a); \
        ErrorCount ++;                                     \
    }                                                      \
    if(f > b){                                             \
        ErrorCount ++;                                     \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %g, which should be <= %g\n",f,b); \
    }
\
#define CELibCheckRangeDoubleMin(f,a)                      \
    if(f < a){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %g, which should be >= %g\n",f,a); \
        ErrorCount ++;                                     \
    }                                                      \

#define CELibCheckRangeDoubleMax(f,a)                      \
    if(f > b){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" is out-of-range, The value is %g, which should be <= %g\n",f,b); \
        ErrorCount ++;                                     \
    }                                                      \

#define CELibCheckLargeSmallRelation(f,g)                  \
    if(f > g){                                             \
        fprintf(stderr,"//\tCELib RunParameter check: "#f" (%g) is larger than "#g" (%g), They should be "#f" < "#g"\n",f,g); \
        ErrorCount ++;                                     \
    }                                                      \



/*
 * All RunParameters are checked in this function. This function should
 * be called at the beginning of the CELib initialization phase.
 */
void CheckRunParameters(void){

    // General
    CELibCheckRangeInteger(CELibRunParameters.IMFType,0,CELibIMF_NTypes-1);
    CELibCheckRangeInteger(CELibRunParameters.SolarAbundancePatternType,0,2);


    // LifeTime
    CELibCheckRangeInteger(CELibRunParameters.LifeTimeType,0,1);


    // SNe II
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIIEnergy,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIIEnergyEfficiency,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIIUpperMass,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIILowerMass,0.0);
    CELibCheckLargeSmallRelation(CELibRunParameters.SNIILowerMass,CELibRunParameters.SNIIUpperMass);

    CELibCheckRangeInteger(CELibRunParameters.SNIIYieldsTableID,0,1);
    CELibCheckRangeInteger(CELibRunParameters.SNIIYieldsModificationP98,0,1);
    CELibCheckRangeDouble(CELibRunParameters.SNIIHyperNovaFraction,0.0,1.0);
    CELibCheckRangeInteger(CELibRunParameters.SNIIExplosionTimeFollowYields,0,1);


    // SNe Ia
    CELibCheckRangeInteger(CELibRunParameters.SNIaType,0,2);
    CELibCheckRangeInteger(CELibRunParameters.SNIaYieldsTableID,0,3);
    if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableModelID_I99_Number){
        CELibCheckRangeInteger(CELibRunParameters.SNIaYieldsTableModelID,0,CELibSNIaYieldsTableModelID_I99_Number-1);
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableModelID_M10_Number){
        CELibCheckRangeInteger(CELibRunParameters.SNIaYieldsTableModelID,0,CELibSNIaYieldsTableModelID_M10_Number-1);
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableModelID_S13_Number){
        // Note that count the metal depending model.
        CELibCheckRangeInteger(CELibRunParameters.SNIaYieldsTableModelID,0,CELibSNIaYieldsTableModelID_S13_Number); 
    } else if(CELibRunParameters.SNIaYieldsTableID == CELibSNIaYieldsTableModelID_T04_Number){
        CELibCheckRangeInteger(CELibRunParameters.SNIaYieldsTableModelID,0,CELibSNIaYieldsTableModelID_T04_Number-1);
    }
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIaEnergy,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIaEnergyEfficiency,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIaUpperMassBinary,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.SNIaLowerMassBinary,0.0);
    CELibCheckLargeSmallRelation(CELibRunParameters.SNIaLowerMassBinary,CELibRunParameters.SNIaUpperMassBinary);


    // AGB
    CELibCheckRangeInteger(CELibRunParameters.AGBYieldsTableID,0,2);
    CELibCheckRangeIntegerMin(CELibRunParameters.AGBBinNumber,0);
    CELibCheckRangeInteger(CELibRunParameters.AGBBinType,0,1);
    CELibCheckRangeDouble(CELibRunParameters.AGBBinTimeInterval,0.0,1.38e10);
    CELibCheckRangeDouble(CELibRunParameters.AGBBinUpperAge,0.0,1.38e10);
    CELibCheckRangeDouble(CELibRunParameters.AGBBinLowerAge,0.0,1.38e10);
    CELibCheckLargeSmallRelation(CELibRunParameters.AGBBinLowerAge,CELibRunParameters.AGBBinUpperAge);
    CELibCheckRangeDoubleMin(CELibRunParameters.AGBUpperMass,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.AGBLowerMass,0.0);
    CELibCheckLargeSmallRelation(CELibRunParameters.AGBLowerMass,CELibRunParameters.AGBUpperMass);


    // NSM
    CELibCheckRangeInteger(CELibRunParameters.NSMType,0,1);
    CELibCheckRangeDoubleMin(CELibRunParameters.NSMNumberPerMass,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.NSMUpperMass,0.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.NSMLowerMass,0.0);
    CELibCheckLargeSmallRelation(CELibRunParameters.NSMLowerMass,CELibRunParameters.NSMUpperMass);
    CELibCheckRangeDouble(CELibRunParameters.NSMFraction,0.0,1.0);
    CELibCheckRangeDoubleMin(CELibRunParameters.NSMDTDOffsetForPower,0.0);
    CELibCheckRangeIntegerMin(CELibRunParameters.NSMNassociation,1);
    CELibCheckRangeIntegerMin(CELibRunParameters.NSMYieldsTableID,0);


    // PopIII
    CELibCheckRangeInteger(CELibRunParameters.PopIIIIMF,0,1);
    CELibCheckRangeInteger(CELibRunParameters.PopIIISNe,0,1);
    CELibCheckRangeInteger(CELibRunParameters.PopIIIAGB,0,1);
    CELibCheckRangeInteger(CELibRunParameters.PopIIIAGBYieldsTableID,0,1);
    CELibCheckRangeDouble(CELibRunParameters.PopIIIMetallicity,0.0,1.0);
    CELibCheckRangeInteger(CELibRunParameters.PopIIILifeTime,0,1);

    if(ErrorCount>0){
        fprintf(stderr,"%d errors are found in CELibRunParamters.\n",ErrorCount);
    }

    return ;
}
