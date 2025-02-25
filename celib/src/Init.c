#include "config.h"
#include "CheckRunParameters.h"
#include "IMFFunctions.h"
#include "SNIIYields.h"
#include "SNIIRate.h"
#include "SNIaYields.h"
#include "SNIaRate.h"
#include "NSMYields.h"
#include "NSMRate.h"
#include "InitialMetallicity.h"
#include "AGBMassLoss.h"
#include "YieldsInfo.h"

/*!
 * This function calls all of the initialization routines used in CELib.
 * Therefore, the user should call this function at the beginning of a
 * simulation. Otherwise, CELib returns incorrect values.
 */
void CELibInit(void){

    CheckRunParameters();

    if(CELibRunParameters.TestMode){
        MakeDir("./CELib");
    }
    CELibInitSolarAbundances();
    CELibInitYieldNamesIDs();

    CELibInitLifeTime();
    CELibInitIMF();

    CELibInitSNIIYields();
    CELibInitSNIIRate();

    CELibInitSNIaYields();
    CELibInitSNIaRate();

    CELibInitAGBYields();

    CELibInitNSMYields();
    CELibInitNSMRate();

    return ;
}

