#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "Utilities.h"
#include "CELib.h"

static void pWriteAllIMFdata(void);

int main(int argc, char *argv[]){

    CELibShowVersion();
    CELibShowCurrentStatus();

    CELibSetRunParameterTestMode(true);
    CELibSetRunParameterIMFType(CELibIMF_Chabrier);
    CELibSetRunParameterSNIIYieldsTableID(CELibSNIIYieldsTableID_N13);
    CELibSetRunParameterSNIaType(1);
    CELibSetRunParameterSNIaYieldsTableID(CELibSNIaYieldsTableID_S13);
    CELibSetRunParameterSNIaYieldsTableModelID(CELibSNIaYieldsTableModelID_S13_N100ZDepend);
    CELibSetRunParameterAGBYieldsTableID(CELibAGBYieldsTableID_K10D14);
    CELibSetRunParameterPopIIIIMF(1);
    CELibSetRunParameterPopIIISNe(1);
    CELibSetRunParameterPopIIIAGB(1);
    CELibSetRunParameterPopIIIAGB(CELibAGBZ0YieldsTableID_CL08G13);

    CELibInit();

    pWriteAllIMFdata();
    CELibWriteIMFCumlative("./CELib/CELibIMF");

    CELibShowAllElementName();
    CELibShowAllElementID();

    CELibShowCurrentRunParameters();

    return EXIT_SUCCESS;
}

static void pWriteIMF(const int Type, const char OutputDir[]){

    CELibWriteIMFData(Type,OutputDir);

    return ;
}

static void pWriteAllIMFdata(void){

    MakeDir("./CELib");
    MakeDir("./CELib/CELibIMF");
    for(int i=0;i<CELibIMF_NTypes;i++){
        pWriteIMF(i,"./CELib/CELibIMF");
    }
    return ;
}

