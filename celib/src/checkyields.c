#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "Utilities.h"
#include "CELib.h"


static void CheckSNIIYieldsTableInfo(void){

    fprintf(stderr,"//\t SNII Yields table ID (RunParameters.SNIIYieldsTableID)\n");
    for(int i=0;i<CELibSNIIYieldsTableID_Number;i++){
        CELibSetRunParameterSNIIYieldsTableID(i);
        if(i==CELibSNIIYieldsTableID_P98){
            fprintf(stderr,"//\t\t CELibSNIIYieldsTableID_P98: %d\n",CELibRunParameters.SNIIYieldsTableID);
        } else if(i==CELibSNIIYieldsTableID_N13){ 
            fprintf(stderr,"//\t\t CELibSNIIYieldsTableID_N13: %d\n",CELibRunParameters.SNIIYieldsTableID);
        }
    }

    return ;
}

static void CheckSNIaYieldsTableInfo(void){

    fprintf(stderr,"//\t SNIa Yields table ID (RunParameters.SNIaYieldsTableID) \n");
    for(int i=0;i<CELibSNIaYieldsTableID_Number;i++){
        if(i == CELibSNIaYieldsTableID_I99){
            fprintf(stderr,"//\t\t CELibSNIaYieldsTableID_I99: %d\n",CELibRunParameters.SNIaYieldsTableID);
            fprintf(stderr,"//\t\t SNIa Yields table model ID (RunParameters.SNIaYieldsTableModelID) \n");
            for(int k=0;k<CELibSNIaYieldsTableModelID_I99_Number;k++){
                CELibSetRunParameterSNIaYieldsTableID(i);
                CELibSetRunParameterSNIaYieldsTableModelID(k);

                if(k==CELibSNIaYieldsTableModelID_I99_W7){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_W7: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_I99_W70){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_W70: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_I99_WDD1){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_WDD1: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_I99_WDD2){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_WDD2: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_I99_WDD3){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_WDD3: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_I99_CDD1){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_CDD1: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_I99_CDD2){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_I99_CDD2: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                }
            }
        } else if(i == CELibSNIaYieldsTableID_M10){
            fprintf(stderr,"//\t\t CELibSNIaYieldsTableID_M10: %d\n",CELibRunParameters.SNIaYieldsTableID);
            fprintf(stderr,"//\t\t SNIa Yields table model ID (RunParameters.SNIaYieldsTableModelID) \n");
            for(int k=0;k<CELibSNIaYieldsTableModelID_M10_Number;k++){
                CELibSetRunParameterSNIaYieldsTableID(i);
                CELibSetRunParameterSNIaYieldsTableModelID(k);

                if(k==CELibSNIaYieldsTableModelID_M10_W7){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_M10_W7: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_M10_CDEF){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_M10_CDEF: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_M10_CDDT){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_M10_CDDT: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_M10_ODDT){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_M10_ODDT: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                }
            }
        } else if(i == CELibSNIaYieldsTableID_S13){
            fprintf(stderr,"//\t\t CELibSNIaYieldsTableID_S13: %d\n",CELibRunParameters.SNIaYieldsTableID);
            fprintf(stderr,"//\t\t SNIa Yields table model ID (RunParameters.SNIaYieldsTableModelID) \n");
            for(int k=0;k<CELibSNIaYieldsTableModelID_S13_Number;k++){
                CELibSetRunParameterSNIaYieldsTableID(i);
                CELibSetRunParameterSNIaYieldsTableModelID(k);

                if(k==CELibSNIaYieldsTableModelID_S13_N1){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N1: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N3){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N3: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N5){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N5: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N10){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N10: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N20){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N20: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N40){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N40: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100H){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100H: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100L){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100L: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N150){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N150: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N200){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N200: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N300C){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N300C: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N1600){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N1600: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N1600C){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N1600C: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100Z05){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100Z05: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100Z01){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100Z01: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100Z001){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100Z001: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_S13_N100ZDepend){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_S13_N100ZDepend: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                }
            }
        } else if(i == CELibSNIaYieldsTableID_T04){
            fprintf(stderr,"//\t\t CELibSNIaYieldsTableID_T04: %d\n",CELibRunParameters.SNIaYieldsTableID);
            fprintf(stderr,"//\t\t SNIa Yields table model ID (RunParameters.SNIaYieldsTableModelID) \n");
            for(int k=0;k<CELibSNIaYieldsTableModelID_T04_Number;k++){
                CELibSetRunParameterSNIaYieldsTableID(i);
                CELibSetRunParameterSNIaYieldsTableModelID(k);

                if(k==CELibSNIaYieldsTableModelID_T04_W7){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_T04_W7: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_T04_c32d512){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_T04_c32d512: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_T04_c33d256a){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_T04_c33d256a: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_T04_c33d256b){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_T04_c33d256b: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_T04_b53d256){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_T04_b53d256: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                } else if(k==CELibSNIaYieldsTableModelID_T04_b303d768){
                    fprintf(stderr,"//\t\t\t CELibSNIaYieldsTableModelID_T04_b303d768: %d\n",CELibRunParameters.SNIaYieldsTableModelID);
                }
            }
        }
    }

    return ;
}

static void CheckAGBYieldsTableInfo(void){

    fprintf(stderr,"//\t AGB Yields table ID (RunParameters.AGBYieldsTableID)\n");
    for(int i=0;i<CELibAGBYieldsTableID_Number;i++){
        CELibSetRunParameterAGBYieldsTableID(i);
        if(i==CELibAGBYieldsTableID_K10){
            fprintf(stderr,"//\t\t CELibAGBYieldsTableID_K10: %d\n",CELibRunParameters.AGBYieldsTableID);
        } else if(i==CELibAGBYieldsTableID_K10D14){ 
            fprintf(stderr,"//\t\t CELibAGBYieldsTableID_K10D14: %d\n",CELibRunParameters.AGBYieldsTableID);
        } else if(i==CELibAGBYieldsTableID_vdHG97){ 
            fprintf(stderr,"//\t\t CELibAGBYieldsTableID_vdHG97: %d\n",CELibRunParameters.AGBYieldsTableID);
        }
    }

    return ;
}

static void CheckAGBZ0YieldsTableInfo(void){

    fprintf(stderr,"//\t AGB Yields table ID (RunParameters.PopIIIAGBYieldsTableID)\n");
    for(int i=0;i<CELibAGBZ0YieldsTableID_Number;i++){
        CELibSetRunParameterPopIIIAGBYieldsTableID(i);
        if(i==CELibAGBZ0YieldsTableID_CL08){
            fprintf(stderr,"//\t\t CELibAGBZ0YieldsTableID_CL08: %d\n",CELibRunParameters.PopIIIAGBYieldsTableID);
        } else if(i==CELibAGBYieldsTableID_K10D14){ 
            fprintf(stderr,"//\t\t CELibAGBZ0YieldsTableID_CL08G13: %d\n",CELibRunParameters.PopIIIAGBYieldsTableID);
        }
    }

    return ;
}

static void CheckNSMYieldsTableInfo(void){

    fprintf(stderr,"//\t NSM Yields table ID (RunParameters.NSMYieldsTableID)\n");
    for(int i=0;i<CELibNSMYieldsTableID_Number;i++){
        CELibSetRunParameterNSMYieldsTableID(i);
        if(i==CELibNSMYieldsTableID_W14){
            fprintf(stderr,"//\t\t CELibNSMYieldsTableID_W14: %d\n",CELibRunParameters.NSMYieldsTableID);
        }
    }

    return ;
}

int main(int argc, char *argv[]){

    CELibShowVersion();
    CELibShowCurrentStatus();

    CheckSNIIYieldsTableInfo();
    CheckSNIaYieldsTableInfo();
    CheckAGBYieldsTableInfo();
    CheckAGBZ0YieldsTableInfo();
    CheckNSMYieldsTableInfo();

    return EXIT_SUCCESS;
}
