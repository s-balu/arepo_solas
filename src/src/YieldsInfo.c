#include "config.h"
#include "RunParameters.h"
#include "IMFFunctions.h"
#include "SNIIYields.h"
#include "SNIaYields.h"

/*! \file YieldsInfo.c
 * \brief This file has functions handling yield information.
 */

static char CELibYieldNameTags[CELibYield_Number][MaxCharactersInLine];
static char CELibSNIIYieldNameTags[CELibYield_Number][MaxCharactersInLine];
static char CELibSNIaYieldNameTags[CELibYield_Number][MaxCharactersInLine];
static char CELibAGBYieldNameTags[CELibYield_Number][MaxCharactersInLine];
static char CELibNSMYieldNameTags[CELibYield_Number][MaxCharactersInLine];

/*!
 * Set yield names.
 */
void CELibInitYieldNamesIDs(void){

    //Yields
    strcpy(CELibYieldNameTags[CELibYield_H],    "H");
    strcpy(CELibYieldNameTags[CELibYield_He],   "He");
    strcpy(CELibYieldNameTags[CELibYield_C],    "C");
    strcpy(CELibYieldNameTags[CELibYield_N],    "N");
    strcpy(CELibYieldNameTags[CELibYield_O],    "O");
    strcpy(CELibYieldNameTags[CELibYield_Ne],   "Ne");
    strcpy(CELibYieldNameTags[CELibYield_Mg],   "Mg");
    strcpy(CELibYieldNameTags[CELibYield_Si],   "Si");
    strcpy(CELibYieldNameTags[CELibYield_S],    "S");
    strcpy(CELibYieldNameTags[CELibYield_Ca],   "Ca");
    strcpy(CELibYieldNameTags[CELibYield_Fe],   "Fe");
    strcpy(CELibYieldNameTags[CELibYield_Ni],   "Ni");
    strcpy(CELibYieldNameTags[CELibYield_Eu],   "Eu");

    //SNII Yields
    strcpy(CELibSNIIYieldNameTags[CELibYield_H],    "H");
    strcpy(CELibSNIIYieldNameTags[CELibYield_He],   "He");
    strcpy(CELibSNIIYieldNameTags[CELibYield_C],    "C");
    strcpy(CELibSNIIYieldNameTags[CELibYield_N],    "N");
    strcpy(CELibSNIIYieldNameTags[CELibYield_O],    "O");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Ne],   "Ne");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Mg],   "Mg");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Si],   "Si");
    strcpy(CELibSNIIYieldNameTags[CELibYield_S],    "S");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Ca],   "Ca");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Fe],   "Fe");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Ni],   "Ni");
    strcpy(CELibSNIIYieldNameTags[CELibYield_Eu],   "Eu");

    //SNIa Yields
    strcpy(CELibSNIaYieldNameTags[CELibYield_H],    "H");
    strcpy(CELibSNIaYieldNameTags[CELibYield_He],   "He");
    strcpy(CELibSNIaYieldNameTags[CELibYield_C],    "C");
    strcpy(CELibSNIaYieldNameTags[CELibYield_N],    "N");
    strcpy(CELibSNIaYieldNameTags[CELibYield_O],    "O");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Ne],   "Ne");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Mg],   "Mg");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Si],   "Si");
    strcpy(CELibSNIaYieldNameTags[CELibYield_S],    "S");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Ca],   "Ca");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Fe],   "Fe");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Ni],   "Ni");
    strcpy(CELibSNIaYieldNameTags[CELibYield_Eu],   "Eu");

    //AGB Yields
    strcpy(CELibAGBYieldNameTags[CELibYield_H],    "H");
    strcpy(CELibAGBYieldNameTags[CELibYield_He],   "He");
    strcpy(CELibAGBYieldNameTags[CELibYield_C],    "C");
    strcpy(CELibAGBYieldNameTags[CELibYield_N],    "N");
    strcpy(CELibAGBYieldNameTags[CELibYield_O],    "O");
    strcpy(CELibAGBYieldNameTags[CELibYield_Ne],   "Ne");
    strcpy(CELibAGBYieldNameTags[CELibYield_Mg],   "Mg");
    strcpy(CELibAGBYieldNameTags[CELibYield_Si],   "Si");
    strcpy(CELibAGBYieldNameTags[CELibYield_S],    "S");
    strcpy(CELibAGBYieldNameTags[CELibYield_Ca],   "Ca");
    strcpy(CELibAGBYieldNameTags[CELibYield_Fe],   "Fe");
    strcpy(CELibAGBYieldNameTags[CELibYield_Ni],   "Ni");
    strcpy(CELibAGBYieldNameTags[CELibYield_Eu],   "Eu");

    //NSM
    strcpy(CELibNSMYieldNameTags[CELibYield_H],    "H");
    strcpy(CELibNSMYieldNameTags[CELibYield_He],   "He");
    strcpy(CELibNSMYieldNameTags[CELibYield_C],    "C");
    strcpy(CELibNSMYieldNameTags[CELibYield_N],    "N");
    strcpy(CELibNSMYieldNameTags[CELibYield_O],    "O");
    strcpy(CELibNSMYieldNameTags[CELibYield_Ne],   "Ne");
    strcpy(CELibNSMYieldNameTags[CELibYield_Mg],   "Mg");
    strcpy(CELibNSMYieldNameTags[CELibYield_Si],   "Si");
    strcpy(CELibNSMYieldNameTags[CELibYield_S],    "S");
    strcpy(CELibNSMYieldNameTags[CELibYield_Ca],   "Ca");
    strcpy(CELibNSMYieldNameTags[CELibYield_Fe],   "Fe");
    strcpy(CELibNSMYieldNameTags[CELibYield_Ni],   "Ni");
    strcpy(CELibNSMYieldNameTags[CELibYield_Eu],   "Eu");

    return ;
}

/*!
 * This function retunrs the "ElementID" of the yield "name[]".
 */
int CELibGetSNIIYieldElementID(const char *name){
    for(int i=0;i<CELibYield_Number;i++){
        int ret = strcmp(name,CELibSNIIYieldNameTags[i]);
        if(ret == 0){
            return i;
        }
    }
    return NONE;
}

/*!
 * This function retunrs the "ElementID" of the yield "name[]".
 */
int CELibGetSNIaYieldElementID(const char *name){
    for(int i=0;i<CELibYield_Number;i++){
        int ret = strcmp(name,CELibSNIaYieldNameTags[i]);
        if(ret == 0){
            return i;
        }
    }
    return NONE;
}

/*!
 * This function retunrs the "ElementID" of the yield "name[]".
 */
int CELibGetAGBYieldElementID(const char *name){
    for(int i=0;i<CELibYield_Number;i++){
        int ret = strcmp(name,CELibAGBYieldNameTags[i]);
        if(ret == 0){
            return i;
        }
    }
    return NONE;
}

/*!
 * This function retunrs the "ElementID" of the yield "name[]".
 */
int CELibGetNSMYieldElementID(const char *name){
    for(int i=0;i<CELibYield_Number;i++){
        int ret = strcmp(name,CELibNSMYieldNameTags[i]);
        if(ret == 0){
            return i;
        }
    }
    return NONE;
}

/*!
 * This function retunrs the "ElementID" of the yield "name[]".
 */
int CELibGetYieldElementID(const char *name){

    for(int i=0;i<CELibYield_Number;i++){
        int ret = strcmp(name,CELibYieldNameTags[i]);
        if(ret == 0){
            return i;
        }
    }
    return NONE;

    return 0;
}

/*!
 * This function retunrs the "ElementName" of the yield ID, "ID".
 */
void CELibGetSNIIYieldElementName(const int ID, char *name){

    if(ID<CELibYield_Number){
        strncpy(name,CELibSNIIYieldNameTags[ID],MaxCharactersInLine);
    } else {
        strncpy(name,"",MaxCharactersInLine);
    }

    return ;
}

/*!
 * This function retunrs the "ElementName" of the yield ID, "ID".
 */
void CELibGetSNIaYieldElementName(const int ID, char *name){

    if(ID<CELibYield_Number){
        strncpy(name,CELibSNIaYieldNameTags[ID],MaxCharactersInLine);
    } else {
        strncpy(name,"",MaxCharactersInLine);
    }

    return ;
}

/*!
 * This function retunrs the "ElementName" of the yield ID, "ID".
 */
void CELibGetAGBYieldElementName(const int ID, char *name){

    if(ID<CELibYield_Number){
        strncpy(name,CELibAGBYieldNameTags[ID],MaxCharactersInLine);
    } else {
        strncpy(name,"",MaxCharactersInLine);
    }

    return ;
}

/*!
 * This function retunrs the "ElementName" of the yield ID, "ID".
 */
void CELibGetNSMYieldElementName(const int ID, char *name){

    if(ID<CELibYield_Number){
        strncpy(name,CELibNSMYieldNameTags[ID],MaxCharactersInLine);
    } else {
        strncpy(name,"",MaxCharactersInLine);
    }

    return ;
}

/*!
 * This function retunrs the "ElementName" of the yield ID, "ID".
 */
void CELibGetYieldElementName(const int ID, char *name){

    if(ID<CELibYield_Number){
        strncpy(name,CELibYieldNameTags[ID],MaxCharactersInLine);
    } else {
        strncpy(name,"",MaxCharactersInLine);
    }

    return ;
}

/*!
 * This function shows all elements names.
 */
void CELibShowAllElementName(void){

    fprintf(stderr,"//\t Show element name:\n");
    for(int i=0;i<CELibYield_Number;i++){
        char Name[MaxCharactersInLine];
        CELibGetYieldElementName(i,Name);
        fprintf(stderr,"//\t\t Element[%d] = %s\n",i,Name);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show SNII element name:\n");
    for(int i=0;i<CELibYield_Number;i++){
        char Name[MaxCharactersInLine];
        CELibGetSNIIYieldElementName(i,Name);
        fprintf(stderr,"//\t\t SNII element[%d] = %s\n",i,Name);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show SNIa element name:\n");
    for(int i=0;i<CELibYield_Number;i++){
        char Name[MaxCharactersInLine];
        CELibGetSNIaYieldElementName(i,Name);
        fprintf(stderr,"//\t\t SNIa element[%d] = %s\n",i,Name);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show AGB element name:\n");
    for(int i=0;i<CELibYield_Number;i++){
        char Name[MaxCharactersInLine];
        CELibGetAGBYieldElementName(i,Name);
        fprintf(stderr,"//\t\t AGB element[%d] = %s\n",i,Name);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show NSM element name:\n");
    for(int i=0;i<CELibYield_Number;i++){
        char Name[MaxCharactersInLine];
        CELibGetNSMYieldElementName(i,Name);
        fprintf(stderr,"//\t\t NSM element[%d] = %s\n",i,Name);
    }
    fprintf(stderr,"\n");
    return ;
}

/*!
 * This function shows all elements IDs.
 */
void CELibShowAllElementID(void){

    fprintf(stderr,"//\t Show element ID:\n");
    for(int i=0;i<CELibYield_Number;i++){
        int ID = CELibGetYieldElementID(CELibYieldNameTags[i]);
        fprintf(stderr,"//\t\t The ID of %s = %d\n",CELibYieldNameTags[i],ID);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show SNII element ID:\n");
    for(int i=0;i<CELibYield_Number;i++){
        int ID = CELibGetSNIIYieldElementID(CELibSNIIYieldNameTags[i]);
        fprintf(stderr,"//\t\t The ID of SNII %s = %d\n",CELibSNIIYieldNameTags[i],ID);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show SNIa element ID:\n");
    for(int i=0;i<CELibYield_Number;i++){
        int ID = CELibGetSNIaYieldElementID(CELibSNIaYieldNameTags[i]);
        fprintf(stderr,"//\t\t The ID of SNIa %s = %d\n",CELibSNIaYieldNameTags[i],ID);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show AGB element ID:\n");
    for(int i=0;i<CELibYield_Number;i++){
        int ID = CELibGetAGBYieldElementID(CELibAGBYieldNameTags[i]);
        fprintf(stderr,"//\t\t The ID of AGB %s = %d\n",CELibAGBYieldNameTags[i],ID);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"//\t Show NSM element ID:\n");
    for(int i=0;i<CELibYield_Number;i++){
        int ID = CELibGetNSMYieldElementID(CELibNSMYieldNameTags[i]);
        fprintf(stderr,"//\t\t The ID of NSM %s = %d\n",CELibNSMYieldNameTags[i],ID);
    }
    fprintf(stderr,"\n");

    return ;
}



