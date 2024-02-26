#include "config.h"
#include "../data/Portinari+1998/Portinari+1998_Struct.h"

double MetallicitySNII[5] = {0.0004,0.004,0.008,0.02,0.05};

/*! \file SNIIYields.c
 * \brief This file has functions regarding SNII yields.
 */

/*!
 * Following variables hold the grid size of the adopted yields table and values
 * of metallicity.
 */
int CELibSNIIYields_Metallicity; 
int CELibSNIIYields_Mass; 

struct CELibStructSNIIYields *pCELibSNIIYields; 
#define CELibSNIIYields(x,y) (pCELibSNIIYields[CELibSNIIYields_Mass*(x)+(y)])


static char SNIIYieldsElementName[CELibYield_Number][MaxCharactersInLine]; 

/*!
 * Set yield names.
 */
static void pCELibSNIIYieldsSetElementName(void){

    strcpy(SNIIYieldsElementName[CELibYield_H],    "H"); 
    strcpy(SNIIYieldsElementName[CELibYield_He],   "He"); 
    strcpy(SNIIYieldsElementName[CELibYield_C],    "C"); 
    strcpy(SNIIYieldsElementName[CELibYield_N],    "N"); 
    strcpy(SNIIYieldsElementName[CELibYield_O],    "O"); 
    strcpy(SNIIYieldsElementName[CELibYield_Ne],   "Ne"); 
    strcpy(SNIIYieldsElementName[CELibYield_Mg],   "Mg"); 
    strcpy(SNIIYieldsElementName[CELibYield_Si],   "Si"); 
    strcpy(SNIIYieldsElementName[CELibYield_S],    "S"); 
    strcpy(SNIIYieldsElementName[CELibYield_Ca],   "Ca"); 
    strcpy(SNIIYieldsElementName[CELibYield_Fe],   "Fe"); 
    strcpy(SNIIYieldsElementName[CELibYield_Ni],   "Ni"); 
    strcpy(SNIIYieldsElementName[CELibYield_Eu],   "Eu"); 

    return ;
}

/*!
 * This function initializes the yield table of type II SNe. 
 */
void CELibInitSNIIYields(void){

    if(CELibRunParameters.SNIIYieldsTableID == CELibSNIIYieldsTableID_P98){
        CELibSNIIYields_Metallicity = CELibSNIIYields_Metallicity_P98; 
        CELibSNIIYields_Mass = CELibSNIIYields_Mass_P98;  
        pCELibSNIIYields = CELibSNIIYieldsP98[0]; 

    }

    pCELibSNIIYieldsSetElementName();

    return ;
}

/*! 
 * This function retuns an interpolated yeild of the element "TargetElement"
 * from a star of which mass is "Mass".
 */
static double pCELibGetSNIIYieldsInterpolatedYieldMass(const int IndexMetal, const int TargetElement, const double Mass){

    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    int IndexMass = 0;
    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Elements[TargetElement]-CELibSNIIYields(IndexMetal,IndexMass).Elements[TargetElement])/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Elements[TargetElement];
}

/*!
 *  This function retuns an interpolated remnant mass from a given mass star.
 */
static double pCELibGetSNIIYieldsInterpolatedRemnantMass(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Mr-CELibSNIIYields(IndexMetal,IndexMass).Mr)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Mr;
}

/*!
 * This function retuns an interpolated mass released via stellar winds from a
 * given mass star.
 */
static double pCELibGetSNIIYieldsInterpolatedStellarWindMass(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Msw-CELibSNIIYields(IndexMetal,IndexMass).Msw)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Msw;
}

/*! 
 * This function retuns an interpolated mass released via SN from a given mass
 * star.
 */
static double pCELibGetSNIIYieldsInterpolatedSuperNovaMass(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Msn-CELibSNIIYields(IndexMetal,IndexMass).Msn)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Msn;
}

/*! 
 * This function retuns the interpolated energy released from a given mass star.
 */
static double pCELibGetSNIIYieldsInterpolatedEnergy(const int IndexMetal, const double Mass){

    int IndexMass = 0;
    if(CELibSNIIYields(IndexMetal,0).Mass > Mass)
        return 0.e0;
    if(CELibSNIIYields(IndexMetal,CELibSNIIYields_Mass-1).Mass < Mass)
        return 0.e0;

    for(int i=1;i<CELibSNIIYields_Mass;i++){
        if(CELibSNIIYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }
    
    return (CELibSNIIYields(IndexMetal,IndexMass+1).Energy-CELibSNIIYields(IndexMetal,IndexMass).Energy)/
        (CELibSNIIYields(IndexMetal,IndexMass+1).Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibSNIIYields(IndexMetal,IndexMass).Mass)+CELibSNIIYields(IndexMetal,IndexMass).Energy;
}

/*!
 * This function returns energy, masses of individual elements, and the sum of
 * them (ejecta mass) released from a single SNII.  It also returns the remnant
 * mass.
 */
struct CELibStructFeedbackStarbyStarOutput CELibGetSNIIFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input){

    struct CELibStructFeedbackStarbyStarOutput SNIIStarbyStarFeedback;

    double Mstar = Input.Mass*Input.MassConversionFactor;

    if(Input.Metallicity <= MetallicitySNII[0]){
        int TableID = 0;

        double Msw = pCELibGetSNIIYieldsInterpolatedStellarWindMass(TableID,Mstar);
        double Msn = pCELibGetSNIIYieldsInterpolatedSuperNovaMass(TableID,Mstar);

        SNIIStarbyStarFeedback.EjectaMass = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            SNIIStarbyStarFeedback.Elements[i] = pCELibGetSNIIYieldsInterpolatedYieldMass(TableID,i,Mstar)
                                      +(Msw+Msn)*(Input.Elements[i]/Input.Mass);
            SNIIStarbyStarFeedback.EjectaMass += SNIIStarbyStarFeedback.Elements[i];
        }
        SNIIStarbyStarFeedback.Energy = pCELibGetSNIIYieldsInterpolatedEnergy(TableID,Mstar);

    } else if(Input.Metallicity >= MetallicitySNII[4]){
        double Msw = pCELibGetSNIIYieldsInterpolatedStellarWindMass(CELibSNIIYields_Metallicity-1,Mstar);
        double Msn = pCELibGetSNIIYieldsInterpolatedSuperNovaMass(CELibSNIIYields_Metallicity-1,Mstar);

        SNIIStarbyStarFeedback.EjectaMass = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            SNIIStarbyStarFeedback.Elements[i] = pCELibGetSNIIYieldsInterpolatedYieldMass(CELibSNIIYields_Metallicity-1,i,Mstar)
                                      +(Msw+Msn)*(Input.Elements[i]/Input.Mass);
            SNIIStarbyStarFeedback.EjectaMass += SNIIStarbyStarFeedback.Elements[i];
        }
        SNIIStarbyStarFeedback.Energy = pCELibGetSNIIYieldsInterpolatedEnergy(CELibSNIIYields_Metallicity-1,Mstar);
    } else {
        int IndexMetallicity = 0;
        for(int i=1;i<CELibSNIIYields_Metallicity;i++){
            if(MetallicitySNII[i] > Input.Metallicity){
                IndexMetallicity = i-1;
                break;
            }
        }

        double MetallicityEdges[] = {MetallicitySNII[IndexMetallicity], 
                                     MetallicitySNII[IndexMetallicity+1]};

        double MswEdges[] = {pCELibGetSNIIYieldsInterpolatedStellarWindMass(IndexMetallicity,Mstar),
                             pCELibGetSNIIYieldsInterpolatedStellarWindMass(IndexMetallicity+1,Mstar)};
        double MsnEdges[] = {pCELibGetSNIIYieldsInterpolatedSuperNovaMass(IndexMetallicity,Mstar),
                             pCELibGetSNIIYieldsInterpolatedSuperNovaMass(IndexMetallicity+1,Mstar)};

        double GradMsw = (MswEdges[1]-MswEdges[0])/(MetallicityEdges[1]-MetallicityEdges[0]);
        double GradMsn = (MsnEdges[1]-MsnEdges[0])/(MetallicityEdges[1]-MetallicityEdges[0]);

        double Msw = GradMsw*(Input.Metallicity-MetallicityEdges[0])+MswEdges[0];
        double Msn = GradMsn*(Input.Metallicity-MetallicityEdges[0])+MsnEdges[0];


        SNIIStarbyStarFeedback.EjectaMass = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            double GradYield = (pCELibGetSNIIYieldsInterpolatedYieldMass(IndexMetallicity+1,i,Mstar)-
                                        pCELibGetSNIIYieldsInterpolatedYieldMass(IndexMetallicity,i,Mstar))/
                                     (MetallicityEdges[1]-MetallicityEdges[0]);

            SNIIStarbyStarFeedback.Elements[i] = GradYield*(Input.Metallicity-MetallicityEdges[0]) + 
                                        pCELibGetSNIIYieldsInterpolatedYieldMass(IndexMetallicity,i,Mstar)
                                        + (Msw+Msn)*(Input.Elements[i]/Input.Mass);
            SNIIStarbyStarFeedback.EjectaMass += SNIIStarbyStarFeedback.Elements[i];
        }

        double GradEnergy = (pCELibGetSNIIYieldsInterpolatedEnergy(IndexMetallicity+1,Mstar)-
                                pCELibGetSNIIYieldsInterpolatedEnergy(IndexMetallicity,Mstar))/
                                 (MetallicityEdges[1]-MetallicityEdges[0]);
        SNIIStarbyStarFeedback.Energy = GradEnergy*(Input.Metallicity-MetallicityEdges[0]) + 
                                    pCELibGetSNIIYieldsInterpolatedEnergy(IndexMetallicity,Mstar);
    }

    SNIIStarbyStarFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor - SNIIStarbyStarFeedback.EjectaMass;
    SNIIStarbyStarFeedback.Energy *= Input.Mass*Input.MassConversionFactor;

    return SNIIStarbyStarFeedback;
}

/*!
 * This function retunrs the star-by-star feedback results. 
 */
struct CELibStructFeedbackStarbyStarOutput CELibGetFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input, const int Type){

    switch (Type){
        case CELibFeedbackType_SNII:
            return CELibGetSNIIFeedbackStarbyStar((struct CELibStructFeedbackStarbyStarInput){
                    .Mass = Input.Mass,
                    .Metallicity = Input.Metallicity,
                    .MassConversionFactor = Input.MassConversionFactor,
                    .Elements = Input.Elements,
                    });
            break;
        default:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            break;
    }

    //return ;
    return (struct CELibStructFeedbackStarbyStarOutput){
        .Energy = 0,
        .EjectaMass = 0,
        .RemnantMass = Input.Mass*Input.MassConversionFactor,
        .Elements[0] = 0,
        .Elements[1] = 0,
        .Elements[2] = 0,
        .Elements[3] = 0,
        .Elements[4] = 0,
        .Elements[5] = 0,
        .Elements[6] = 0,
        .Elements[7] = 0,
        .Elements[8] = 0,
        .Elements[9] = 0,
        .Elements[10] = 0,
        .Elements[11] = 0,
        .Elements[12] = 0,
        };
}