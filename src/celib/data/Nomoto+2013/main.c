#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "Astro.h"
#include "../../src/CELib.h"

#define InsertBoundaryMass (40.1)

/*
 * The ejecta mass of each element i can be expressed as
 *   E_ej,i = y_i + X_0,i E_sn + X_0,i E_sw
 *   where y_i is the yield of the i-th element, X_0,i is the mass fraction of
 *   i-th element of the projenitar star, E_sn and E_sw is the return masses via
 *   supernova and stellar wind.
 */


/*
 * In the case of OutputY == OFF, this code calculates the yield y assuming the solar
 * metal abundance pattern. The ejecta mass can be obtained by using the following
 * equation:
 *  E_ej,i = y_i + X_0,i (E_sn+E_sw), where E_sw = 0.
 *
 * On the other hand, in the case of OutputY == ON, this code uses Y_sn,i = y_i +X_0,i
 * E_sn. E_sw is also evaluated in this code. The ejecta mass is 
 *  E_ej,i = Y_sn,i + X_0,i E_sw.
 * Since E_sn = 0 in this mode, one can write 
 *  E_ej,i = Y_sn,i + X_0,i (E_sn + E_sw).
 *
 */
#define OutputY (OFF) 

#define NEntry_Zero  36
#define NEntry_Normal  27
#define NMetal  6   // 1+5
#define NElements  83 

#define CELibElements (13)

double Metal[NMetal];
double Mass[NMetal][NEntry_Zero];
double Energy[NMetal][NEntry_Zero];
double Mrem[NMetal][NEntry_Zero];
double Elements[NMetal][NEntry_Zero][NElements];
char Tags[NElements][MaxCharactersInLine];

static void ReadZero(FILE *fp){

    fscanf(fp,"%*s %le",&(Metal[0]));
    fprintf(stderr,"Z = %g\n",Metal[0]);
    
    fscanf(fp,"%*s");
    for(int k=0;k<NEntry_Zero;k++){
        fscanf(fp,"%le",Mass[0]+k);
    }

    fscanf(fp,"%*s");
    for(int k=0;k<NEntry_Zero;k++){
        fscanf(fp,"%le",Energy[0]+k);
    }

    fscanf(fp,"%*s");
    for(int k=0;k<NEntry_Zero;k++){
        fscanf(fp,"%le",Mrem[0]+k);
    }

    for(int i=0;i<NElements;i++){
        fscanf(fp,"%s",Tags[i]);
        fscanf(fp,"%*d");
        for(int k=0;k<NEntry_Zero;k++){
            fscanf(fp,"%le",Elements[0][k]+i);
        }
    }

    for(int k=0;k<NEntry_Zero;k++){
        fprintf(stderr,"[%d] %g %g %g\n",k,Mass[0][k],Mrem[0][k],Energy[0][k]);
    }

    for(int i=0;i<NElements;i++){
        fprintf(stderr,"%s %g\n",Tags[i],Elements[0][0][i]);
    }

    return;
}

static void ReadNormal(FILE *fp){

    for(int l=1;l<NMetal;l++){
        fscanf(fp,"%*s %le",&(Metal[l]));
        fprintf(stderr,"Z = %g\n",Metal[l]);
    
        fscanf(fp,"%*s");
        for(int k=0;k<NEntry_Normal;k++){
            fscanf(fp,"%le",Mass[l]+k);
        }

        fscanf(fp,"%*s");
        for(int k=0;k<NEntry_Normal;k++){
            fscanf(fp,"%le",Energy[l]+k);
        }

        fscanf(fp,"%*s");
        for(int k=0;k<NEntry_Normal;k++){
            fscanf(fp,"%le",Mrem[l]+k);
        }

        for(int i=0;i<NElements;i++){
            fscanf(fp,"%s",Tags[i]);
            fscanf(fp,"%*d");
            for(int k=0;k<NEntry_Normal;k++){
                fscanf(fp,"%le",Elements[l][k]+i);
            }
        }

        for(int k=0;k<NEntry_Normal;k++){
            fprintf(stderr,"[%d] %g %g %g\n",k,Mass[l][k],Mrem[l][k],Energy[l][k]);
        }

        for(int i=0;i<NElements;i++){
            fprintf(stderr,"%s %g %g\n",Tags[i],Elements[l][0][i],Elements[l][1][i]);
        }
    }

    return;
}

#define NEntry_Zero_HN 6
#define NEntry_HN 4
#define NMetal_HN  6   
#define NElements_HN  83 

double MetalHN[NMetal_HN];
double MassHN[NMetal_HN][NEntry_Zero_HN];
double EnergyHN[NMetal_HN][NEntry_Zero_HN];
double MremHN[NMetal_HN][NEntry_Zero_HN];
double ElementsHN[NMetal_HN][NEntry_Zero_HN][NElements_HN];
char TagsHN[NElements_HN][MaxCharactersInLine];

static void ReadHyperNovae(FILE *fp, const int mode){

    if(mode == 0){ // for Zero HyperNovae
        fscanf(fp,"%*s %le",&(MetalHN[0]));
        fprintf(stderr,"Z = %g\n",MetalHN[0]);

        fscanf(fp,"%*s");
        for(int k=0;k<NEntry_Zero_HN;k++){
            fscanf(fp,"%le",MassHN[0]+k);
        }

        fscanf(fp,"%*s");
        for(int k=0;k<NEntry_Zero_HN;k++){
            fscanf(fp,"%le",EnergyHN[0]+k);
        }

        fscanf(fp,"%*s");
        for(int k=0;k<NEntry_Zero_HN;k++){
            fscanf(fp,"%le",MremHN[0]+k);
        }

        for(int i=0;i<NElements_HN;i++){
            fscanf(fp,"%s",TagsHN[i]);
            fscanf(fp,"%*d");
            for(int k=0;k<NEntry_Zero_HN;k++){
                fscanf(fp,"%le",ElementsHN[0][k]+i);
            }
        }

        for(int k=0;k<NEntry_Zero_HN;k++){
            fprintf(stderr,"[%d] %g %g %g\n",k,MassHN[0][k],MremHN[0][k],EnergyHN[0][k]);
        }

        for(int i=0;i<NElements_HN;i++){
            fprintf(stderr,"%s %g %g\n",TagsHN[i],ElementsHN[0][0][i],ElementsHN[0][1][i]);
        }

    } else if(mode == 1){ // for Normal HyperNovae

        for(int l=1;l<NMetal_HN;l++){
            fscanf(fp,"%*s %le",&(MetalHN[l]));
            fprintf(stderr,"Z = %g\n",MetalHN[l]);
        
            fscanf(fp,"%*s");
            for(int k=0;k<NEntry_HN;k++){
                fscanf(fp,"%le",MassHN[l]+k);
            }

            fscanf(fp,"%*s");
            for(int k=0;k<NEntry_HN;k++){
                fscanf(fp,"%le",EnergyHN[l]+k);
            }

            fscanf(fp,"%*s");
            for(int k=0;k<NEntry_HN;k++){
                fscanf(fp,"%le",MremHN[l]+k);
            }

            for(int i=0;i<NElements_HN;i++){
                fscanf(fp,"%s",TagsHN[i]);
                fscanf(fp,"%*d");
                for(int k=0;k<NEntry_HN;k++){
                    fscanf(fp,"%le",ElementsHN[l][k]+i);
                }
            }

            for(int k=0;k<NEntry_HN;k++){
                fprintf(stderr,"[%d] %g %g %g\n",k,MassHN[l][k],MremHN[l][k],EnergyHN[l][k]);
            }

            for(int i=0;i<NElements_HN;i++){
                fprintf(stderr,"%s %g %g\n",TagsHN[i],ElementsHN[l][0][i],ElementsHN[l][1][i]);
            }
        }
    }

    return;
}



// Yields
#if 0
/*
 * The mass fractions of H(X), He(Y), and metals(Z) obtained from table 4 in
 * Asplund et al. ARAA, 2009.
 */
static double MassFraction[3] = {0.7381,0.2485,0.0134};
static double MassFractionFromTable[3];
// All data except Eu are obtained by Asplund et al. (2009).
// The Eu data is obtained from the table 2 in Anders & Grevesse (1989).
static double ElementAbundances[CELibElements] = {
    12.00, // H
    10.93, // He
    8.43, // C
    7.83, // N
    8.69, // O
    7.93, // Ne
    7.60, // Mg
    7.51, // Si
    7.12, // S
    6.34, // Ca
    7.50, // Fe
    6.22, // Ni
    0.51, // Eu // Derived from Anders & Grevesse 1989 table 2.
};
#elif 1
/*
 * The mass fractions of H(X), He(Y), and metals(Z) obtained from 
 * Anders & Grevesse 1989.
 */
static double MassFraction[3] = {0.70683,0.27431,0.01886};
static double MassFractionFromTable[3];
// All data except Eu are obtained from table 2 in Anders & Grevesse (1989).
static double ElementAbundances[CELibElements] = {
    12.00, // H
    10.99, // He
    8.56, // C
    8.05, // N
    8.93, // O
    8.09, // Ne
    7.58, // Mg
    7.55, // Si
    7.21, // S
    6.36, // Ca
    7.67, // Fe
    6.25, // Ni
    0.51, // Eu 
};
#else
/*
 * The mass fractions of H(X), He(Y), and metals(Z) obtained from 
 * Grevesse & Sauval 1998.
 */
static double MassFraction[3] = {0.735,0.248,0.017};
static double MassFractionFromTable[3];
// All data except Eu are obtained from table 2 in Grevesse & Sauval(1998).
static double ElementAbundances[CELibElements] = {
    12.00, // H
    10.93, // He
    8.52, // C
    7.92, // N
    8.83, // O
    8.08, // Ne
    7.58, // Mg
    7.55, // Si
    7.33, // S
    6.36, // Ca
    7.50, // Fe
    6.25, // Ni
    0.51, // Eu 
};
#endif

static double ElementZ[CELibElements] = {
     1, // H
     4, // He
    12, // C
    14, // N
    16, // O
    20, // Ne
    24, // Mg
    28, // Si
    32, // S
    40, // Ca
    56, // Fe
    58, // Ni, 58.69
    153, // Eu, 151,153, 153 is more abundant
};


static double ElementNumberDensities[CELibElements];
static double ElementMasses[CELibElements];
static double ElementNormalizedMasses[CELibElements];

static bool first = true;
/*
 * In this function, the mass fractions of the material with the solar
 * abundances is calculated. The mass fraction is little bit different from the
 * literature (Table 4, Asplund et al. ARAA. 2009).
 */
void CELibInitSolarAbundances(void){

    if(first == false) 
        return;

    double Mass = 0.e0;
    for(int i=0;i<CELibElements;i++){
        ElementNumberDensities[i] = pow(10.0,ElementAbundances[i]-12);
        ElementMasses[i] = ElementZ[i]*ElementNumberDensities[i];
        Mass += ElementMasses[i];
    }


    MassFractionFromTable[2] = 0.e0;
    for(int i=0;i<CELibElements;i++){
        ElementNormalizedMasses[i] = ElementMasses[i]/Mass;
        if(i>1)
            MassFractionFromTable[2] += ElementNormalizedMasses[i];
    }
    MassFractionFromTable[0] = ElementNormalizedMasses[0];
    MassFractionFromTable[1] = ElementNormalizedMasses[1];

    fprintf(stderr,"X,Y,Z = %g %g %g | %g\n",
            MassFractionFromTable[0],
            MassFractionFromTable[1],
            MassFractionFromTable[2],
            MassFractionFromTable[0]+MassFractionFromTable[1]+MassFractionFromTable[2]);

    first = false;

    return ;
}

/*
 * This function returns the mass of each element. The total metallicity is set
 * to "Metallicity". The abundance pattern is based on Asploud et al. Annu.
 * Rev. Astron. Astrophys. 2009. 47:481–522 Rescale factors are multipled to H
 * and Zs abundances.
 */
void CELibSetMetallicityWithSolarAbundancePattern(const double Mass, double Elements[restrict], const double Metallicity){ 

    if(first == true)
        CELibInitSolarAbundances();

    double Scale = Metallicity*(MassFractionFromTable[0]+MassFractionFromTable[1])/((1-Metallicity)*MassFractionFromTable[2]);
    Elements[CELibYield_H]  = MassFractionFromTable[0]*Mass;
    Elements[CELibYield_He] = MassFractionFromTable[1]*Mass;
    Elements[CELibYield_C]  = Scale*Mass*ElementNormalizedMasses[CELibYield_C];
    Elements[CELibYield_N]  = Scale*Mass*ElementNormalizedMasses[CELibYield_N];
    Elements[CELibYield_O]  = Scale*Mass*ElementNormalizedMasses[CELibYield_O];
    Elements[CELibYield_Ne] = Scale*Mass*ElementNormalizedMasses[CELibYield_Ne];
    Elements[CELibYield_Mg] = Scale*Mass*ElementNormalizedMasses[CELibYield_Mg];
    Elements[CELibYield_Si] = Scale*Mass*ElementNormalizedMasses[CELibYield_Si];
    Elements[CELibYield_S]  = Scale*Mass*ElementNormalizedMasses[CELibYield_S];
    Elements[CELibYield_Ca] = Scale*Mass*ElementNormalizedMasses[CELibYield_Ca];
    Elements[CELibYield_Fe] = Scale*Mass*ElementNormalizedMasses[CELibYield_Fe];
    Elements[CELibYield_Ni] = Scale*Mass*ElementNormalizedMasses[CELibYield_Ni];
    //Elements[CELibYield_Eu] = Scale*Mass*ElementNormalizedMasses[CELibYield_Eu];
    Elements[CELibYield_Eu] = 0.e0;

    double Mtotal = 0.e0;
    for(int i=0;i<CELibElements;i++)
        Mtotal += Elements[i];
    double InvMtotal = Mass/Mtotal;
    for(int i=0;i<CELibElements;i++)
        Elements[i] *= InvMtotal;

#if 0
    fprintf(stderr," Check %g %g\n",Mass,
            Elements[CELibYield_H]  +
            Elements[CELibYield_He] +
            Elements[CELibYield_C]  +
            Elements[CELibYield_N]  +
            Elements[CELibYield_O]  +
            Elements[CELibYield_Ne] +
            Elements[CELibYield_Mg] +
            Elements[CELibYield_Si] +
            Elements[CELibYield_S]  +
            Elements[CELibYield_Ca] +
            Elements[CELibYield_Fe] +
            Elements[CELibYield_Ni] +
            Elements[CELibYield_Eu]);
#endif

    return ;
}


/*
 * Cosmic helium abundance obtained by Planck paper XVI(Planck Collaboration,
 * arXiv:1303.5076). This value is found in Table 2 of this paper.
 */
#define CELibCosmicHeliumAbundance  (0.247695)


/*
 * This function returns the mass of each element based on the primodrial
 * abandance pattern.
 */
void CELibSetPrimordialMetallicity(const double Mass, double Elements[restrict]){

    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] = 0.e0;
    }
    Elements[CELibYield_H] = (1.0-CELibCosmicHeliumAbundance)*Mass;
    Elements[CELibYield_He] = CELibCosmicHeliumAbundance*Mass;

    return ;
}



static bool FlagMetal[] = {true,true,true,false,true,true};
static double FlagMass[] = {6,7,10,13,15,20,24,30,40};
static char TagNames[][MaxCharactersInLine] = {"H","He","C","N","O","Ne","Mg","Si","S","Ca","Fe","Ni","Eu"};
static bool FlagMatchZero[] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true};
static bool FlagMatch[] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    true,true,true,true,true,true,true};
static int Length_Zero,Length_Normal;


#define MassMergin (0.05)
static bool MatchingMass(const double mass){
    size_t NM = sizeof(FlagMass)/sizeof(double);
    for(int i=0;i<NM;i++){
        //if(mass == FlagMass[i]){
            //return true;
        //}
        if((FlagMass[i]*(1.0-MassMergin) < mass)&&(mass < FlagMass[i]*(1.0+MassMergin))){
            return true;
        }
    }
    return false;
}

static double MatchingElement(const int MetalID, const int MassID, const char ElementName[]){
    double Mass = 0.e0;
    for(int i=0;i<NElements;i++){
        int ret = strcmp(ElementName,Tags[i]);
        if(ret == 0){
            Mass += Elements[MetalID][MassID][i];
        }
    }
    return Mass;
}

static double CountAllMatchingElement(const int MetalID, const int MassID){

    bool Flags[NElements];
    for(int k=0;k<NElements;k++) Flags[k] = true;
    
    for(int k=0;k<CELibElements;k++){
        for(int i=0;i<NElements;i++){
            int ret = strcmp(TagNames[k],Tags[i]);
            if(ret == 0){
                Flags[i] = false;
            }
        }
    }

    double Mass = 0.e0;
    for(int i=0;i<NElements;i++){
        if(Flags[i] == false){
            Mass += Elements[MetalID][MassID][i];
        }
    }

    return Mass;
}

static double MatchingNoUseElement(const int MetalID, const int MassID){

    bool Flags[NElements];
    for(int k=0;k<NElements;k++) Flags[k] = true;
    
    for(int k=0;k<CELibElements;k++){
        for(int i=0;i<NElements;i++){
            int ret = strcmp(TagNames[k],Tags[i]);
            if(ret == 0){
                Flags[i] = false;
            }
        }
    }

    double Mass = 0.e0;
    for(int i=0;i<NElements;i++){
        if(Flags[i] == true){
            Mass += Elements[MetalID][MassID][i];
        }
    }
    return Mass;
}

static double CountAllElement(const int MetalID, const int MassID){

    double Mass = 0.e0;
    for(int i=0;i<NElements;i++){
        Mass += Elements[MetalID][MassID][i];
    }
    return Mass;
}

static void WriteData(FILE *fp){
    
    for(int i=0;i<NMetal;i++){  
        if(FlagMetal[i]){
            if(i == 0){
                int count_output = 0;
                for(int k=0;k<NEntry_Zero;k++){
                    if(FlagMatchZero[k]){
                        fprintf(fp,"%g\n",Metal[i]);
                        fprintf(fp,"%g\n",Mass[i][k]);
                        fprintf(fp,"%g\n",Mrem[i][k]);
                        fprintf(fp,"%g\n",Mass[i][k]-(CountAllMatchingElement(i,k)+Mrem[i][k]));
                        fprintf(fp,"%g\n",Energy[i][k]);
                        double mass_total = 0.e0;
                        for(int l=0;l<CELibElements;l++){
                            fprintf(fp,"%g\n",MatchingElement(i,k,TagNames[l]));
                            mass_total += MatchingElement(i,k,TagNames[l]);
                        }
                        fprintf(fp,"\n");
                        count_output ++;
                    }
                }
                Length_Zero = count_output;
            } else {
                int count_output = 0;
                for(int k=0;k<NEntry_Normal;k++){
                    if(FlagMatch[k]){
                        fprintf(fp,"%g\n",Metal[i]);
                        fprintf(fp,"%g\n",Mass[i][k]);
                        fprintf(fp,"%g\n",Mrem[i][k]);
                        fprintf(fp,"%g\n",Mass[i][k]-(CountAllMatchingElement(i,k)+Mrem[i][k]));
                        fprintf(fp,"%g\n",Energy[i][k]);
                        double mass_total = 0.e0;
                        for(int l=0;l<CELibElements;l++){
                            fprintf(fp,"%g\n",MatchingElement(i,k,TagNames[l]));
                            mass_total += MatchingElement(i,k,TagNames[l]);
                        }
                        fprintf(fp,"\n");
                        count_output ++;
                    }
                }
                Length_Normal = count_output;
            }
        }
    }

    return ;
}


static double CountAllMatchingElementHN(const int MetalID, const int MassID){

    bool Flags[NElements_HN];
    for(int k=0;k<NElements_HN;k++) Flags[k] = true;
    
    for(int k=0;k<CELibElements;k++){
        for(int i=0;i<NElements_HN;i++){
            int ret = strcmp(TagNames[k],TagsHN[i]);
            if(ret == 0){
                Flags[i] = false;
            }
        }
    }

    double Mass = 0.e0;
    for(int i=0;i<NElements_HN;i++){
        if(Flags[i] == false){
            Mass += ElementsHN[MetalID][MassID][i];
        }
    }

    return Mass;
}

static double MatchingElementHN(const int MetalID, const int MassID, const char ElementName[]){
    double Mass = 0.e0;
    for(int i=0;i<NElements_HN;i++){
        int ret = strcmp(ElementName,TagsHN[i]);
        if(ret == 0){
            Mass += ElementsHN[MetalID][MassID][i];
        }
    }
    return Mass;
}

static void WriteDataHN(FILE *fp){
    
    for(int i=0;i<NMetal_HN;i++){  
        if(FlagMetal[i]){
            if(i == 0){
                int count_output = 0;
                for(int k=0;k<NEntry_Zero_HN;k++){
                    fprintf(fp,"%g\n",MetalHN[i]);
                    fprintf(fp,"%g\n",MassHN[i][k]);
                    fprintf(fp,"%g\n",MremHN[i][k]);
                    fprintf(fp,"%g\n",MassHN[i][k]-(CountAllMatchingElementHN(i,k)+MremHN[i][k]));
                    fprintf(fp,"%g\n",EnergyHN[i][k]);
                    double mass_total = 0.e0;
                    for(int l=0;l<CELibElements;l++){
                        fprintf(fp,"%g\n",MatchingElementHN(i,k,TagNames[l]));
                        mass_total += MatchingElementHN(i,k,TagNames[l]);
                    }
                    fprintf(fp,"\n");
                    count_output ++;
                }
                Length_Zero = count_output;
            } else {
                int count_output = 0;
                for(int k=0;k<NEntry_HN;k++){
                    fprintf(fp,"%g\n",MetalHN[i]);
                    fprintf(fp,"%g\n",MassHN[i][k]);
                    fprintf(fp,"%g\n",MremHN[i][k]);
                    fprintf(fp,"%g\n",MassHN[i][k]-(CountAllMatchingElementHN(i,k)+MremHN[i][k]));
                    fprintf(fp,"%g\n",EnergyHN[i][k]);
                    double mass_total = 0.e0;
                    for(int l=0;l<CELibElements;l++){
                        fprintf(fp,"%g\n",MatchingElementHN(i,k,TagNames[l]));
                        mass_total += MatchingElementHN(i,k,TagNames[l]);
                    }
                    // Check others
                    fprintf(fp,"\n");
                    count_output ++;
                }
                Length_Normal = count_output;
            }
        }
    }

    return ;
}


static char CELibSNIILabel[][MaxCharactersInLine] = {
    "CELibYield_H","CELibYield_He","CELibYield_C",
    "CELibYield_N","CELibYield_O","CELibYield_Ne",
    "CELibYield_Mg","CELibYield_Si","CELibYield_S",
    "CELibYield_Ca","CELibYield_Fe","CELibYield_Ni",
    "CELibYield_Eu"};

static void WriteStructure(FILE *fp){

    fprintf(fp,"#define CELibSNIIYields_Metallicity_N13 (%d)\n",NMetal-1);

    int counter_mass = 0;
    for(int k=0;k<NEntry_Zero;k++){
        //if(MatchingMass(Mass[0][k])){
        if(FlagMatchZero[k]){
            counter_mass ++;
        }
    }
    fprintf(fp,"#define CELibSNIIYields_Mass_N13 (%d)\n",counter_mass+2+1); 
    // +2 is used for the outer boundary.
    // +1 is used for the boundary for 40 Msun.

    fprintf(fp,"static struct CELibStructSNIIYields CELibSNIIYieldsN13[CELibSNIIYields_Metallicity_N13][CELibSNIIYields_Mass_N13] = \{\n");
    fprintf(fp,"\t\{\n");

    ////// For Z==0
    double EmptyMassZero[] = {6,10};
    for(int i=0;i<2;i++){
        fprintf(fp,"\t\{\n");
        fprintf(fp,"\t.Metallicity = %g,\n",Metal[0]);
        fprintf(fp,"\t.Mass = %g,\n",EmptyMassZero[i]);
        fprintf(fp,"\t.Mr = %g,\n",EmptyMassZero[i]);
        fprintf(fp,"\t.Msw = %g,\n",0.0);
        fprintf(fp,"\t.Energy = %g,\n",0.0);
        for(int l=0;l<CELibElements;l++){
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],0.0);
        }
        fprintf(fp,"\t},\n");
    }


    int count_output = 0;
    for(int k=0;k<NEntry_Zero;k++){
        if((Mass[0][k] == 140.0)&&(Energy[0][k] == 1.0)){
            Mass[0][k] = 139.5;
        }
        if(FlagMatchZero[k]){
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",Metal[0]);
            fprintf(fp,"\t.Mass = %g,\n",Mass[0][k]);
            if(Mass[0][k] <= 10.e0){
                fprintf(fp,"\t.Mr = %g,\n",Mass[0][k]);
            } else {
                fprintf(fp,"\t.Mr = %g,\n",Mrem[0][k]);
            }
#if OutputY // Esn = 0 and Esw has a certain value.
            if(Mass[0][k] <= 10.e0){
                fprintf(fp,"\t.Msw = %g,\n",0.e0);
                fprintf(fp,"\t.Msn = %g,\n",0.e0);
            } else {
                fprintf(fp,"\t.Msw = %g,\n",0.e0);
                fprintf(fp,"\t.Msn = %g,\n",Mass[0][k]-(CountAllMatchingElement(0,k)+Mrem[0][k]));
            }
#else       // Esn has a certain value while Esw = 0
            if(Mass[0][k] <= 10.e0){
                fprintf(fp,"\t.Msw = %g,\n",0.e0);
                fprintf(fp,"\t.Msn = %g,\n",0.e0);
            } else {
                fprintf(fp,"\t.Msw = %g,\n",Mass[0][k]-(CountAllMatchingElement(0,k)+Mrem[0][k]));
                fprintf(fp,"\t.Msn = %g,\n",CountAllMatchingElement(0,k));
                //fprintf(fp,"\t.Msn = %g,\n",Mass[i][k]-Mrem[i][k]);
            }
#endif
            fprintf(fp,"\t.Energy = %g,\n",Energy[0][k]*1.e+51);
#if OutputY // Esn = 0 and Esw has a certain value.
            double mass_total = 0.e0;
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],MatchingElement(0,k,TagNames[l]));
                mass_total += MatchingElement(0,k,TagNames[l]);
            }
#else       // Esn has a certain value while Esw = 0

            // Get solar abundance pattern
            double X0[CELibElements];
            CELibSetPrimordialMetallicity(1.0,X0);

            double Msw = Mass[0][k]-(CountAllMatchingElement(0,k)+Mrem[0][k]);
            double Eejecta_i[CELibElements];
            for(int l=0;l<CELibElements;l++){
                Eejecta_i[l] = MatchingElement(0,k,TagNames[l]) + X0[l]*Msw;
            }

            double Eejecta = Mass[0][k]-Mrem[0][k];
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],Eejecta_i[l]-X0[l]*Eejecta);
            }

#endif
            fprintf(fp,"\t},\n");
            count_output ++;
        }
        if(Mass[0][k] == 40.0){ // InsertBoundaryMass 
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",Metal[0]);
            fprintf(fp,"\t.Mass = %g,\n",InsertBoundaryMass);
            fprintf(fp,"\t.Mr = %g,\n",Mrem[0][k]);
            fprintf(fp,"\t.Msw = %g,\n",Mass[0][k]-(CountAllMatchingElement(0,k)+Mrem[0][k]));
            fprintf(fp,"\t.Msn = %g,\n",CountAllMatchingElement(0,k));
            fprintf(fp,"\t.Energy = %g,\n",Energy[0][k]*1.e+51);
            double X0[CELibElements];
            CELibSetPrimordialMetallicity(1.0,X0);

            double Msw = Mass[0][k]-(CountAllMatchingElement(0,k)+Mrem[0][k]);
            double Eejecta_i[CELibElements];
            for(int l=0;l<CELibElements;l++){
                Eejecta_i[l] = MatchingElement(0,k,TagNames[l]) + X0[l]*Msw;
            }

            double Eejecta = Mass[0][k]-Mrem[0][k];
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],Eejecta_i[l]-X0[l]*Eejecta);
            }
            fprintf(fp,"\t},\n");
            count_output ++;
        }
    }
    fprintf(fp,"\t},\n");
    fprintf(fp,"\t\{\n");

    ////// For Z>0
    for(int i=1;i<NMetal;i++){  
        if(Metal[i] == 0.008) continue;
        double EmptyMassNormal[] = {6,10,11};
        for(int ll=0;ll<3;ll++){
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",Metal[i]);
            fprintf(fp,"\t.Mass = %g,\n",EmptyMassNormal[ll]);
            fprintf(fp,"\t.Mr = %g,\n",EmptyMassNormal[ll]);
            fprintf(fp,"\t.Msw = %g,\n",0.0);
            fprintf(fp,"\t.Energy = %g,\n",0.0);
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],0.0);
            }
            fprintf(fp,"\t},\n");
        }

        for(int k=0;k<NEntry_Normal;k++){
            if(FlagMatch[k]){
                fprintf(fp,"\t\{\n");
                fprintf(fp,"\t.Metallicity = %g,\n",Metal[i]);
                fprintf(fp,"\t.Mass = %g,\n",Mass[i][k]);
                if(Mass[i][k] <= 10.e0){
                    fprintf(fp,"\t.Mr = %g,\n",Mass[i][k]);
                } else {
                    fprintf(fp,"\t.Mr = %g,\n",Mrem[i][k]);
                }
#if OutputY // Esn = 0 and Esw has a certain value.
                if(Mass[i][k] <= 10.e0){
                    fprintf(fp,"\t.Msw = %g,\n",0.e0);
                    fprintf(fp,"\t.Msn = %g,\n",0.e0);
                } else {
                    fprintf(fp,"\t.Msw = %g,\n",0.e0);
                    fprintf(fp,"\t.Msn = %g,\n",Mass[i][k]-(CountAllMatchingElement(i,k)+Mrem[i][k]));
                }
#else       // Esn has a certain value while Esw = 0
                if(Mass[i][k] <= 10.e0){
                    fprintf(fp,"\t.Msw = %g,\n",0.e0);
                    fprintf(fp,"\t.Msn = %g,\n",0.e0);
                } else {
                    fprintf(fp,"\t.Msw = %g,\n",Mass[i][k]-(CountAllMatchingElement(i,k)+Mrem[i][k]));
                    fprintf(fp,"\t.Msn = %g,\n",CountAllMatchingElement(i,k));
                }
#endif
                fprintf(fp,"\t.Energy = %g,\n",Energy[i][k]*1.e+51);
#if OutputY // Esn = 0 and Esw has a certain value.
                double mass_total = 0.e0;
                for(int l=0;l<CELibElements;l++){
                    fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],MatchingElement(i,k,TagNames[l]));
                    mass_total += MatchingElement(i,k,TagNames[l]);
                }
#else       // Esn has a certain value while Esw = 0

                // Get solar abundance pattern
                double X0[CELibElements];
                CELibSetMetallicityWithSolarAbundancePattern(1.0,X0,Metal[i]);

                double Msw = Mass[i][k]-(CountAllMatchingElement(i,k)+Mrem[i][k]);
                double Eejecta_i[CELibElements];
                for(int l=0;l<CELibElements;l++){
                    Eejecta_i[l] = MatchingElement(i,k,TagNames[l]) + X0[l]*Msw;
                }

                double Eejecta = Mass[i][k]-Mrem[i][k];
                for(int l=0;l<CELibElements;l++){
                    fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],Eejecta_i[l]-X0[l]*Eejecta);
                }

#endif
                fprintf(fp,"\t},\n");
                count_output ++;
            }
        }

        double EmptyMassNormal2[] = {InsertBoundaryMass,100,139.5,140,150,170,200,270,300};
        for(int ll=0;ll<9;ll++){
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",Metal[i]);
            fprintf(fp,"\t.Mass = %g,\n",EmptyMassNormal2[ll]);
            fprintf(fp,"\t.Mr = %g,\n",EmptyMassNormal2[ll]);
            fprintf(fp,"\t.Msw = %g,\n",0.0);
            fprintf(fp,"\t.Energy = %g,\n",0.0);
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],0.0);
            }
            fprintf(fp,"\t},\n");
        }
        fprintf(fp,"\t},\n");
        if(i!=NMetal-1){
            fprintf(fp,"\t\{\n");
        } else {
            fprintf(fp,"};\n");
        }
    }

    return ;
}

static void WriteStructureHN(FILE *fp){

    fprintf(fp,"#define CELibSNIIYields_Metallicity_N13HN (%d)\n",NMetal_HN-1);
    fprintf(fp,"#define CELibSNIIYields_Mass_N13HN (%d)\n",NEntry_Zero_HN);
    fprintf(fp,"#define CELibSNIIYields_Mass_N13HN_Zero (%d)\n",NEntry_Zero_HN);
    fprintf(fp,"#define CELibSNIIYields_Mass_N13HN_Normal (%d)\n",NEntry_HN);

    fprintf(fp,"static struct CELibStructSNIIYields CELibSNIIYieldsN13HN[CELibSNIIYields_Metallicity_N13HN][CELibSNIIYields_Mass_N13HN] = \{\n");
    fprintf(fp,"\t\{\n");

    ////// For Z==0
    for(int k=0;k<NEntry_Zero_HN;k++){
        fprintf(fp,"\t\{\n");
        fprintf(fp,"\t.Metallicity = %g,\n",MetalHN[0]);
        fprintf(fp,"\t.Mass = %g,\n",MassHN[0][k]);
        fprintf(fp,"\t.Mr = %g,\n",MremHN[0][k]);
#if OutputY // Esn = 0 and Esw has a certain value.
        fprintf(fp,"\t.Msw = %g,\n",0.e0);
        fprintf(fp,"\t.Msn = %g,\n",MassHN[0][k]-(CountAllMatchingElementHN(0,k)+MremHN[0][k]));
#else       // Esn has a certain value while Esw = 0
        fprintf(fp,"\t.Msw = %g,\n",MassHN[0][k]-(CountAllMatchingElementHN(0,k)+MremHN[0][k]));
        fprintf(fp,"\t.Msn = %g,\n",CountAllMatchingElementHN(0,k));
#endif
        fprintf(fp,"\t.Energy = %g,\n",EnergyHN[0][k]*1.e+51);

#if OutputY // Esn = 0 and Esw has a certain value.
        double mass_total = 0.e0;
        for(int l=0;l<CELibElements;l++){
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],MatchingElementHN(0,k,TagNames[l]));
            mass_total += MatchingElementHN(0,k,TagNames[l]);
        }
#else       // Esn has a certain value while Esw = 0

        // Get solar abundance pattern
        double X0[CELibElements];
        CELibSetPrimordialMetallicity(1.0,X0);

        double Msw = MassHN[0][k]-(CountAllMatchingElementHN(0,k)+MremHN[0][k]);
        double Eejecta_i[CELibElements];
        for(int l=0;l<CELibElements;l++){
            Eejecta_i[l] = MatchingElementHN(0,k,TagNames[l]) + X0[l]*Msw;
        }

        double Eejecta = MassHN[0][k]-MremHN[0][k];
        for(int l=0;l<CELibElements;l++){
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],Eejecta_i[l]-X0[l]*Eejecta);
        }
#endif
        fprintf(fp,"\t},\n");
    }
    fprintf(fp,"\t},\n");
    fprintf(fp,"\t\{\n");

    ////// For Z>0
    for(int i=1;i<NMetal_HN;i++){  
        if(MetalHN[i] == 0.008) continue;
        for(int k=0;k<NEntry_HN;k++){
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",MetalHN[i]);
            fprintf(fp,"\t.Mass = %g,\n",MassHN[i][k]);
            fprintf(fp,"\t.Mr = %g,\n",MremHN[i][k]);
#if OutputY // Esn = 0 and Esw has a certain value.
            fprintf(fp,"\t.Msw = %g,\n",0.e0);
            fprintf(fp,"\t.Msn = %g,\n",MassHN[i][k]-(CountAllMatchingElementHN(i,k)+MremHN[i][k]));
#else       // Esn has a certain value while Esw = 0
            fprintf(fp,"\t.Msw = %g,\n",MassHN[i][k]-(CountAllMatchingElementHN(i,k)+MremHN[i][k]));
            fprintf(fp,"\t.Msn = %g,\n",CountAllMatchingElementHN(i,k));
#endif
            fprintf(fp,"\t.Energy = %g,\n",EnergyHN[i][k]*1.e+51);
#if OutputY // Esn = 0 and Esw has a certain value.
            double mass_total = 0.e0;
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],MatchingElementHN(i,k,TagNames[l]));
                mass_total += MatchingElementHN(i,k,TagNames[l]);
            }
#else       // Esn has a certain value while Esw = 0

            // Get solar abundance pattern
            double X0[CELibElements];
            CELibSetMetallicityWithSolarAbundancePattern(1.0,X0,Metal[i]);

            double Msw = MassHN[i][k]-(CountAllMatchingElementHN(i,k)+MremHN[i][k]);
            double Eejecta_i[CELibElements];
            for(int l=0;l<CELibElements;l++){
                Eejecta_i[l] = MatchingElementHN(i,k,TagNames[l]) + X0[l]*Msw;
            }

            double Eejecta = MassHN[i][k]-MremHN[i][k];
            for(int l=0;l<CELibElements;l++){
                fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[l],Eejecta_i[l]-X0[l]*Eejecta);
            }

#endif
            fprintf(fp,"\t},\n");
            //fprintf(fp,"\t\{\n");
        }
        fprintf(fp,"\t},\n");
        if(i!=NMetal_HN-1){
            fprintf(fp,"\t\{\n");
        } else {
            fprintf(fp,"};\n");
        }
    }

    return ;
}

static void DumpAllData(void){
    FILE *fp;
    FileOpen(fp,"./dump.dat","w");

    fprintf(fp,"Z= %g\n",Metal[0]);
    fprintf(fp,"M= ");
    for(int i=0;i<NEntry_Zero;i++)
        fprintf(fp,"%g ",Mass[0][i]);
    fprintf(fp,"\n");

    fprintf(fp,"E= ");
    for(int i=0;i<NEntry_Zero;i++)
        fprintf(fp,"%g ",Energy[0][i]);
    fprintf(fp,"\n");

    fprintf(fp,"Mrem= ");
    for(int i=0;i<NEntry_Zero;i++)
        fprintf(fp,"%g ",Mrem[0][i]);
    fprintf(fp,"\n");

    for(int k=0;k<NElements;k++){
        fprintf(fp,"%s ",Tags[k]);
        for(int i=0;i<NEntry_Zero;i++){
            fprintf(fp,"%g ",Elements[0][i][k]);
        }
        fprintf(fp,"\n");
    }

    
    for(int l=1;l<NMetal;l++){

        fprintf(fp,"Z= %g\n",Metal[l]);
        fprintf(fp,"M= ");
        for(int i=0;i<NEntry_Normal;i++)
            fprintf(fp,"%g ",Mass[l][i]);
        fprintf(fp,"\n");

        fprintf(fp,"E= ");
        for(int i=0;i<NEntry_Normal;i++)
            fprintf(fp,"%g ",Energy[l][i]);
        fprintf(fp,"\n");

        fprintf(fp,"Mrem= ");
        for(int i=0;i<NEntry_Normal;i++)
            fprintf(fp,"%g ",Mrem[l][i]);
        fprintf(fp,"\n");

        for(int k=0;k<NElements;k++){
            fprintf(fp,"%s ",Tags[k]);
            for(int i=0;i<NEntry_Normal;i++){
                fprintf(fp,"%g ",Elements[l][i][k]);
            }
            fprintf(fp,"\n");
        }
    }
    
    fclose(fp);


    int counter = 0;
    for(int i=0;i<NEntry_Zero;i++){
        if(Energy[0][i] > 0.0){
            fprintf(stderr,"Mass %g, Energy = %g\n",Mass[0][i],Energy[0][i]);
            counter ++;
        }
    }
    fprintf(stderr,"Total E>0 for Z=0 is %d\n",counter);
    counter = 0;
    for(int i=0;i<NEntry_Normal;i++){
        if(Energy[1][i] > 0.0){
            fprintf(stderr,"Mass %g, Energy = %g\n",Mass[1][i],Energy[1][i]);
            counter ++;
        }
    }
    fprintf(stderr,"Total E>0 for Z>0 is %d\n",counter);

    return;
}


double ElementsSolar[NMetal][NEntry_Zero][CELibElements];
static void ConvertYield(void){

    for(int i=0;i<NMetal;i++){
        double Z = Metal[i];

        int MaxEntry = i==0?NEntry_Zero:NEntry_Normal;
        int count_output = 0;
        //for(int k=0;k<MaxEntry;k++){
        //for(int j=0;j<NEntry;j++){
        for(int j=0;j<MaxEntry;j++){
            //double CurrentMass = Mass[i][j];
            double CurrentMass = Mass[i][j]-(CountAllMatchingElement(i,j)+Mrem[i][j]);
            double Elements[CELibElements];
            if(i==0){
                CELibSetPrimordialMetallicity(CurrentMass,Elements);
            } else {
                CELibSetMetallicityWithSolarAbundancePattern(CurrentMass,Elements,Z);
            }
            //CELibSetMetallicityWithSolarAbundancePattern(CurrentMass,Elements,Z);
            // fprintf(stderr," Check Y = %g %g Z = %g, Mass = %g\n",Elements[CELibYield_He],Elements[CELibYield_He]/CurrentMass,Z,Mass[i][j]);
            for(int k=0;k<CELibElements;k++){
                ElementsSolar[i][j][k] = Elements[k];
            }
        }
    }

    return ;
}


double ElementsSolarHN[NMetal_HN][NEntry_Zero_HN][CELibElements];
static void ConvertYieldHN(void){

    for(int i=0;i<NMetal_HN;i++){
        double Z = MetalHN[i];

        int MaxEntry = i==0?NEntry_Zero_HN:NEntry_HN;
        int count_output = 0;
        for(int j=0;j<MaxEntry;j++){
            double CurrentMass = MassHN[i][j]-(CountAllMatchingElementHN(i,j)+MremHN[i][j]);
            double Elements[CELibElements];
            if(i==0){
                CELibSetPrimordialMetallicity(CurrentMass,Elements);
            } else {
                CELibSetMetallicityWithSolarAbundancePattern(CurrentMass,Elements,Z);
            }
            for(int k=0;k<CELibElements;k++){
                ElementsSolarHN[i][j][k] = Elements[k];
            }
        }
    }

    return ;
}

static void WriteYield(void){

    for(int i=0;i<NMetal;i++){

        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"Yields.%02d",i);
        FileOpen(fp,fname,"w");

        double EejectaLog_i[CELibElements]; 
        int MaxEntry = i==0?NEntry_Zero:NEntry_Normal;
        int count_output = 0;
        for(int k=0;k<MaxEntry;k++){
            if(MatchingMass(Mass[i][k])){

                // Get solar abundance pattern
                double X0[CELibElements];
                CELibSetMetallicityWithSolarAbundancePattern(1.0,X0,Metal[i]);

                double Msw = Mass[i][k]-(CountAllMatchingElement(i,k)+Mrem[i][k]);
                double Eejecta_i[CELibElements];
                for(int l=0;l<CELibElements;l++){
                    Eejecta_i[l] = MatchingElement(i,k,TagNames[l]) + X0[l]*Msw;
                }


                double Eejecta = Mass[i][k]-Mrem[i][k];
                for(int l=0;l<CELibElements;l++){
                    EejectaLog_i[l] = Eejecta_i[l]-X0[l]*Eejecta;
                }


                if(Mass[i][k]>=10){
                    fprintf(fp,"%g %g %g %g %g %g %g\n",Mass[i][k],Mrem[i][k]/Mass[i][k],
                            (EejectaLog_i[CELibYield_H])/Mass[i][k],
                            (EejectaLog_i[CELibYield_He])/Mass[i][k],
                            (EejectaLog_i[CELibYield_C])/Mass[i][k],
                            (EejectaLog_i[CELibYield_O])/Mass[i][k],
                            ((EejectaLog_i[CELibYield_Si])+EejectaLog_i[CELibYield_S]+EejectaLog_i[CELibYield_Ca]+EejectaLog_i[CELibYield_Fe])/Mass[i][k]);
                }
            }

        }
        fclose(fp);
    }

    return ;
}




int main(int argc, char **argv){

    FILE *fp;
    FileOpen(fp,"./YIELD_CK13.DAT","r");

    // Read type II SNe data
    ReadZero(fp);
    ReadNormal(fp);
    ConvertYield();

    // Read hypernova data
    ReadHyperNovae(fp,0);
    ReadHyperNovae(fp,1);
    ConvertYieldHN();

    fclose(fp);

    // Write yields data
    FileOpen(fp,"./Nomoto+2013.dat","w");
    WriteData(fp);
    WriteDataHN(fp);
    fclose(fp);

    // Write a type II SNe structure
    FileOpen(fp,"./Nomoto+2013_Struct.h","w");
    WriteStructure(fp);
    fclose(fp);

    // Write a hyper nova structure
    FileOpen(fp,"./Nomoto+2013HN_Struct.h","w");
    WriteStructureHN(fp);
    fclose(fp);

    // dump all data
    DumpAllData();

    // Write all yields data.
    WriteYield();

    return EXIT_SUCCESS;
}
