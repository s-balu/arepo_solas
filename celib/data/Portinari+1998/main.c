#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "Astro.h"
#include "../../src/CELib.h"


#define CELibElements (13)
#define NEntry  9      
#define NMetal  5 
#define NElements  17 

double Metal[NMetal];
double Mass[NMetal][NEntry];
double Energy[NMetal][NEntry];
double Mrem[NMetal][NEntry];
double Elements[NMetal][NEntry][NElements];
char Tags[NElements][MaxCharactersInLine];

static void ReadNormal(FILE *fp){

    for(int l=0;l<NMetal;l++){
        fscanf(fp,"%*s %le",&(Metal[l]));
        fprintf(stderr,"Z = %g\n",Metal[l]);
 

        if(l==0){
            fscanf(fp,"%*s");
            for(int i=0;i<NElements;i++){
                fscanf(fp,"%s",Tags[i]);
            }
            fscanf(fp,"%*s");
        } else {
            fscanf(fp,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
        }

        for(int j=0;j<NEntry;j++){
            fscanf(fp,"%le",Mass[l]+j);

            for(int i=0;i<NElements;i++){
                fscanf(fp,"%le",Elements[l][j]+i);
            }
            fscanf(fp,"%le",Mrem[l]+j);
        }
    }

    return;
}

static void FillEmptylines(void){

    for(int l=0;l<NMetal;l++){
        for(int j=0;j<NEntry;j++){
            if(Mrem[l][j] == 0){
                fprintf(stderr,"Interpolation start! %d %d | %g < %g < %g\n",l,j,Mass[l][j-1],Mass[l][j],Mass[l][j+1]);

                for(int i=0;i<NElements;i++){
                    double grad = (Elements[l][j+1][i]-Elements[l][j-1][i])/(Mass[l][j+1]-Mass[l][j-1]);
                    Elements[l][j][i] = grad*(Mass[l][j]-Mass[l][j-1])+Elements[l][j-1][i];
                    fprintf(stderr,"-- %g %g %g\n",Elements[l][j-1][i],Elements[l][j][i],Elements[l][j+1][i]);
                }
                double grad = (Mrem[l][j+1]-Mrem[l][j-1])/(Mass[l][j+1]-Mass[l][j-1]);
                Mrem[l][j] = grad*(Mass[l][j]-Mass[l][j-1])+Mrem[l][j-1];
            }
        }
    }

    return ;
}

static void WriteData(FILE *fp){
    
    for(int i=0;i<NMetal;i++){  
        fprintf(fp,"%g\n",Metal[i]);
        fprintf(stderr,"%g\n",Metal[i]);
        for(int k=0;k<NEntry;k++){
            fprintf(fp,"%g ",Mass[i][k]);
            fprintf(fp,"%g ",Mrem[i][k]);
            for(int j=0;j<NElements;j++){
                fprintf(fp,"%g ",Elements[i][k][j]);
            }
            fprintf(fp,"\n");
        }
    }

    return ;
}

static char CELibSNIILabel[][MaxCharactersInLine] = {
    "CELibYield_H","CELibYield_He","CELibYield_C","CELibYield_N","CELibYield_O","CELibYield_Ne",
    "CELibYield_Mg","CELibYield_Si","CELibYield_S","CELibYield_Ca","CELibYield_Fe","CELibYield_Ni","CELibYield_Eu"};

static void WriteStructure(FILE *fp){

    fprintf(fp,"struct StructCELibSNIIYields CELibSNIIYieldsP98[CELibSNIIYields_Metallicity][CELibSNIIYields_Mass] = \{\n");
    fprintf(fp,"\t\{\n");

    
    for(int i=0;i<NMetal;i++){  
        for(int k=0;k<NEntry;k++){
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",Metal[i]);
            fprintf(fp,"\t.Mass = %g,\n",Mass[i][k]);
            fprintf(fp,"\t.Mr = %g,\n",Mrem[i][k]);
            fprintf(fp,"\t.Msw = %g,\n",0.0);
            fprintf(fp,"\t.Energy = %g,\n",1.e+51);

            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[0],Elements[i][k][0]);  //H
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[1],Elements[i][k][1]+Elements[i][k][2]); //He
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[2],Elements[i][k][3]+Elements[i][k][4]); //C
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[3],Elements[i][k][5]+Elements[i][k][6]); //N
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[4],Elements[i][k][7]+Elements[i][k][8]+Elements[i][k][9]); //O
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[5],Elements[i][k][10]+Elements[i][k][11]); //Ne
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[6],Elements[i][k][12]); //Mg
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[7],Elements[i][k][13]); //Si
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[8],Elements[i][k][14]); //S
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[9],Elements[i][k][15]); //Ca
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[10],Elements[i][k][16]); //Fe
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[11],0.e0); //Ni
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[12],0.e0); //Ni

            fprintf(fp,"\t},\n");
        }

        //fprintf(fp,"\t};,\n");
        fprintf(fp,"\t},\n");
        if(i!=NMetal-1){
            fprintf(fp,"\t\{\n");
        } else {
            fprintf(fp,"};\n");
        }
    }

    return ;
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
#elif 0
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
 * Anders & Grevesse 1989.
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
 * This function calculates the mass fractions of the material with the solar
 * abundances.
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
    Elements[CELibYield_Eu] = Scale*Mass*ElementNormalizedMasses[CELibYield_Eu];

    double Mtotal = 0.e0;
    for(int i=0;i<CELibElements;i++)
        Mtotal += Elements[i];
    double InvMtotal = Mass/Mtotal;
    for(int i=0;i<CELibElements;i++)
        Elements[i] *= InvMtotal;

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

    return ;
}

double ElementsSolar[NMetal][NEntry][CELibElements];
static void ConvertYield(void){

    for(int i=0;i<NMetal;i++){
        double Z = Metal[i];
        for(int j=0;j<NEntry;j++){
            double CurrentMass = Mass[i][j];
            double Elements[CELibElements];
            CELibSetMetallicityWithSolarAbundancePattern(CurrentMass,Elements,Z);
            fprintf(stderr," Check Y = %g %g\n",Elements[CELibYield_He],Elements[CELibYield_He]/CurrentMass);
            for(int k=0;k<CELibElements;k++){
                ElementsSolar[i][j][k] = Elements[k];
            }
        }
    }

    return ;
}

static void WriteYieldStructure(FILE *fp){

    
    fprintf(fp,"#define CELibSNIIYields_Metallicity_P98 (%d)\n",NMetal);
    fprintf(fp,"#define CELibSNIIYields_Mass_P98 (%d)\n",NEntry);

    fprintf(fp,"static struct StructCELibSNIIYields CELibSNIIYieldsP98[CELibSNIIYields_Metallicity_P98][CELibSNIIYields_Mass_P98] = \{\n");
    fprintf(fp,"\t\{\n");

    
    for(int i=0;i<NMetal;i++){  
        for(int k=0;k<NEntry;k++){
            double EjectaMassFraction = (Mass[i][k]-Mrem[i][k])/Mass[i][k];
            fprintf(fp,"\t\{\n");
            fprintf(fp,"\t.Metallicity = %g,\n",Metal[i]);
            fprintf(fp,"\t.Mass = %g,\n",Mass[i][k]);
            fprintf(fp,"\t.Mr = %g,\n",Mrem[i][k]);
            fprintf(fp,"\t.Msn = %g,\n",Mass[i][k]-Mrem[i][k]);
            fprintf(fp,"\t.Msw = %g,\n",0.0);
            fprintf(fp,"\t.Energy = %g,\n",1.e+51);

            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[0],Elements[i][k][0]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_H]);  //H
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[1],Elements[i][k][1]+Elements[i][k][2]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_He]); //He
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[2],Elements[i][k][3]+Elements[i][k][4]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_C]); //C
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[3],Elements[i][k][5]+Elements[i][k][6]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_N]); //N
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[4],Elements[i][k][7]+Elements[i][k][8]+Elements[i][k][9]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_O]); //O
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[5],Elements[i][k][10]+Elements[i][k][11]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_Ne]); //Ne
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[6],Elements[i][k][12]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_Mg]); //Mg
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[7],Elements[i][k][13]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_Si]); //Si
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[8],Elements[i][k][14]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_S]); //S
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[9],Elements[i][k][15]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_Ca]); //Ca
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[10],Elements[i][k][16]-EjectaMassFraction*ElementsSolar[i][k][CELibYield_Fe]); //Fe
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[11],0.e0); //Ni
            fprintf(fp,"\t.Elements[%s] = %g,\n",CELibSNIILabel[12],0.e0); //Eu

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

static void WriteYield(void){

    for(int i=0;i<NMetal;i++){  

        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"Yields.%02d",i);
        FileOpen(fp,fname,"w");


        for(int j=0;j<NEntry;j++){
            double EjectaMassFraction = (Mass[i][j]-Mrem[i][j])/Mass[i][j];
            fprintf(fp,"%g %g %g %g %g %g %g\n",Mass[i][j],Mrem[i][j]/Mass[i][j],
                    (Elements[i][j][0]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_H])/Mass[i][j],
                    (Elements[i][j][1]+Elements[i][j][2]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_He])/Mass[i][j],
                    (Elements[i][j][3]+Elements[i][j][4]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_C])/Mass[i][j],
                    (Elements[i][j][7]+Elements[i][j][8]+Elements[i][j][9]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_O])/Mass[i][j],
                    ((Elements[i][j][13]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_Si])+
                     (Elements[i][j][14]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_S])+
                     (Elements[i][j][15]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_Ca])+
                     (Elements[i][j][16]-EjectaMassFraction*ElementsSolar[i][j][CELibYield_Fe]))/Mass[i][j]);
        }
        fclose(fp);
    }

    return ;
}


int main(int argc, char **argv){

    FILE *fp;
    // read table 10 of Portinari et al. 1998.
    FileOpen(fp,"./table_10.txt","r");
    ReadNormal(fp);
    fclose(fp);

    FillEmptylines();
    ConvertYield();
    WriteYield();

    // Write data
    FileOpen(fp,"./Portinari+1998.dat","w");
    WriteData(fp);
    fclose(fp);

    // Write structure
    FileOpen(fp,"./PrePortinari+1998_Struct.h","w");
    WriteStructure(fp);
    fclose(fp);

    // Write structure
    FileOpen(fp,"./Portinari+1998_Struct.h","w");
    WriteYieldStructure(fp);
    fclose(fp);

    return EXIT_SUCCESS;
}
