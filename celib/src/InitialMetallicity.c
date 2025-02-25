#include "config.h"
#include "SNIIYields.h"

/*!
 * Cosmic helium abundance evaluated in Planck paper XVI(Planck Collaboration,
 * arXiv:1303.5076). This value is found in Table 2 of their paper.
 */
#define CELibCosmicHeliumAbundance  (0.247695)

/*!
 * The solar metallicity based on table 4 of Asplund et al. ARAA, 2009.
 * Note that this is the present day value.
 */
#define CELibSolarMetallicity   (0.0134) 

/*!
 * These element abundances are given in the Table 1 of Asplund et al. ARAA,
 * 2009.  They are defined as the logarithmic values and the H abundance is
 * defined to be log e_H = 12.00. Thus, the number denisty of the other element
 * X is log e_X = log(N_X/N_H) + 12.00, where N_X and N_H are the number density
 * of elements X and H. Since the Eu data is not provided in their paper, the
 * value given by Anders & Grevesse (1989) table 2 is adopted.
 */
static double MassFractionA09[3] = {0.7381,0.2485,0.0134};
static double ElementAbundancesA09[] = {
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
     0.51, // Eu // Adopted from Anders & Grevesse 1989 table 2.
};

/*!
 * Same as ElementAbundancesA09, but obtained from Grevesse & Sauval (1998).
 */
static double MassFractionGS98[3] = {0.735,0.248,0.017};
static double ElementAbundancesGS98[] = {
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

/*!
 * Same as ElementAbundancesA09, but obtained from Anders & Grevesse (1989).
 */
static double MassFractionAG89[3] = {0.70683,0.27431,0.01886};
static double ElementAbundancesAG89[] = {
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


/*!
 * Mass fractions X, Y,and Z. These values are evaluated by the adopted solar
 * metallicity table.
 */
static double MassFractionFromTable[3];

/*!
 * The mass of elements obtained from table 4 in Asplund et al. ARAA, 2009. For
 * simplicity, we adopt the value of the most frequent isotopes.
 */
static double ElementZ[CELibYield_Number] = {
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

static double ElementNumberDensities[CELibYield_Number];
static double ElementMasses[CELibYield_Number];
static double ElementNormalizedMasses[CELibYield_Number];

static bool first = true; //!< First call flag.


/*!
 * In this function, the mass fractions of the material with the solar
 * abundances is calculated. There are three solar abundance patterns and one
 * should select one of them. The default one is that of Asplund et al. (2009).
 * if CELibRunParameters.SolarAbundancePatternType == 0, the solar abundance pattern of Asplund et al. (2009) is used.
 * if CELibRunParameters.SolarAbundancePatternType == 1, the solar abundance pattern of Grevesse & Sauval (2009) is used.
 * if CELibRunParameters.SolarAbundancePatternType == 2, the solar abundance pattern of Anders & Grevesse (1989) is used.
 */
void CELibInitSolarAbundances(void){

    if(first == false) 
        return;

    double Mass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        if(CELibRunParameters.SolarAbundancePatternType == CELibSolarAbundancePattern_A09){
            ElementNumberDensities[i] = pow(10.0,ElementAbundancesA09[i]-12);
            ElementMasses[i] = ElementZ[i]*ElementNumberDensities[i];
        } else if(CELibRunParameters.SolarAbundancePatternType == CELibSolarAbundancePattern_GS98){
            ElementNumberDensities[i] = pow(10.0,ElementAbundancesGS98[i]-12);
            ElementMasses[i] = ElementZ[i]*ElementNumberDensities[i];
        } else if(CELibRunParameters.SolarAbundancePatternType == CELibSolarAbundancePattern_AG89){
            ElementNumberDensities[i] = pow(10.0,ElementAbundancesAG89[i]-12);
            ElementMasses[i] = ElementZ[i]*ElementNumberDensities[i];
        } else {
            fprintf(stderr,"Incorrect CELibRunParameters.SolarAbundancePatternType");
        }
    
        Mass += ElementMasses[i];
    }


    MassFractionFromTable[2] = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        ElementNormalizedMasses[i] = ElementMasses[i]/Mass;
        if(i>1)
            MassFractionFromTable[2] += ElementNormalizedMasses[i];
    }
    MassFractionFromTable[0] = ElementNormalizedMasses[0];
    MassFractionFromTable[1] = ElementNormalizedMasses[1];

    first = false;

    return ;
}

/*!
 * This function returns the mass of each element based on the primodrial
 * abandance pattern. Here, X_i is the mass fraction of an element i, and the
 * following condition is satisfied: 1 = \sum_i X_i. The return value of
 * Elements[i] is X_i * Mass.
 */
void CELibSetPrimordialMetallicity(const double Mass, double Elements[restrict]){ 
    
    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] = 0.e0;
    }
    Elements[CELibYield_H] = (1.0-CELibCosmicHeliumAbundance)*Mass;
    Elements[CELibYield_He] = CELibCosmicHeliumAbundance*Mass;

    return ;
}


/*!
 * This function returns the mass of each element with the adopted solar
 * abudance pattern. The total metallicity is set to "Metallicity". 
 */
void CELibSetMetallicityWithSolarAbundancePattern(const double Mass, double Elements[restrict], const double Metallicity){ 

    if(first == true)
        CELibInitSolarAbundances();

    double Factor_X = 1.0-MassFractionFromTable[1]-Metallicity;
    double Factor_Z = Metallicity/MassFractionFromTable[2];

    Elements[CELibYield_H]  = Mass*Factor_X;
    Elements[CELibYield_He] = Mass*ElementNormalizedMasses[CELibYield_He];
    Elements[CELibYield_C]  = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_C];
    Elements[CELibYield_N]  = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_N];
    Elements[CELibYield_O]  = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_O];
    Elements[CELibYield_Ne] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Ne];
    Elements[CELibYield_Mg] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Mg];
    Elements[CELibYield_Si] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Si];
    Elements[CELibYield_S]  = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_S];
    Elements[CELibYield_Ca] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Ca];
    Elements[CELibYield_Fe] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Fe];
    Elements[CELibYield_Ni] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Ni];
    Elements[CELibYield_Eu] = Factor_Z*Mass*ElementNormalizedMasses[CELibYield_Eu];

    return ;
}

/*!
 * This function returns the mass of each element based on the adopted solar
 * abundance pattern.
 */
void CELibSetSolarMetallicity(const double Mass, double Elements[restrict]){ 

    CELibSetMetallicityWithSolarAbundancePattern(Mass,Elements,CELibSolarMetallicity);

    return ;
}


/*!
 * This function just returns the metallicity ``Z'' of the solar abundance
 * pattern which is defined in CElib. The value of Asplund et al. (2009) is
 * adopted. Hence, this value is inconsistent when one adopts different solar
 * abundance patterns such as those of Grevesse & Sauval (2009) and Anders &
 * Grevesse (1989).
 */
double CELibGetMetalFractionForSolarChemicalComposision(void){
    return CELibSolarMetallicity;
}

/*!
 * This function returns the weight number of the given element.
 */
double CELibGetElementWeight(const int ElementIndex){
    return ElementZ[ElementIndex];
}

/*!
 * This function returns the number density of the solar abundance of the given
 * element. Values are evaluated under the adopted solar abundance pattern.
 */
double CELibGetElementNumberDensitySolar(const int ElementIndex){
    if(CELibRunParameters.SolarAbundancePatternType == CELibSolarAbundancePattern_A09){
        return pow(10.0,ElementAbundancesA09[ElementIndex]-12.0);
    } else if(CELibRunParameters.SolarAbundancePatternType == CELibSolarAbundancePattern_GS98){
        return pow(10.0,ElementAbundancesGS98[ElementIndex]-12.0);
    } else if(CELibRunParameters.SolarAbundancePatternType == CELibSolarAbundancePattern_AG89){
        return pow(10.0,ElementAbundancesAG89[ElementIndex]-12.0);
    } else {
        fprintf(stderr,"\t Incorrect CELibRunParameters.SolarAbundancePatternType is used.");
        return 0.0;
    }
}


/*!
 * This function returns the current adopted solar abundance pattern type ID.
 * If Asplund et al. (2009) is adopted, this function returns 0.
 * If Grevesse & Sauval (2009) is adopted, this function returns 1.
 * If Anders & Grevesse (1989) is adopted, this function returns 2.
 */
int CELibGetSolarAbundancePatternType(void){
    return CELibRunParameters.SolarAbundancePatternType;
}

/*!
 * This function sets the solar abundance pattern type ID.
 */
void CELibSetSolarAbundancePatternType(const int ID){
    CELibRunParameters.SolarAbundancePatternType = ID;
}
