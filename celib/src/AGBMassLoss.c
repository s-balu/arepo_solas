#include "config.h"
#include "AGBMassLoss.h"
#include "../data/Karakas2010/Karakas2010_Struct.h"
#include "../data/Doherty+2014/Doherty+2014_Struct.h"
#include "../data/CampbellLattanzio2008/CampbellLattanzio2008_Struct.h"
#include "../data/Gil-Pons+2013/Gil-Pons+2013_Struct.h"
#include "../data/vandenHoekGroenewegen1997/vandenHoekGroenewegen1997_Struct.h"

/*! \file AGBMassLoss.c
 * \brief This file has functions regarding the AGB mass loss.
 */

/*
 * Yields for AGB stars.
 * The yeilds provided by Karakas (2010) and Doherty et al. (2014).  This table includes the
 * amount of 11 elements, H, He C, N, O, Ne, Mg, Si, S, Fe, Ni, in Solar mass
 * for 4 metallicities from 0.001 to 0.02 and 16 masses from 1 Msun to 6 (6.5)
 * Msun.  When the metallicity of an SSP particle is less than 0.001, this
 * library adopts the yield for the metallicity of 0.001.
 * This table also includes the remnant mass for each star after its evolution.
 * Thus, the total ejecta mass can be caluculated easily.  The data is
 * rearranged by a ruby program.
 */


static int CELibAGBYields_Metallicity;
static int CELibAGBYields_Mass;

static int pCELibAGBYieldsCurrentIMFType;
static int pCELibAGBYieldsCurrentMetallicity;
static int pCELibAGBYieldsCurrentElement;

/* This is the quick look up table of metallicties. */
static double *CELibAGBYieldsMetallicityLookUpTable; 

struct CELibStructAGBYields *pCELibAGBYields; 
struct CELibStructAGBYields *pCELibAGBYieldsCombined; 
#define CELibAGBYields(x,y) (pCELibAGBYields[CELibAGBYields_Mass*(x)+(y)])

/*!
 * Make a combined yields table (Karakas 2010 and Deherty+2014)
 */
static void pCELibMakeAGBYieldsCombinedYieldsTable(void){

    pCELibAGBYieldsCombined = realloc(pCELibAGBYieldsCombined,
            sizeof(struct CELibStructAGBYields)*CELibAGBYields_Metallicity*CELibAGBYields_Mass);

    int counter = 0;
    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        for(int k=0;k<CELibAGBYields_Mass;k++){
            if(k<CELibAGBYields_Mass_K10){
                pCELibAGBYieldsCombined[counter] = CELibAGBYieldsK10[i][k];
            } else {
                if(i==0){
                    pCELibAGBYieldsCombined[counter] = CELibAGBYieldsD14[0][k-CELibAGBYields_Mass_K10];
                    pCELibAGBYieldsCombined[counter].Metallicity = CELibAGBYieldsK10[0][0].Metallicity;
                } else {
                    pCELibAGBYieldsCombined[counter] = CELibAGBYieldsD14[i-1][k-CELibAGBYields_Mass_K10];
                }
            }
            counter ++;
        }
    }

    return ;
}

/*!
 * This function returns interpolated yields based on the Campbell and Lattanzioi
 * (2008)'s yields table.
 */
static struct CELibStructAGBYields pCELibGetAGBYieldsIntepolatedValuesFromCL08Table(const double Mass){

    struct CELibStructAGBYields tmpAGB;

    if((Mass < CELibAGBYieldsCL08[0][0].Mass)||(Mass > CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Mass)){
        tmpAGB.Metallicity = 0.e0;
        tmpAGB.Mass = Mass;
        tmpAGB.Mr = Mass;
        tmpAGB.Mej = 0.e0;
        tmpAGB.Erf = 0.e0;
        tmpAGB.Energy = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            tmpAGB.Elements[i] = 0.e0;
        }
        return tmpAGB;
    }

    int IndexMass = (int)Mass-1;
    assert(IndexMass>=0);
    assert(IndexMass<=2);

    tmpAGB.Metallicity = 0.e0;
    tmpAGB.Mass = Mass;

    double dMass = Mass-CELibAGBYieldsCL08[0][IndexMass].Mass;
    double MassInterval = CELibAGBYieldsCL08[0][IndexMass+1].Mass-CELibAGBYieldsCL08[0][IndexMass].Mass;

    tmpAGB.Mr = (CELibAGBYieldsCL08[0][IndexMass+1].Mr-CELibAGBYieldsCL08[0][IndexMass].Mr)*dMass/MassInterval + CELibAGBYieldsCL08[0][IndexMass].Mr;
    //tmpAGB.Mej = (CELibAGBYieldsCL08[0][IndexMass+1].Mej-CELibAGBYieldsCL08[0][IndexMass].Mej)*dMass/MassInterval + CELibAGBYieldsCL08[0][IndexMass].Mej;
    tmpAGB.Mej = Mass-tmpAGB.Mr;
    tmpAGB.Erf = (CELibAGBYieldsCL08[0][IndexMass+1].Erf-CELibAGBYieldsCL08[0][IndexMass].Erf)*dMass/MassInterval + CELibAGBYieldsCL08[0][IndexMass].Erf;
    tmpAGB.Energy = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        tmpAGB.Elements[i] = (CELibAGBYieldsCL08[0][IndexMass+1].Elements[i]-CELibAGBYieldsCL08[0][IndexMass].Elements[i])*dMass/MassInterval
            + CELibAGBYieldsCL08[0][IndexMass].Elements[i];
    }

    return tmpAGB;
}


/*!
 * This function returns interpolated yields based on the Gil-Pons et al.
 * (2013)'s yields table.
 */
static struct CELibStructAGBYields pCELibGetAGBYieldsIntepolatedValuesFromG13Table(const double Mass){

    struct CELibStructAGBYields tmpAGB;
    if((CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Mass < Mass)&&(Mass < CELibAGBYieldsG13[0][0].Mass)){
        tmpAGB.Metallicity = 0.e0;
        tmpAGB.Mass = Mass;

        double MassMin = CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Mass;
        double MassMax = CELibAGBYieldsG13[0][0].Mass;
        
        double dMass = Mass-MassMin;
        double MassInterval = MassMax-MassMin;

        tmpAGB.Mr = (CELibAGBYieldsG13[0][0].Mr-CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Mr)*dMass/MassInterval 
                    + CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Mr;
        tmpAGB.Mej = Mass-tmpAGB.Mr;
        tmpAGB.Erf = (CELibAGBYieldsG13[0][0].Erf-CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Erf)*dMass/MassInterval 
                    + CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Erf;
        tmpAGB.Energy = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            tmpAGB.Elements[i] = (CELibAGBYieldsG13[0][0].Elements[i]-CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Elements[i])*dMass/MassInterval
                                 + CELibAGBYieldsCL08[0][CELibAGBYields_Mass_CL08-1].Elements[i];
        }
        return tmpAGB;
    } else if((Mass < CELibAGBYieldsG13[0][0].Mass)||(Mass > CELibAGBYieldsG13[0][CELibAGBYields_Mass_G13-1].Mass)){
        tmpAGB.Metallicity = 0.e0;
        tmpAGB.Mass = Mass;
        tmpAGB.Mr = Mass;
        tmpAGB.Mej = 0.e0;
        tmpAGB.Erf = 0.e0;
        tmpAGB.Energy = 0.e0;
        for(int i=0;i<CELibYield_Number;i++){
            tmpAGB.Elements[i] = 0.e0;
        }
        return tmpAGB;
    }

    int IndexMass = (int)Mass;
    assert(IndexMass>=4);
    assert(IndexMass<=9);
    IndexMass -= 4;

    tmpAGB.Metallicity = 0.e0;
    tmpAGB.Mass = Mass;

    double dMass = Mass-CELibAGBYieldsG13[0][IndexMass].Mass;
    double MassInterval = CELibAGBYieldsG13[0][IndexMass+1].Mass-CELibAGBYieldsG13[0][IndexMass].Mass;

    tmpAGB.Mr = (CELibAGBYieldsG13[0][IndexMass+1].Mr-CELibAGBYieldsG13[0][IndexMass].Mr)*dMass/MassInterval + CELibAGBYieldsG13[0][IndexMass].Mr;
    tmpAGB.Mej = Mass-tmpAGB.Mr;
    tmpAGB.Erf = (CELibAGBYieldsG13[0][IndexMass+1].Erf-CELibAGBYieldsG13[0][IndexMass].Erf)*dMass/MassInterval + CELibAGBYieldsG13[0][IndexMass].Erf;
    tmpAGB.Energy = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        tmpAGB.Elements[i] = (CELibAGBYieldsG13[0][IndexMass+1].Elements[i]-CELibAGBYieldsG13[0][IndexMass].Elements[i])*dMass/MassInterval
            + CELibAGBYieldsG13[0][IndexMass].Elements[i];
    }

    return tmpAGB;
}

/*!
 * Make a combined yields table (Karakas 2010, Deherty+2014, Campbell and Lattanzio
 * 2008, and Gil-Pons+2013)
 */
static void pCELibMakeAGBYieldsCombinedYieldsTablePopIII(void){

    pCELibAGBYieldsCombined = realloc(pCELibAGBYieldsCombined,
            sizeof(struct CELibStructAGBYields)*CELibAGBYields_Metallicity*CELibAGBYields_Mass);
    memset(pCELibAGBYieldsCombined,0,
            sizeof(struct CELibStructAGBYields)*CELibAGBYields_Metallicity*CELibAGBYields_Mass);

    // for zero metal
    int counter = 0;
    for(int k=0;k<CELibAGBYields_Mass;k++){
        double Mass;
        if(k<CELibAGBYields_Mass_K10){
            Mass = CELibAGBYieldsK10[0][k].Mass;
        } else {
            Mass = CELibAGBYieldsD14[0][k-CELibAGBYields_Mass_K10].Mass;
        }
        
        if(CELibRunParameters.PopIIIAGBYieldsTableID == CELibAGBZ0YieldsTableID_CL08){
            pCELibAGBYieldsCombined[counter] = pCELibGetAGBYieldsIntepolatedValuesFromCL08Table(Mass);
        } else if(CELibRunParameters.PopIIIAGBYieldsTableID == CELibAGBZ0YieldsTableID_CL08G13){
            pCELibAGBYieldsCombined[counter] = pCELibGetAGBYieldsIntepolatedValuesFromCL08Table(Mass);
            if(pCELibAGBYieldsCombined[counter].Mej == 0.e0){
                pCELibAGBYieldsCombined[counter] = pCELibGetAGBYieldsIntepolatedValuesFromG13Table(Mass);
            }
        }
        pCELibAGBYieldsCombined[counter].Metallicity = CELibRunParameters.PopIIIMetallicity;
        counter ++;
    }

    for(int i=1;i<CELibAGBYields_Metallicity;i++){
        for(int k=0;k<CELibAGBYields_Mass;k++){
            if(k<CELibAGBYields_Mass_K10){
                pCELibAGBYieldsCombined[counter] = CELibAGBYieldsK10[i-1][k];
            } else {
                if(i==1){
                    pCELibAGBYieldsCombined[counter] = CELibAGBYieldsD14[0][k-CELibAGBYields_Mass_K10];
                    pCELibAGBYieldsCombined[counter].Metallicity = CELibAGBYieldsK10[0][0].Metallicity;
                } else {
                    pCELibAGBYieldsCombined[counter] = CELibAGBYieldsD14[i-2][k-CELibAGBYields_Mass_K10];
                }
            }
            counter ++;
        }
    }

    return ;
}

/*!
 * This function sets the lookup table of metallicity.
 */
static void pCELibGetAGBYieldsLookUpTables(void){

    if((CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_K10)
            ||(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_K10D14)
            ||(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_vdHG97)){
        CELibAGBYieldsMetallicityLookUpTable = realloc(CELibAGBYieldsMetallicityLookUpTable,sizeof(double)*CELibAGBYields_Metallicity);
    } 
    if(CELibRunParameters.TestMode){
        fprintf(stderr,"//\t Show AGB metallicity lookup table.\n");
    }
    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        CELibAGBYieldsMetallicityLookUpTable[i] = CELibAGBYields(i,0).Metallicity;
        if(CELibRunParameters.TestMode){
            fprintf(stderr,"//\t\t Z[%d] = %g\n",i,CELibAGBYieldsMetallicityLookUpTable[i]);
        }
    }
    return ;
}


struct CELibStructAGBYields *CELibAGBYieldsIntegrated;

/*!
 * Allocate CELibAGBYieldsIntegrated.
 */
static void pCELibAllocateAGBYieldsAGBYieldsIntegrated(void){
    CELibAGBYieldsIntegrated = realloc(CELibAGBYieldsIntegrated,sizeof(struct CELibStructAGBYields)*CELibAGBYields_Metallicity);
    return ;
}

static char CELibAGBYieldsElementName[CELibYield_Number][MaxCharactersInLine];

/*!
 * Set element names in CELibAGBYieldsElementName[].
 */
static void pCELibSetAGBYieldsElementName(void){

    strcpy(CELibAGBYieldsElementName[CELibYield_H],    "H");
    strcpy(CELibAGBYieldsElementName[CELibYield_He],   "He");
    strcpy(CELibAGBYieldsElementName[CELibYield_C],    "C");
    strcpy(CELibAGBYieldsElementName[CELibYield_N],    "N");
    strcpy(CELibAGBYieldsElementName[CELibYield_O],    "O");
    strcpy(CELibAGBYieldsElementName[CELibYield_Ne],   "Ne");
    strcpy(CELibAGBYieldsElementName[CELibYield_Mg],   "Mg");
    strcpy(CELibAGBYieldsElementName[CELibYield_Si],   "Si");
    strcpy(CELibAGBYieldsElementName[CELibYield_S],    "S");
    strcpy(CELibAGBYieldsElementName[CELibYield_Ca],   "Ca");
    strcpy(CELibAGBYieldsElementName[CELibYield_Fe],   "Fe");
    strcpy(CELibAGBYieldsElementName[CELibYield_Ni],   "Ni");
    strcpy(CELibAGBYieldsElementName[CELibYield_Eu],   "Eu");

    return ;
}

/*!
 * This function writes the AGB yields table. Initial and final masses,
 * metallicities, and the yield of each element which is released to the ISM as
 * a wind.
 */
static void pCELibDumpAGBYieldsTable(const char OutDir[]){

    MakeDir(OutDir);

    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/AGBTable.%02d",OutDir,i);
        FileOpen(fp,fname,"w");
        for(int k=0;k<CELibAGBYields_Mass;k++){
            fprintf(fp,"%g %g %g ",CELibAGBYields(i,k).Mass,
                    CELibAGBYields(i,k).Mr,CELibAGBYields(i,k).Metallicity);
            for(int l=0;l<CELibYield_Number;l++){
                fprintf(fp,"%g ",CELibAGBYields(i,k).Elements[l]);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }

    return ;
}


/*!
 *  This function retuns the interpolated yield for a given element at given
 *  mass and metallicity
 *  y = (M_i,u-M_i,l)/(M,u-M,l)*(M-M,l)+M_i,l
 */
static double pCELibGetAGBYieldsInterpolatedYieldMass(const int IndexMetal, const int TargetElement, const double Mass){

    int IndexMass = 0;
    if(CELibAGBYields(IndexMetal,0).Mass > Mass){
        return 0.e0;
    }
    if(CELibAGBYields(IndexMetal,CELibAGBYields_Mass-1).Mass < Mass){
        return 0.e0;
    }

    for(int i=1;i<CELibAGBYields_Mass;i++){
        if(CELibAGBYields(IndexMetal,i).Mass > Mass){
            IndexMass = i-1;
            break;
        }
    }

    return (CELibAGBYields(IndexMetal,IndexMass+1).Elements[TargetElement]-CELibAGBYields(IndexMetal,IndexMass).Elements[TargetElement])/
        (CELibAGBYields(IndexMetal,IndexMass+1).Mass-CELibAGBYields(IndexMetal,IndexMass).Mass)*
        (Mass-CELibAGBYields(IndexMetal,IndexMass).Mass)+CELibAGBYields(IndexMetal,IndexMass).Elements[TargetElement];
}


static double pCELibAGBYieldsGetYieldMassFunction(const double Mass){
    double ElementMass = pCELibGetAGBYieldsInterpolatedYieldMass(pCELibAGBYieldsCurrentMetallicity,pCELibAGBYieldsCurrentElement,Mass);
    return ElementMass*CELibIMF[pCELibAGBYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}


/*!
 * This function returns the IMF weighted released mass fraction of the element
 * pCELibAGBYieldsCurrentElement during the AGB mass loss phase. 
 *
 *  p(m_l < m < m_u) = \int_m_l^m_u y_i \xi(\log m) dm/m
 *
 * IMF type, metallicity, and element IDs are determined by
 * pCELibAGBYieldsCurrentIMFType, pCELibAGBYieldsCurrentMetallicity and
 * pCELibAGBYieldsCurrentElement, respectively.
 */
static double CELibGetAGBYieldsReturnElementMassFractionInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibAGBYieldsGetYieldMassFunction);
}


/*! 
 *  This function retuns the return mass fraction at a given mass.
 *  m-Mr_interpolated 
 *  where  Mr_interpolated = (Mr,u-Mr,l)/(M,u-M,l)*(M-M,l)+Mr,l.
 */
static double pCELibGetAGBYieldsInterpolatedReturnMass(const int IndexMetal, const double Mass_in_SolarMass){

    int IndexMass = 0;
    if(CELibAGBYields(IndexMetal,0).Mass > Mass_in_SolarMass){
        return 0.e0;
    }
    if(CELibAGBYields(IndexMetal,CELibAGBYields_Mass-1).Mass < Mass_in_SolarMass){
        return 0.e0;
    }

    for(int i=1;i<CELibAGBYields_Mass;i++){
        if(CELibAGBYields(IndexMetal,i).Mass > Mass_in_SolarMass){
            IndexMass = i-1;
            break;
        }
    }

    return Mass_in_SolarMass-((CELibAGBYields(IndexMetal,IndexMass+1).Mr-CELibAGBYields(IndexMetal,IndexMass).Mr)/
        (CELibAGBYields(IndexMetal,IndexMass+1).Mass-CELibAGBYields(IndexMetal,IndexMass).Mass)*
        (Mass_in_SolarMass-CELibAGBYields(IndexMetal,IndexMass).Mass)+CELibAGBYields(IndexMetal,IndexMass).Mr);
}


/*!
 * This function returns the interpolated return mass fraction.
 */
static double pCELibGetAGBYieldsReturnMassFractionFunction(const double Mass){
    double ReturnMass = pCELibGetAGBYieldsInterpolatedReturnMass(pCELibAGBYieldsCurrentMetallicity,Mass);
    return ReturnMass*CELibIMF[pCELibAGBYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 * This function returns the IMF weighted, total return mass fraction during the
 * AGB mass loss phase. 
 *
 *  Erf(m_l < m < m_u) = \int_m_l^m_u (m-Mr) \xi(\log m) dm/m
 *
 * IMF type, metallicity, and element IDs are determined by
 * pCELibAGBYieldsCurrentIMFType, pCELibAGBYieldsCurrentMetallicity and
 * pCELibAGBYieldsCurrentElement, respectively.
 */
static double CELibGetAGBYieldsTotalReturnMassFractionInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetAGBYieldsReturnMassFractionFunction);
}


/*!
 *
 */
double CELibGetAGBYieldsReturnMassInGivenTimeRange(const double Time_start_in_year, const double Time_end_in_year, const double Metallicity){

    double MassMin = CELibGetDyingStellarMass(Time_end_in_year,Metallicity);
    double MassMax = CELibGetDyingStellarMass(Time_start_in_year,Metallicity);
    if(MassMax > 100) MassMax = 100;
    if(MassMin < 0.9) MassMin = 0.1;
    // fprintf(stderr,"%g %g / %g\n",MassMax,MassMin,Time_start_in_year);

    int SaveIMFType = pCELibAGBYieldsCurrentIMFType;
    if((CELibRunParameters.PopIIIIMF == 1)&&(CELibRunParameters.SNIIYieldsTableID == CELibAGBYieldsTableID_K10D14)){
        if(Metallicity < CELibRunParameters.PopIIIMetallicity){
            pCELibAGBYieldsCurrentIMFType = CELibIMF_Susa;
        }
    }
    double Er = IntegralSimpson(MassMin,MassMax,CELibRunParameters.IntegrationSteps,&pCELibGetAGBYieldsReturnMassFractionFunction);
    pCELibAGBYieldsCurrentIMFType = SaveIMFType;

    return Er;
}

/*!
 *  This function retuns the remnant mass fraction at a given mass.
 */
static double pCELibGetAGBYieldsInterpolatedRemnantMass(const int IndexMetal, const double Mass_in_SolarMass){

    int IndexMass = 0;
    if(CELibAGBYields(IndexMetal,0).Mass > Mass_in_SolarMass){
        return 0.e0;
    }
    if(CELibAGBYields(IndexMetal,CELibAGBYields_Mass-1).Mass < Mass_in_SolarMass){
        return 0.e0;
    }

    for(int i=1;i<CELibAGBYields_Mass;i++){
        if(CELibAGBYields(IndexMetal,i).Mass > Mass_in_SolarMass){
            IndexMass = i-1;
            break;
        }
    }

    return (CELibAGBYields(IndexMetal,IndexMass+1).Mr-CELibAGBYields(IndexMetal,IndexMass).Mr)/
        (CELibAGBYields(IndexMetal,IndexMass+1).Mass-CELibAGBYields(IndexMetal,IndexMass).Mass)*
        (Mass_in_SolarMass-CELibAGBYields(IndexMetal,IndexMass).Mass)+CELibAGBYields(IndexMetal,IndexMass).Mr;
}


/*!
 * This function returns the remnant mass fraction at a given mass. 
 */
static double pCELibGetAGBYieldsRemnantMassFractionFunction(const double Mass){
    double ReturnMass = pCELibGetAGBYieldsInterpolatedRemnantMass(pCELibAGBYieldsCurrentMetallicity,Mass);
    return ReturnMass*CELibIMF[pCELibAGBYieldsCurrentIMFType].IMFFunctionPerMass(Mass);
}

/*!
 * This function returns the IMF weighted, total remnant mass fraction during the
 * AGB mass loss phase. 
 *
 *  Mr(m_l < m < m_u) = \int_m_l^m_u Mr(m) \xi(\log m) dm/m
 *
 * IMF type, metallicity, and element IDs are determined by
 * pCELibAGBYieldsCurrentIMFType, pCELibAGBYieldsCurrentMetallicity and
 * pCELibAGBYieldsCurrentElement, respectively.
 */
static double CELibGetAGBYieldsRemnantMassFractionInGivenMassRange(const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){
    return IntegralSimpson(LowerMass_in_SolarMass,UpperMass_in_SolarMass,CELibRunParameters.IntegrationSteps,&pCELibGetAGBYieldsRemnantMassFractionFunction);
}

/*!
 * This function makes the integrated return mass fraction table of all
 * elements. Each element is weighted by the selected IMF and is normalzied by
 * mass in the unit of the solar mass.  Therefore, p_i * SSP_mass provides us
 * the amount of the return fraction for the element i.  In order to get the
 * total amount of the ejecta, 
 *
 * Mej_i = p_i * SSP_mass + Z_i * Erf * SSP_mass,
 *
 *       where Erf is the total return mass fraction.
 * Note that 
 * CELibAGBYieldsIntegrated[].Elements[i] = p_i
 * and 
 * CELibAGBYieldsIntegrated[].Erf = Erf.
 */
static void pCELibMakeAGBYieldsIntegratedTable(void){

    pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;

    if(CELibRunParameters.TestMode)
        fprintf(stderr,"//\t CELib integrated AGB MassLoss Table.\n");

    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        if((CELibRunParameters.PopIIIIMF == 1)&&(CELibRunParameters.PopIIIAGB == 1)&&(i==0)){
            pCELibAGBYieldsCurrentIMFType = CELibIMF_Susa;
        } else {
            pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;
        }

        CELibAGBYieldsIntegrated[i].Metallicity = CELibAGBYields(i,0).Metallicity;
        pCELibAGBYieldsCurrentMetallicity = i;
        CELibAGBYieldsIntegrated[i].Erf = CELibGetAGBYieldsTotalReturnMassFractionInGivenMassRange(CELibRunParameters.AGBLowerMass,CELibRunParameters.AGBUpperMass);
        CELibAGBYieldsIntegrated[i].Mr = CELibGetAGBYieldsRemnantMassFractionInGivenMassRange(CELibRunParameters.AGBLowerMass,CELibRunParameters.AGBUpperMass);
        if(CELibRunParameters.TestMode){
            fprintf(stderr,"//\t AGB Yields\n");
            fprintf(stderr,"//\t\t Metallicity = %g\n",CELibAGBYieldsIntegrated[i].Metallicity);
            fprintf(stderr,"//\t\t Mr = %g\n",CELibAGBYieldsIntegrated[i].Mr);
            fprintf(stderr,"//\t\t Erf = %g\n",CELibAGBYieldsIntegrated[i].Erf);
        }

        if(CELibAGBYieldsIntegrated[i].Erf >= 1.0){
            dprint(pCELibAGBYieldsCurrentIMFType);
            gprint(CELibRunParameters.AGBLowerMass);
            gprint(CELibRunParameters.AGBUpperMass);
            gprint(CELibAGBYieldsIntegrated[i].Erf);
            assert(CELibAGBYieldsIntegrated[i].Erf < 1.0);
        }
        if(CELibAGBYieldsIntegrated[i].Erf < 0.0){
            gprint(CELibAGBYieldsIntegrated[i].Erf);
            assert(CELibAGBYieldsIntegrated[i].Erf >= 0.0);
        }

        assert(CELibAGBYieldsIntegrated[i].Erf < 1.0);
        for(int j=0;j<CELibYield_Number;j++){
            pCELibAGBYieldsCurrentElement = j;
            CELibAGBYieldsIntegrated[i].Elements[j] =
                CELibGetAGBYieldsReturnElementMassFractionInGivenMassRange(CELibRunParameters.AGBLowerMass,CELibRunParameters.AGBUpperMass);
        }
        if(CELibRunParameters.TestMode){
            for(int j=0;j<CELibYield_Number;j++){
                fprintf(stderr,"//\t\t %-2s %g\n",CELibAGBYieldsElementName[j],CELibAGBYieldsIntegrated[i].Elements[j]);
            }
            fprintf(stderr,"\n");
        }
    }

    return ;
}


/*!
 * This function writes the values of integrated yields, p_i, the total return
 * mass fraction and metalicity in files, CELib/CELibAGB/AGBYieldsTableIntegral.??.
 */
static void pCELibWriteAGBYieldsIntegratedTable(const char OutDir[]){

    MakeDir(OutDir);

    pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;

    // if(CELibRunParameters.TestMode)
        // fprintf(stderr,"CELib integrated AGB MassLoss Table.\n");

    for(int i=0;i<CELibAGBYields_Metallicity;i++){

        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"%s/AGBYieldsTableIntegral.%02d",OutDir,i);
        FileOpen(fp,fname,"w");

        double Sum = 0.e0;
        for(int j=0;j<CELibYield_Number;j++){
            Sum += CELibAGBYieldsIntegrated[i].Elements[j];
            if(CELibRunParameters.TestMode)
                fprintf(fp,"%g %s\n",CELibAGBYieldsIntegrated[i].Elements[j],CELibAGBYieldsElementName[j]);
        }
        if(CELibRunParameters.TestMode){
            fprintf(fp,"Mr = %g\n",Sum);
            fprintf(fp,"Erf = %g\n",CELibAGBYieldsIntegrated[i].Erf);
            fprintf(fp,"Metallicity = %g\n",CELibAGBYieldsIntegrated[i].Metallicity);
            fprintf(fp,"\n");
        }

        fclose(fp);
    }

    return ;
}


static int CELibStructAGBYieldsArraySize;
static struct CELibStructAGBYields **CELibAGBYieldsIntegratedBin = NULL;

static struct CELibStructAGBYieldsIntegratedBinInfo{
    double AgeUpper;
    double AgeLower;
    double MassUpper;
    double MassLower;
} **CELibAGBYieldsIntegratedBinInfo = NULL;
//};
//
//static struct CELibStructAGBYieldsIntegratedBinInfo **CELibAGBYieldsIntegratedBinInfo = NULL;

/*!
 * This function allocates the array "CELibAGBYieldsIntegratedBin[][]"
 */
static void pCELibAllocateAGBYieldsIntegratedBin(const int Number){

    CELibAGBYieldsIntegratedBin = malloc(sizeof(struct CELibStructAGBYields *)*CELibAGBYields_Metallicity);
    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        CELibAGBYieldsIntegratedBin[i] = malloc(sizeof(struct CELibStructAGBYields)*Number);
    }
    CELibStructAGBYieldsArraySize = CELibAGBYields_Metallicity;

    return ;
}


/*!
 * This function releases the array "CELibAGBYieldsIntegratedBin[][]".
 */
static void pCELibReleaseAGBYieldsIntegratedBin(void){

    for(int i=0;i<CELibStructAGBYieldsArraySize;i++){
        free(CELibAGBYieldsIntegratedBin[i]);
    }
    free(CELibAGBYieldsIntegratedBin);
    return ;
}


/*!
 * This function allocates the array "CELibAGBYieldsIntegratedBinInfo[][]"
 */
static void pCELibAllocateAGBYieldsIntegratedBinInfo(const int Number){

    CELibAGBYieldsIntegratedBinInfo = malloc(sizeof(struct CELibStructAGBYieldsIntegratedBinInfo *)*CELibAGBYields_Metallicity);
    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        CELibAGBYieldsIntegratedBinInfo[i] = malloc(sizeof(struct CELibStructAGBYieldsIntegratedBinInfo)*Number);
    }

    return ;
}


/*!
 * This function releases the array "CELibAGBYieldsIntegratedBinInfo[][]".
 */
static void pCELibReleaseAGBYieldsIntegratedBinInfo(void){

    for(int i=0;i<CELibStructAGBYieldsArraySize;i++){
        free(CELibAGBYieldsIntegratedBinInfo[i]);
    }
    free(CELibAGBYieldsIntegratedBinInfo);
    return ;
}


/*!
 * This function makes the integrated AGB yields bins. The bin size is
 * given by "RunParameters.AGBBinNumber". If "RunParameters.AGBBinType == 1",
 * the log based bin is used. 
 */
static void pCELibMakeAGBYieldsIntegratedBinTable(void){

    pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;

    // Set bin max and min.
    double BinMin = fmax(CELibRunParameters.AGBBinLowerAge,CELibGetLifeTimeofStar(CELibRunParameters.AGBUpperMass,CELibAGBYieldsMetallicityLookUpTable[0]));
    double BinMax = fmin(CELibRunParameters.AGBBinUpperAge,CELibGetLifeTimeofStar(CELibRunParameters.AGBLowerMass,CELibAGBYieldsMetallicityLookUpTable[0]));

    BinMin = CELibRunParameters.AGBBinLowerAge;
    BinMax = CELibRunParameters.AGBBinUpperAge;

    if(CELibRunParameters.AGBBinType == CELibAGBRateModelID_LinearBin){
        CELibRunParameters.AGBBinNumber = (int)(CELibRunParameters.AGBBinUpperAge/CELibRunParameters.AGBBinTimeInterval)+1;
    }
    if(CELibAGBYieldsIntegratedBin != NULL){
        pCELibReleaseAGBYieldsIntegratedBin();
    }
    if(CELibAGBYieldsIntegratedBinInfo != NULL){
        pCELibReleaseAGBYieldsIntegratedBinInfo();
    }

    pCELibAllocateAGBYieldsIntegratedBin(CELibRunParameters.AGBBinNumber);
    pCELibAllocateAGBYieldsIntegratedBinInfo(CELibRunParameters.AGBBinNumber);

    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        if((CELibRunParameters.PopIIIIMF == 1)&&(CELibRunParameters.PopIIIAGB == 1)&&(i==0)){
            pCELibAGBYieldsCurrentIMFType = CELibIMF_Susa;
        } else {
            pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;
        }

        pCELibAGBYieldsCurrentMetallicity = i;
        // Make bin and put date.
        if(CELibRunParameters.AGBBinType == CELibAGBRateModelID_LogBin){

            double BinMinLog = log10(BinMin);
            double BinMaxLog = log10(BinMax);

            double dBin = (BinMaxLog-BinMinLog)/CELibRunParameters.AGBBinNumber;
            for(int k=0;k<CELibRunParameters.AGBBinNumber;k++){
                double _BinMin = pow(10.0,BinMinLog+dBin*k);
                double _BinMax = pow(10.0,BinMinLog+dBin*(k+1));
                CELibAGBYieldsIntegratedBinInfo[i][k].AgeLower = _BinMin;
                CELibAGBYieldsIntegratedBinInfo[i][k].AgeUpper = _BinMax;

                CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper = CELibGetDyingStellarMass(_BinMin,CELibAGBYields(i,0).Metallicity);
                CELibAGBYieldsIntegratedBinInfo[i][k].MassLower = CELibGetDyingStellarMass(_BinMax,CELibAGBYields(i,0).Metallicity);

                if(CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper == __CELib_NoValue__){
                    CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper = CELibIMF[CELibRunParameters.IMFType].MassMax;
                } 
                if(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower == __CELib_NoValue__){
                    if(CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper >= CELibIMF[CELibRunParameters.IMFType].MassMax){
                        CELibAGBYieldsIntegratedBinInfo[i][k].MassLower = CELibIMF[CELibRunParameters.IMFType].MassMax;
                    } else {
                        CELibAGBYieldsIntegratedBinInfo[i][k].MassLower = CELibIMF[CELibRunParameters.IMFType].MassMin;
                    }
                } 

                CELibAGBYieldsIntegratedBin[i][k].Metallicity = CELibAGBYields(i,0).Metallicity;
                CELibAGBYieldsIntegratedBin[i][k].Erf = 
                    CELibGetAGBYieldsTotalReturnMassFractionInGivenMassRange(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower,
                            CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);
                /////////////////////////////////////////////////
                if(CELibAGBYieldsIntegratedBin[i][k].Erf >= 1.0){ // check range.
                    dprint(pCELibAGBYieldsCurrentIMFType);
                    gprint(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower);
                    gprint(CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);
                    gprint(CELibAGBYieldsIntegratedBin[i][k].Erf);
                    assert(CELibAGBYieldsIntegratedBin[i][k].Erf < 1.0);
                }
                if(CELibAGBYieldsIntegratedBin[i][k].Erf < 0.0){ // check range.
                    gprint(CELibAGBYieldsIntegratedBin[i][k].Erf);
                    assert(CELibAGBYieldsIntegratedBin[i][k].Erf >= 0.0);
                }
                /////////////////////////////////////////////////

                for(int l=0;l<CELibYield_Number;l++){
                    pCELibAGBYieldsCurrentElement = l;
                    CELibAGBYieldsIntegratedBin[i][k].Elements[l] = 
                        CELibGetAGBYieldsReturnElementMassFractionInGivenMassRange(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower,CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);
                }
            }
        } else {

            for(int k=0;k<CELibRunParameters.AGBBinNumber;k++){
                double _BinMin = CELibRunParameters.AGBBinTimeInterval*k;
                double _BinMax = CELibRunParameters.AGBBinTimeInterval*(k+1);
                CELibAGBYieldsIntegratedBinInfo[i][k].AgeLower = fmax(_BinMin,CELibRunParameters.AGBBinLowerAge);
                CELibAGBYieldsIntegratedBinInfo[i][k].AgeUpper = fmax(_BinMax,CELibRunParameters.AGBBinLowerAge);

                CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper = CELibGetDyingStellarMass(CELibAGBYieldsIntegratedBinInfo[i][k].AgeLower,CELibAGBYieldsMetallicityLookUpTable[i]);
                CELibAGBYieldsIntegratedBinInfo[i][k].MassLower = CELibGetDyingStellarMass(CELibAGBYieldsIntegratedBinInfo[i][k].AgeUpper,CELibAGBYieldsMetallicityLookUpTable[i]);

                if((CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper == __CELib_NoValue__)&&(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower == __CELib_NoValue__)){
                    CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper = CELibAGBYieldsIntegratedBinInfo[i][k].MassLower = CELibIMF[CELibRunParameters.IMFType].MassMax;
                } else if (CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper == __CELib_NoValue__){
                    CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper = CELibIMF[CELibRunParameters.IMFType].MassMax;
                } else if (CELibAGBYieldsIntegratedBinInfo[i][k].MassLower == __CELib_NoValue__){
                    CELibAGBYieldsIntegratedBinInfo[i][k].MassLower = CELibIMF[CELibRunParameters.IMFType].MassMin;
                }


                CELibAGBYieldsIntegratedBin[i][k].Metallicity = CELibAGBYields(i,0).Metallicity;
                CELibAGBYieldsIntegratedBin[i][k].Erf = 
                    CELibGetAGBYieldsTotalReturnMassFractionInGivenMassRange(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower,
                            CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);

                /////////////////////////////////////////////////
                if(CELibAGBYieldsIntegratedBin[i][k].Erf >= 1.0){ // check range.
                    dprint(pCELibAGBYieldsCurrentIMFType);
                    gprint(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower);
                    gprint(CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);
                    gprint(CELibAGBYieldsIntegratedBin[i][k].Erf);
                    assert(CELibAGBYieldsIntegratedBin[i][k].Erf < 1.0);
                }
                if(CELibAGBYieldsIntegratedBin[i][k].Erf < 0.0){ // check range.
                    gprint(CELibAGBYieldsIntegratedBin[i][k].Erf);
                    assert(CELibAGBYieldsIntegratedBin[i][k].Erf >= 0.0);
                }
                /////////////////////////////////////////////////

                for(int l=0;l<CELibYield_Number;l++){
                    pCELibAGBYieldsCurrentElement = l;
                    CELibAGBYieldsIntegratedBin[i][k].Elements[l] = 
                        CELibGetAGBYieldsReturnElementMassFractionInGivenMassRange(CELibAGBYieldsIntegratedBinInfo[i][k].MassLower,CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);
                }
            }
        }
    }
    return ;
}


struct CELibStructFeedbackOutput CELibGetAGBYieldsIntegratedInGivenMassRange(const int MetalID, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){

    struct CELibStructFeedbackOutput AGBFeedback;

    pCELibAGBYieldsCurrentMetallicity = MetalID;

    double InitElements[CELibYield_Number];
    if((CELibRunParameters.PopIIIIMF == 1)&&(CELibRunParameters.PopIIIAGB == 1)&&(MetalID == 0)){
        pCELibAGBYieldsCurrentIMFType = CELibIMF_Susa;
        CELibSetPrimordialMetallicity(1.0,InitElements); 
    } else {
        pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;
        CELibSetMetallicityWithSolarAbundancePattern(1.0,InitElements,CELibAGBYieldsMetallicityLookUpTable[MetalID]); 
    }

    double Erf = CELibGetAGBYieldsTotalReturnMassFractionInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    // fprintf(stderr,"//\t\t AGB: Er= %g\n",Erf);
    double Elements[CELibYield_Number];
    double EjectaMass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        pCELibAGBYieldsCurrentElement = i;
        Elements[i] = CELibGetAGBYieldsReturnElementMassFractionInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass)
                +InitElements[i]*Erf;
        EjectaMass += Elements[i];
        char Name[MaxCharactersInLine];
        CELibGetSNIIYieldElementName(i,Name);
        //fprintf(stderr,"//\t\t AGB: Elements %s = %g\n",Name,Elements[i]);
        /// fprintf(stderr,"//\t\t AGB: Elements %s = %g %g %g\n",Name,Elements[i],CELibGetAGBYieldsReturnElementMassFractionInGivenMassRange(LowerMass_in_SolarMass,UpperMass_in_SolarMass),InitElements[i]*Erf);
    }
    // fprintf(stderr,"//\t\t AGB: Ejecta mass = %g\n",EjectaMass);

    AGBFeedback.EjectaMass = 0.0;
    for(int i=0;i<CELibYield_Number;i++){
        AGBFeedback.Elements[i] = Elements[i];
        AGBFeedback.EjectaMass += Elements[i];
    }

    return AGBFeedback;
}

struct CELibStructFeedbackOutput CELibGetAGBYieldsIntegratedInGivenMassRangeZ(const double Metallicity, const double LowerMass_in_SolarMass, const double UpperMass_in_SolarMass){

    struct CELibStructFeedbackOutput AGBFeedback;


    if(Metallicity <= CELibRunParameters.PopIIIMetallicity){
        AGBFeedback = CELibGetAGBYieldsIntegratedInGivenMassRange(0,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    } else if(Metallicity >= CELibAGBYields(CELibAGBYields_Metallicity-1,0).Metallicity){
        AGBFeedback = CELibGetAGBYieldsIntegratedInGivenMassRange(CELibAGBYields_Metallicity-1,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
    } else { // interpolation

        int IndexZ = CELibAGBYields_Metallicity-2;
        for(int i=1;i<CELibAGBYields_Metallicity;i++){
            if(CELibAGBYields(i,0).Metallicity > Metallicity){
                IndexZ = i-1;
                break;
            }
        }

        struct CELibStructFeedbackOutput AGBFeedback1 = 
            CELibGetAGBYieldsIntegratedInGivenMassRange(IndexZ,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
        struct CELibStructFeedbackOutput AGBFeedback2 = 
            CELibGetAGBYieldsIntegratedInGivenMassRange(IndexZ+1,LowerMass_in_SolarMass,UpperMass_in_SolarMass);
#if 0
        AGBFeedback.Energy = (AGBFeedback2.Energy-AGBFeedback1.Energy)/(CELibAGBYields(IndexZ+1,0).Metallicity-CELibAGBYields(IndexZ,0).Metallicity)
                        *(Metallicity-CELibAGBYields(IndexZ,0).Metallicity)+AGBFeedback1.Energy;
#endif
        double EjectaMass = 0.e0;
        for(int k=0;k<CELibYield_Number;k++){
            AGBFeedback.Elements[k] = 
                (AGBFeedback2.Elements[k]-AGBFeedback1.Elements[k])/(CELibAGBYields(IndexZ+1,0).Metallicity-CELibAGBYields(IndexZ,0).Metallicity)
                                *(Metallicity-CELibAGBYields(IndexZ,0).Metallicity)+AGBFeedback1.Elements[k];
            EjectaMass += AGBFeedback.Elements[k];
            // fprintf(stderr,"%g %g | %g %g \n",AGBFeedback2.Elements[k],
                    // AGBFeedback1.Elements[k],CELibAGBYields(IndexZ+1,0).Metallicity,
                    // CELibAGBYields(IndexZ,0).Metallicity);
        }
        AGBFeedback.EjectaMass = EjectaMass;
        // gprint(AGBFeedback1.EjectaMass);
        // gprint(AGBFeedback2.EjectaMass);
        // gprint(AGBFeedback.EjectaMass);
        // for(int i=0;i<CELibAGBYields_Metallicity;i++){
            // gprint(CELibAGBYields(i,0).Metallicity);
            // dprint(IndexZ);
            // gprint(CELibAGBYields(IndexZ,0).Metallicity);
            // gprint(CELibAGBYields(IndexZ+1,0).Metallicity);
        // }
    }

    return AGBFeedback;
}

/*!
 * This function writes the values of integrated AGB mass loss bins in files
 * "CELib/CELibAGB/AGBIntegratedMassLossBin.??.??". 
 */
static void pCELibWriteAGBYieldsIntegratedBinTable(const char OutDir[]){

    MakeDir(OutDir);

    for(int i=0;i<CELibAGBYields_Metallicity;i++){
        for(int k=0;k<CELibRunParameters.AGBBinNumber;k++){
            FILE *fp;
            char fname[MaxCharactersInLine];
            Snprintf(fname,"%s/AGBYieldsIntegratedMassBin.%02d.%02d",OutDir,i,k);
            FileOpen(fp,fname,"w");

            fprintf(fp,"Bin number = %d\n",k);
            fprintf(fp,"Bin Mass Range = %g -- %g\n",CELibAGBYieldsIntegratedBinInfo[i][k].MassLower,CELibAGBYieldsIntegratedBinInfo[i][k].MassUpper);
            fprintf(fp,"Bin Time Range = %g -- %g\n",CELibAGBYieldsIntegratedBinInfo[i][k].AgeLower,CELibAGBYieldsIntegratedBinInfo[i][k].AgeUpper);
            fprintf(fp,"Metallicity = %g\n",CELibAGBYieldsIntegratedBin[i][k].Metallicity);
            fprintf(fp,"Return Mass fraction = %g\n",CELibAGBYieldsIntegratedBin[i][k].Erf);
            for(int l=0;l<CELibYield_Number;l++){
                fprintf(fp,"%g %s\n",CELibAGBYieldsIntegratedBin[i][k].Elements[l],CELibAGBYieldsElementName[l]);
            }

            fclose(fp);
        }
    }

    return ;
}


/*!
 * This function write the return mass fractions from AGBs in a file.
 */
static void pCELibWriteAGBYieldsReturnMassFraction(const char OutDir[]){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"%s/AGBRetrunMassFraction.dat",OutDir);
    FileOpen(fp,fname,"w");
    for(int k=0;k<CELibIMF_NTypes;k++){
        pCELibAGBYieldsCurrentIMFType = k;
        fprintf(fp,"IMF Type: %d\n",k);
        for(int i=0;i<CELibAGBYields_Metallicity;i++){
            double Z =  CELibAGBYields(i,0).Metallicity;
            pCELibAGBYieldsCurrentMetallicity = i;
            double Erf = CELibGetAGBYieldsTotalReturnMassFractionInGivenMassRange(CELibRunParameters.AGBLowerMass,CELibRunParameters.AGBUpperMass);
            fprintf(fp,"Metallicity = %g\n",Z);
            fprintf(fp,"Erf = %g\n",Erf);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);


    pCELibAGBYieldsCurrentIMFType = CELibRunParameters.IMFType;

    return ;
}

/*!
 * The adopted AGB feedback model is tested with various cases in this function.
 */
static void pCELibCheckAGBYieldsFeedBack(const char OutDir[], const int Type){

    MakeDir(OutDir);

    FILE *fp;
    char fname[MaxCharactersInLine];
    if(Type < 0){
        Snprintf(fname,"%s/AGBTests.dat",OutDir);
    } else {
        Snprintf(fname,"%s/AGBTests.dat.%02d",OutDir,Type);
    }
    FileOpen(fp,fname,"w");

    double tmpElements[CELibYield_Number];

    // Test 1: different SSP mass
    double MassArray[] = {1,10,100,1000,3000};
    int MassArraySize = sizeof(MassArray)/sizeof(double);
    // dprintl(MassArraySize);
    fprintf(fp,"IMF type : %s\n",CELibIMF[CELibRunParameters.IMFType].Name);
    fprintf(fp,"Test 1: Different SSP mass\n");
    struct CELibStructFeedbackOutput AGBFeedback;
    for(int i=0;i<MassArraySize;i++){
        for(int k=0;k<CELibYield_Number;k++)
            AGBFeedback.Elements[k] = 0.e0;

        for(int j=0;j<CELibRunParameters.AGBBinNumber;j++){
            struct CELibStructFeedbackOutput TempAGBFeedback =
                CELibGetAGBFeedback((struct CELibStructFeedbackInput){
                        .Mass = MassArray[i],
                        .MassConversionFactor = 1.0,
                        .Metallicity = 0.02,
                        .Count = j,
                        .Elements = tmpElements,
                        });
            for(int k=0;k<CELibYield_Number;k++)
                AGBFeedback.Elements[k] += TempAGBFeedback.Elements[k];
        }
        fprintf(fp," SSP mass = %g Msun\n",MassArray[i]);
        double Sum = 0.e0;
        for(int k=0;k<CELibYield_Number;k++){
            fprintf(fp," %g %s\n",AGBFeedback.Elements[k],CELibAGBYieldsElementName[k]);
            Sum += AGBFeedback.Elements[k];
        }
        fprintf(fp," Total released mass = %g\n",Sum);
        fprintf(fp,"\n");

    }

    fprintf(fp,"\n\n");
    // Test 2: different metallicity
    double MetallicityArray[] = {0,0.0005,0.005,0.02,0.05};
    int MetallicityArraySize = sizeof(MetallicityArray)/sizeof(double);
    fprintf(fp,"Test 2: Different metallicity\n");
    for(int i=0;i<MetallicityArraySize;i++){
        for(int k=0;k<CELibYield_Number;k++)
            AGBFeedback.Elements[k] = 0.e0;

        for(int j=0;j<CELibRunParameters.AGBBinNumber;j++){
            struct CELibStructFeedbackOutput TempAGBFeedback =
                CELibGetAGBFeedback((struct CELibStructFeedbackInput){
                        .Mass = 1.0,
                        .MassConversionFactor = 1.0,
                        .Metallicity = MetallicityArray[i],
                        .Count = j,
                        .Elements = tmpElements,
                        });
            for(int k=0;k<CELibYield_Number;k++)
                AGBFeedback.Elements[k] += TempAGBFeedback.Elements[k];
        }
        double Sum = 0.e0;
        for(int k=0;k<CELibYield_Number;k++){
            fprintf(fp," %g %s\n",AGBFeedback.Elements[k],CELibAGBYieldsElementName[k]);
            Sum += AGBFeedback.Elements[k];
        }
        fprintf(fp," Total released mass = %g\n",Sum);
        fprintf(fp,"\n");
    }

    fprintf(fp,"\n\n");

    // Test 3: different explosion time
    fprintf(fp,"Test 3: Different explosion time\n");
    for(int i=0;i<CELibRunParameters.AGBBinNumber+1;i++){
        struct CELibStructFeedbackOutput AGBFeedback =
            CELibGetAGBFeedback((struct CELibStructFeedbackInput){
                    .Mass = 1.0,
                    .Metallicity = 0.02,
                    .MassConversionFactor = 1.0,
                    .Count = i,
                    .Elements = tmpElements,
                    });
        fprintf(fp," Count = %d\n",i);
        for(int k=0;k<CELibYield_Number;k++)
            fprintf(fp," %g %s\n",AGBFeedback.Elements[k],CELibAGBYieldsElementName[k]);
        fprintf(fp,"\n");
    }


    fclose(fp);

    return ;
}

/*!
 * This function initializes all functions used for the AGB feedback.
 */
void CELibInitAGBYields(void){

    if(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_K10){
        CELibAGBYields_Mass = CELibAGBYields_Mass_K10;
        if(CELibRunParameters.PopIIIAGB == 0){
            CELibAGBYields_Metallicity = CELibAGBYields_Metallicity_K10;
            pCELibAGBYields = CELibAGBYieldsK10[0]; 
        } else {
            CELibAGBYields_Metallicity = CELibAGBYields_Metallicity_K10+1;
            pCELibMakeAGBYieldsCombinedYieldsTablePopIII();
            pCELibAGBYields = pCELibAGBYieldsCombined; 
        }
    } else if(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_K10D14){
        CELibAGBYields_Mass = CELibAGBYields_Mass_K10+CELibAGBYields_Mass_D14;
        if(CELibRunParameters.PopIIIAGB == 0){
            CELibAGBYields_Metallicity = CELibAGBYields_Metallicity_K10;
            pCELibMakeAGBYieldsCombinedYieldsTable();
        } else {
            CELibAGBYields_Metallicity = CELibAGBYields_Metallicity_K10+1;
            pCELibMakeAGBYieldsCombinedYieldsTablePopIII();
        }
        pCELibAGBYields = pCELibAGBYieldsCombined; 
    } else if(CELibRunParameters.AGBYieldsTableID == CELibAGBYieldsTableID_vdHG97){
        CELibAGBYields_Mass = CELibAGBYields_Mass_vdHG97;
        CELibAGBYields_Metallicity = CELibAGBYields_Metallicity_vdHG97;
        pCELibAGBYields = CELibAGBYieldsvdHG97[0]; 
        CELibRunParameters.PopIIIAGB = 0;
    }

    pCELibSetAGBYieldsElementName();

    pCELibGetAGBYieldsLookUpTables(); 
    pCELibAllocateAGBYieldsAGBYieldsIntegrated();
    pCELibMakeAGBYieldsIntegratedTable();
    pCELibMakeAGBYieldsIntegratedBinTable();

    if(CELibRunParameters.TestMode == true){
        pCELibDumpAGBYieldsTable("./CELib/CELibAGB");
        pCELibWriteAGBYieldsIntegratedTable("./CELib/CELibAGB");
        pCELibWriteAGBYieldsIntegratedBinTable("./CELib/CELibAGB");
        pCELibCheckAGBYieldsFeedBack("./CELib/CELibAGB",-1);
        pCELibWriteAGBYieldsReturnMassFraction("./CELib/CELibAGB");
    }

    return ;
}


/*!
 * This function returns the next mass loss time based on the random value of ``Rate'',
 * ``metalicity'', and the count of AGB ``Count''.
 */
double CELibGetAGBFeedbackTime(const double Rate, const double Metallicity, const int Count){

    if(Count >= CELibRunParameters.AGBBinNumber){
        return 10*CELibRunParameters.AGBBinUpperAge;
    }

    double TimeInterval = CELibAGBYieldsIntegratedBinInfo[0][Count].AgeUpper-CELibAGBYieldsIntegratedBinInfo[0][Count].AgeLower;
    return TimeInterval*Rate + CELibAGBYieldsIntegratedBinInfo[0][Count].AgeLower; 
}


/*!
 * This function returns the released mass of each element according to the
 * input parametes.
 */
struct CELibStructFeedbackOutput CELibGetAGBFeedback(struct CELibStructFeedbackInput Input){

    struct CELibStructFeedbackOutput AGBFeedback;

    AGBFeedback.Energy = 0.e0;
    if((Input.Count >= CELibRunParameters.AGBBinNumber)||(Input.Mass==0.e0)){
        AGBFeedback.EjectaMass = 0.e0;
        AGBFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor;
        for(int i=0;i<CELibYield_Number;i++)
            AGBFeedback.Elements[i] = 0.e0;
        return AGBFeedback;
    }


    if(Input.Metallicity <= CELibAGBYieldsMetallicityLookUpTable[0]){
        int TableID = 0;
        if(Input.noPopIII == 1){
            TableID ++;
        }
        for(int i=0;i<CELibYield_Number;i++){
            AGBFeedback.Elements[i] = CELibAGBYieldsIntegratedBin[TableID][Input.Count].Elements[i]
                +(Input.Elements[i]/Input.Mass)*CELibAGBYieldsIntegratedBin[TableID][Input.Count].Erf;
        }
    } else if (Input.Metallicity >= CELibAGBYieldsMetallicityLookUpTable[CELibAGBYields_Metallicity-1]){
        for(int i=0;i<CELibYield_Number;i++){
            AGBFeedback.Elements[i] = CELibAGBYieldsIntegratedBin[CELibAGBYields_Metallicity-1][Input.Count].Elements[i]
                +(Input.Elements[i]/Input.Mass)*CELibAGBYieldsIntegratedBin[CELibAGBYields_Metallicity-1][Input.Count].Erf;
        }
    } else {
        int IndexMetallicity = 0;
        for(int i=1;i<CELibAGBYields_Metallicity;i++){
            if(CELibAGBYieldsMetallicityLookUpTable[i] > Input.Metallicity){
                IndexMetallicity = i-1;
                break;
            }
        }
        double MetallicityEdges[] = {CELibAGBYieldsMetallicityLookUpTable[IndexMetallicity], 
                                     CELibAGBYieldsMetallicityLookUpTable[IndexMetallicity+1]};

        double GradErfLoss = (CELibAGBYieldsIntegratedBin[IndexMetallicity+1][Input.Count].Erf-CELibAGBYieldsIntegratedBin[IndexMetallicity][Input.Count].Erf)/
                            (MetallicityEdges[1]-MetallicityEdges[0]);

        double Erf = GradErfLoss*(Input.Metallicity-MetallicityEdges[0])+CELibAGBYieldsIntegratedBin[IndexMetallicity][Input.Count].Erf;

        for(int i=0;i<CELibYield_Number;i++){
            double GradYields = (CELibAGBYieldsIntegratedBin[IndexMetallicity+1][Input.Count].Elements[i]-
                                CELibAGBYieldsIntegratedBin[IndexMetallicity][Input.Count].Elements[i])/
                                (MetallicityEdges[1]-MetallicityEdges[0]);

            AGBFeedback.Elements[i] = (GradYields*(Input.Metallicity-MetallicityEdges[0]) + 
                                        CELibAGBYieldsIntegratedBin[IndexMetallicity][Input.Count].Elements[i])
                                    +(Input.Elements[i]/Input.Mass)*Erf;
        }
    }


    AGBFeedback.EjectaMass = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        AGBFeedback.Elements[i] *= Input.Mass*Input.MassConversionFactor;
        AGBFeedback.EjectaMass += AGBFeedback.Elements[i];
    }
    AGBFeedback.RemnantMass = Input.Mass*Input.MassConversionFactor - AGBFeedback.EjectaMass;

    return AGBFeedback;
}

