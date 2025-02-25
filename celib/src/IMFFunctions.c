#include "config.h"
#include "Integral.h"
#include "IMFFunctions.h"


static bool FirstCall = true; //!< First call flag.


/*!
 * The function of the Salpeter IMF is \phi(m) = m^{-1.35}, where \phi(m) =
 * dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */
static double pCELibGetSalpeterIMF(const double x){
    return CELibIMF[CELibIMF_Salpeter].Normalization[0]*pow(x,CELibIMF[CELibIMF_Salpeter].Power[0]);
}


/*!
 * This function returns the Salpeter IMF divided by mass, i.e.,
 * \phi(m)_Salpeter/m.
 */
static double pCELibGetSalpeterIMFPerMass(const double x){
    return pCELibGetSalpeterIMF(x)/x;
}


/*!
 * The function of the Diet Salpeter IMF is 
 * \phi(m) = m^{-0} for 0.1<m<0.6, and
 *           m^{-1.35} for 0.6<m<120
 * \phi(m) = dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */
static double pCELibGetDietSalpeterIMF(const double x){
    
    if(x>CELibIMF[CELibIMF_DietSalpeter].Lmass[1]){
        return CELibIMF[CELibIMF_DietSalpeter].Normalization[1]*pow(x,CELibIMF[CELibIMF_DietSalpeter].Power[1]);
    } else {
        return CELibIMF[CELibIMF_DietSalpeter].Normalization[0]*pow(x,CELibIMF[CELibIMF_DietSalpeter].Power[0]);
    }

}


/*!
 * This function returns \phi(m)_Diet-Salpeter/m.
 */
static double pCELibGetDietSalpeterIMFPerMass(const double x){
    return pCELibGetDietSalpeterIMF(x)/x;
}


/*!
 * The function of the Miller-Scalo IMF is 
 * \phi(m) = m^{-0.4} for 0.1<m<1.0,
 *           m^{-1.5} for 1.0<m<10, and
 *           m-{-2.3} for 10<m<120.
 * \phi(m) = dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */
static double pCELibGetMillerScaloIMF(const double x){

    if(x>CELibIMF[CELibIMF_MillerScalo].Lmass[2]){
        return CELibIMF[CELibIMF_MillerScalo].Normalization[2]*pow(x,CELibIMF[CELibIMF_MillerScalo].Power[2]);
    } else if(x>CELibIMF[CELibIMF_MillerScalo].Lmass[1]){
        return CELibIMF[CELibIMF_MillerScalo].Normalization[1]*pow(x,CELibIMF[CELibIMF_MillerScalo].Power[1]);
    } else {
        return CELibIMF[CELibIMF_MillerScalo].Normalization[0]*pow(x,CELibIMF[CELibIMF_MillerScalo].Power[0]);
    }

}


/*!
 * This function returns \phi(m)_Miller-Scalo/m.
 */
static double pCELibGetMillerScaloIMFPerMass(const double x){
    return pCELibGetMillerScaloIMF(x)/x;
}


/*!
 * The function of the Kroupa (2001) IMF is 
 * \phi(m) = m^{-0.3} for 0.1<m<0.5,
 *           m^{-1.3} for 0.5<m<120.
 * \phi(m) = dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */

static double pCELibGetKroupaIMF(const double x){
    
    if(x>CELibIMF[CELibIMF_Kroupa].Lmass[1]){
        return CELibIMF[CELibIMF_Kroupa].Normalization[1]*pow(x,CELibIMF[CELibIMF_Kroupa].Power[1]);
    } else {
        return CELibIMF[CELibIMF_Kroupa].Normalization[0]*pow(x,CELibIMF[CELibIMF_Kroupa].Power[0]);
    }

}


/*!
 * This function returns \phi(m)_Kroupa/m.
 */
static double pCELibGetKroupaIMFPerMass(const double x){
    return pCELibGetKroupaIMF(x)/x;
}


/*!
 * The function of the Kroupa (1993) IMF is 
 * \phi(m) = m^{-0.3} for 0.1<m<0.5,
 *           m^{-1.2} for 0.5<m<1.0, and
 *           m-{-1.7} for 1.0<m<120.
 * \phi(m) = dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */

static double pCELibGetKroupa1993IMF(const double x){
    
    if(x>CELibIMF[CELibIMF_Kroupa1993].Lmass[2]){
        return CELibIMF[CELibIMF_Kroupa1993].Normalization[2]*pow(x,CELibIMF[CELibIMF_Kroupa1993].Power[2]);
    }else if(x>CELibIMF[CELibIMF_Kroupa1993].Lmass[1]){
        return CELibIMF[CELibIMF_Kroupa1993].Normalization[1]*pow(x,CELibIMF[CELibIMF_Kroupa1993].Power[1]);
    } else {
        return CELibIMF[CELibIMF_Kroupa1993].Normalization[0]*pow(x,CELibIMF[CELibIMF_Kroupa1993].Power[0]);
    }

}


/*!
 * This function returns \phi(m)_Kroupa/m.
 */
static double pCELibGetKroupa1993IMFPerMass(const double x){
    return pCELibGetKroupa1993IMF(x)/x;
}


/*!
 * The function of the Kennicutt IMF is 
 * \phi(m) = m^{-0.4} for 0.1<m<1.0, and
 *           m^{-1.5} for 1.0<m<120.
 * \phi(m) = dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */
static double pCELibGetKennicuttIMF(const double x){
    
    if(x>CELibIMF[CELibIMF_Kennicutt].Lmass[1]){
        return CELibIMF[CELibIMF_Kennicutt].Normalization[1]*pow(x,CELibIMF[CELibIMF_Kennicutt].Power[1]);
    } else {
        return CELibIMF[CELibIMF_Kennicutt].Normalization[0]*pow(x,CELibIMF[CELibIMF_Kennicutt].Power[0]);
    }

}


/*!
 * This function returns \phi(m)_Kennicutt/m.
 */
static double pCELibGetKennicuttIMFPerMass(const double x){
    return pCELibGetKennicuttIMF(x)/x;
}


/*!
 * The functional form of the Chabrier IMF is 
 * \phi(m) = m^{-0.0}exp(-[log10(m/CriticalMass)]^2/(2*sigma^2)) for 0.1<m<1.0, and
 *           m^{-1.3} for 1.0<m<100.
 * \phi(m) = dN/d\log10(m) = m dN/dm \propto m^{-x}.
 */
static double pCELibGetChabrierIMF(const double x){
    
    if(x>CELibIMF[CELibIMF_Chabrier].Lmass[1]){
        return CELibIMF[CELibIMF_Chabrier].Normalization[1]*pow(x,CELibIMF[CELibIMF_Chabrier].Power[1]);
    } else {
        return CELibIMF[CELibIMF_Chabrier].Normalization[0]*pow(x,CELibIMF[CELibIMF_Chabrier].Power[0])*
            exp(-SQ(log10(x/CELibIMF[CELibIMF_Chabrier].CriticalMass))/(2*SQ(CELibIMF[CELibIMF_Chabrier].Sigma)));
    }

}


/*!
 * This function returns \phi(m)_Chabrier/m.
 */
static double pCELibGetChabrierIMFPerMass(const double x){
    return pCELibGetChabrierIMF(x)/x;
}


/*!
 * The functional form of the Susa IMF is 
 * \phi(m) = (A/sigma)*exp(-[log10(m/CriticalMass)]^2/(2*sigma^2)) for 0.7<m<300.
 */
static double pCELibGetSusaIMF(const double x){

    if(x < CELibIMF[CELibIMF_Susa].MassMin) return 0.e0;
    if(x > CELibIMF[CELibIMF_Susa].MassMax) return 0.e0;
    
    return CELibIMF[CELibIMF_Susa].Normalization[0]*
        exp(-SQ(log10(x/CELibIMF[CELibIMF_Susa].CriticalMass))/(2*SQ(CELibIMF[CELibIMF_Susa].Sigma)));
}


/*!
 * This function returns \phi(m)_Susa/m.
 */
static double pCELibGetSusaIMFPerMass(const double x){
    return pCELibGetSusaIMF(x)/x;
}


/*!
 * This function returns the value of the adopted IMF at the given mass of x.
 * The IMF type is given by the input parameter, IMFType.  Note that the mass of
 * x should be in the unit of solar mass.
 */
double CELibGetIMFValue(const double x, const int IMFType){
    
    if(IMFType == CELibIMF_Salpeter){
        return pCELibGetSalpeterIMF(x);
    } else if(IMFType == CELibIMF_DietSalpeter){
        return pCELibGetDietSalpeterIMF(x);
    } else if(IMFType == CELibIMF_MillerScalo){
        return pCELibGetMillerScaloIMF(x);
    } else if(IMFType == CELibIMF_Kroupa){
        return pCELibGetKroupaIMF(x);
    } else if(IMFType == CELibIMF_Kroupa1993){
        return pCELibGetKroupa1993IMF(x);
    } else if(IMFType == CELibIMF_Kennicutt){
        return pCELibGetKennicuttIMF(x);
    } else if(IMFType == CELibIMF_Chabrier){
        return pCELibGetChabrierIMF(x);
    } else if(IMFType == CELibIMF_Susa){
        return pCELibGetSusaIMF(x);
    }

    fprintf(stderr,"IMF type is wrong. %s:%d\n",__FUNCTION__,__LINE__);
    abort();
}

/*!
 * This function returns the value of the adopted IMF over mass at the given
 * mass of x.  The IMF type is given by the input parameter, IMFType.  Note that
 * the mass of x should be in the unit of solar mass.
 */
double CELibGetIMFValuePerMass(const double x, const int IMFType){
    return CELibGetIMFValue(x,IMFType)/x;
}

#define SHOWENUM(label) #label

/*!
 * This function computes normalizations of all IMFs. If TestMode == true, these
 * values are shown on the screen.
 */
static void pCELibSetIMFNormalization(const int IMFType){

    if(CELibRunParameters.TestMode){
        fprintf(stderr,"//\t Compute IMF normalization for %d.\n",IMFType);
    }

    if(IMFType == CELibIMF_Salpeter){

        CELibIMF[IMFType].IMFFunction = pCELibGetSalpeterIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetSalpeterIMFPerMass;

        CELibIMF[IMFType].Normalization[0] /= 
            (IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                                                  CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));
        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_DietSalpeter){

        CELibIMF[IMFType].IMFFunction = pCELibGetDietSalpeterIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetDietSalpeterIMFPerMass;

        for(int i=1;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] = CELibIMF[IMFType].Normalization[i-1]*
                            pow(CELibIMF[IMFType].Umass[i-1],CELibIMF[IMFType].Power[i-1])/
                            pow(CELibIMF[IMFType].Lmass[i],CELibIMF[IMFType].Power[i]);

        double Norm = 1.0/(IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                        CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));
        for(int i=0;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] *= Norm;

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_MillerScalo){

        CELibIMF[IMFType].IMFFunction = pCELibGetMillerScaloIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetMillerScaloIMFPerMass;

        for(int i=1;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] = CELibIMF[IMFType].Normalization[i-1]*
                            pow(CELibIMF[IMFType].Umass[i-1],CELibIMF[IMFType].Power[i-1])/
                            pow(CELibIMF[IMFType].Lmass[i],CELibIMF[IMFType].Power[i]);


        double Norm = 1.0/(IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                    CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));

        for(int i=0;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] *= Norm;

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_Kroupa){

        CELibIMF[IMFType].IMFFunction = pCELibGetKroupaIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetKroupaIMFPerMass;

        for(int i=1;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] = CELibIMF[IMFType].Normalization[i-1]*
                            pow(CELibIMF[IMFType].Umass[i-1],CELibIMF[IMFType].Power[i-1])/
                            pow(CELibIMF[IMFType].Lmass[i],CELibIMF[IMFType].Power[i]);


        double Norm = 1.0/(IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                    CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));

        for(int i=0;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] *= Norm;

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_Kroupa1993){

        CELibIMF[IMFType].IMFFunction = pCELibGetKroupa1993IMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetKroupa1993IMFPerMass;

        for(int i=1;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] = CELibIMF[IMFType].Normalization[i-1]*
                            pow(CELibIMF[IMFType].Umass[i-1],CELibIMF[IMFType].Power[i-1])/
                            pow(CELibIMF[IMFType].Lmass[i],CELibIMF[IMFType].Power[i]);

        double Norm = 1.0/(IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                    CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));

        for(int i=0;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] *= Norm;

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_Kennicutt){

        CELibIMF[IMFType].IMFFunction = pCELibGetKennicuttIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetKennicuttIMFPerMass;

        for(int i=1;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] = CELibIMF[IMFType].Normalization[i-1]*
                            pow(CELibIMF[IMFType].Umass[i-1],CELibIMF[IMFType].Power[i-1])/
                            pow(CELibIMF[IMFType].Lmass[i],CELibIMF[IMFType].Power[i]);


        double Norm = 1.0/(IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                    CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));

        for(int i=0;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] *= Norm;

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_Chabrier){

        CELibIMF[IMFType].IMFFunction = pCELibGetChabrierIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetChabrierIMFPerMass;

        for(int i=1;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] = CELibIMF[IMFType].Normalization[i-1]*
                pow(CELibIMF[IMFType].Umass[i-1],CELibIMF[IMFType].Power[i-1])*
                exp(-SQ(log10(CELibIMF[IMFType].Umass[i-1]/CELibIMF[IMFType].CriticalMass))/(2*SQ(CELibIMF[IMFType].Sigma)))/
                            pow(CELibIMF[IMFType].Lmass[i],CELibIMF[IMFType].Power[i]);

        double Norm = 1.0/(IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                    CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));

        for(int i=0;i<CELibIMF[IMFType].Ndomain;i++)
            CELibIMF[IMFType].Normalization[i] *= Norm;

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    } else if(IMFType == CELibIMF_Susa){

        CELibIMF[IMFType].IMFFunction = pCELibGetSusaIMF;
        CELibIMF[IMFType].IMFFunctionPerMass = pCELibGetSusaIMFPerMass;

        CELibIMF[IMFType].Normalization[0] /= 
            (IntegralSimpsonOneExtraParameterType(CELibIMF[IMFType].MassMin,CELibIMF[IMFType].MassMax,
                    CELibRunParameters.IntegrationSteps,&CELibGetIMFValue,IMFType));

        CELibIMF[IMFType].RSNII = IntegralSimpsonOneExtraParameterType(8.0,CELibIMF[IMFType].MassMax,
                CELibRunParameters.IntegrationSteps,&CELibGetIMFValuePerMass,IMFType);

    }

    return ;
}

/*!
 * This function computes normalizations of all IMFs. If TestMode == true, these
 * values are shown on the screen.
 */
static void pCELibInitIMFNormalization(void){
    for(int i=0;i<CELibIMF_NTypes;i++){
        pCELibSetIMFNormalization(i);
    }
    return ;
}

/*!
 * This function initializes all IMFs used in CELib.
 */
void CELibInitIMF(void){

    if(FirstCall == true){
        pCELibInitIMFNormalization();
        FirstCall = false;
    }

    return ;
}

/*!
 * This function sets the current IMF name to "IMFFname".
 */
void CELibGetIMFName(char *IMFName){
    Snprintf(IMFName,"%s",CELibIMF[CELibRunParameters.IMFType].Name);
    return ;
}

/*!
 * This function returns the current IMF Type ID.
 */
int CELibGetIMFID(void){
    return CELibRunParameters.IMFType;
}

