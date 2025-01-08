#include "config.h"
#include "IMFFunctions.h"

#define NBin    (100)

/*!
 * This function writes information of all pre-implemented IMFs.
 */
void CELibWriteIMFData(const int Type, const char OutputDir[]){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./%s/%s.dat",OutputDir,CELibIMF[Type].Name);

    FileOpen(fp,fname,"w");

    fprintf(fp,"#IMFName:%s\n",CELibIMF[Type].Name);
    fprintf(fp,"#MassRange:%g-%g\n",CELibIMF[Type].MassMin,CELibIMF[Type].MassMax);
    fprintf(fp,"#IMFDomainNumber:%d\n",CELibIMF[Type].Ndomain);
    for(int i=0;i<CELibIMF[Type].Ndomain;i++){
        fprintf(fp,"#%g-%g #%g\n",CELibIMF[Type].Lmass[i],CELibIMF[Type].Umass[i],CELibIMF[Type].Power[i]);
    }

    fclose(fp);

    return ;
}


/*!
 * This function writes cumlative mass functions of all pre-implemented IMFs.
 */
void CELibWriteIMFCumlative(const char OutputDir[]){

    MakeDir(OutputDir);

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./%s/Cumulative.dat",OutputDir);

    FileOpen(fp,fname,"w");

    fprintf(fp,"#Mass ");
    for(int i=0;i<CELibIMF_NTypes;i++){
        fprintf(fp,"#%s ",CELibIMF[i].Name);
    }
    fprintf(fp,"\n");

    int Nbin = 100;
    double dmass = (log10(120)-log10(0.1))/Nbin;
    for(int i=0;i<Nbin;i++){
        double Mu = pow(10.0,log10(120)-dmass*i);
        double Ml = pow(10.0,log10(120)-dmass*(i+1));

        fprintf(fp,"%g %g %g ",0.5*(Ml+Mu),Ml,Mu);
        for(int k=0;k<CELibIMF_NTypes;k++){
            if(Mu > CELibIMF[k].MassMax) Mu = CELibIMF[k].MassMax;
            if(Ml > CELibIMF[k].MassMax) Ml = CELibIMF[k].MassMax;
            if(Mu < CELibIMF[k].MassMin) Mu = CELibIMF[k].MassMin;
            if(Ml < CELibIMF[k].MassMin) Ml = CELibIMF[k].MassMin;
            fprintf(fp,"%g ",IntegralSimpson(Ml,CELibIMF[k].MassMax,CELibRunParameters.IntegrationSteps,CELibIMF[k].IMFFunction));
        }
        fprintf(fp,"\n");

    }

    fclose(fp);

    return ;
}
