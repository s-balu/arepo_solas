#pragma once

#define MaxCharacterForIMFName  (20)
struct CELibStructIMF{
    int Ndomain;     //!< Number of domains for the target IMF.
    double MassMin;  //!< Minimum mass of the IMF.
    double MassMax;  //!< Maximum mass of the IMF.
    double Power[3]; //!< Power law index of each domain.
    double Lmass[3]; //!< Lower mass of the i-th domain's power law index.
    double Umass[3]; //!< Upper mass of the i-th domain's power law index.
    double Sigma;            //!< Mass dispersion used in log-normal type IMFs. 
    double CriticalMass;     //!< Critical mass used in log-normal type IMFs.
    double Normalization[3]; //!< Normalization of each domain.
    double (*IMFFunction)(const double);        //!< Function pointer to the functional form of the IMF.
    double (*IMFFunctionPerMass)(const double); //!< Function pointer to the functional form of the IMF(x)/x.
    double RSNII;                      //!< Type II rate (will be removed).
    char Name[MaxCharacterForIMFName]; //!< IMF name.
};

extern struct CELibStructIMF CELibIMF[CELibIMF_NTypes];
