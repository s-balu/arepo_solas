#pragma once

struct CELibStructSNIIYields{
    double Metallicity; // Initial metallicity.
    double Mass; // Initial mass.
    double Mr;  // Remnant mass.
    double Msn; // Mass released via supernova
    double Msw; // Mass released via stellar wind
    double Energy; // Energy by a single SN
    double Elements[CELibYield_Number];
};

extern int CELibSNIIYields_Metallicity;
extern int CELibSNIIYields_Mass;
extern double *CELibSNIIYieldsZ; 
extern struct CELibStructSNIIYields *CELibSNIIYieldsIntegrated;

void CELibInitSNIIYields(void);
