#pragma once 

struct CELibStructAGBYields{
    double Metallicity; // Initial metallicity.
    double Mass;        // Initial mass.
    double Mr;          // Remnant mass.
    double Mej;         // Ejecta mass.
    double Erf;         // Return mass fraction.
    double Energy;      // Energy.
    double Elements[CELibYield_Number];
};

void CELibInitAGBYields(void);

