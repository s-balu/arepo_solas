/*! \file CELib.h
 *  \brief This file has all prototype declarations.
 */

/*!
 * @enum CELib Lifetime types
 * Lifetime types.
 */
enum {
    CELibLifeTime_Original, // Original data with linear interpolation.
    CELibLifeTime_LSF,      // Least square fitting
    CELibLifeTime_Number,
};


/*!
 * @enum CELib SN II yields table IDs
 * SN II yields table IDs.
 */
enum {
    CELibSNIIYieldsTableID_P98, // Portinari et al. 1998
    CELibSNIIYieldsTableID_Number,
};

/*!
 * @enum CELib feedback types
 * Feedback types.
 */
enum {
    CELibFeedbackType_SNII,
    CELibFeedbackType_Number,
};

/*!
 * @enum CELibYields
 * Elements used in CELib.
 */
enum {
    CELibYield_H,      //! Hydrogen
    CELibYield_He,     //! Helium
    CELibYield_C,      //! Carbon
    CELibYield_N,      //! Nitrogen
    CELibYield_O,      //! Oxygen
    CELibYield_Ne,     //! Neon
    CELibYield_Mg,     //! Magnesium
    CELibYield_Si,     //! Silicon
    CELibYield_S,      //! Sulfur
    CELibYield_Ca,     //! Calcium
    CELibYield_Fe,     //! Iron
    CELibYield_Ni,     //! Nickel
    CELibYield_Eu,     //! Europium
    CELibYield_Number, //! Total number of elements.
};

/*! @struct CELibStructRunParameters
 * This structure contains all control parameters of CELib.
 */
struct CELibStructRunParameters{
    int IntegrationSteps;           //!< Number of integration steps
    int LifeTimeType;               //!< Type of stellar life time.

    /* Parameters for Type II SNe */
    int SNIIYieldsTableID;           //!< Yields table ID. 0=Portinari+1998

};

extern struct CELibStructRunParameters CELibRunParameters;

/*! @struct CELibStructFeedbackStarbyStarOutput
 * This structure contains all quantities evaluated by the star-by-star feedback functions.
 */
struct CELibStructFeedbackStarbyStarOutput{
    double Energy;             // Energy in erg.
    double EjectaMass;         // Total ejecta mass in Msun.
    double RemnantMass;        // Remnant mass in Msun.
    double Elements[CELibYield_Number]; // Return mass to the ISM (Y, not y) in Msun.
};

/*! @struct CELibStructFeedbackStarbyStarInput
 * This structure contains all parameters used in the star-by-star feedback.
 */
struct CELibStructFeedbackStarbyStarInput{
    double Mass;                 // Mass of the star in simulation unit.
    double Metallicity;          // Metallicity of the star Z
    double MassConversionFactor; // A factor to convert Elements[] from the simulation mass unit to Msun.
    double *Elements;            // Star particle's elements composition in simulation unit.
};

/*! @struct CELibStructNextEventTimeStarbyStarInput
 * This structure contains all parameters used to evaluate next event time for
 * the star-by-star case.
 */
struct CELibStructNextEventTimeStarbyStarInput{
    double InitialMass_in_Msun;  // Mass of the star
    double Metallicity;          // Metallicity of the star, Z
};

//_______Functions for stellar lifetime in LifeTime.c
#define CELIB_LIFETIME_Z_P98  5
#define CELIB_LIFETIME_M_P98  30

#define CELIB_LIFETIME_Z_S02  1
#define CELIB_LIFETIME_M_S02  4

extern double *CELibLifeTimeZ;
extern double *CELibLifeTimeLogZ;
extern double CELibLifeTimeMass[CELIB_LIFETIME_M_P98];
extern double CELibLifeTimeZMF[CELIB_LIFETIME_Z_P98][CELIB_LIFETIME_M_P98];

void CELibInitLifeTime(void);
double CELibGetLifeTimeofStarTable(const double Mass, const double Metallicity);
double CELibGetDyingStellarMassTable(const double LifeTime, const double Metallicity);
double CELibGetLifeTimeofStarLSF(const double Mass, const double Metallicity);
double CELibGetDyingStellarMassLSF(const double LifeTime, const double Metallicity);
int CELibGetLifeTimeTableSizeMetallicity(void);

//_______Functions in SNIIYields.c
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

void CELibInitSNIIYields(void);

extern double CELibSetSNIIEmptyElements[CELibYield_Number];
struct CELibStructFeedbackStarbyStarOutput CELibGetSNIIFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input);

//_______Functions in UnifiedAPIs.c
struct CELibStructFeedbackStarbyStarOutput CELibGetFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input, const int Type);
double CELibGetNextEventTimeStarbyStar(struct CELibStructNextEventTimeStarbyStarInput Input, const int Type);