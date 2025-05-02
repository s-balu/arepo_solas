#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef USE_GRACKLE
#include <grackle.h>

//Initialize Grackle
void InitGrackle(void)
{
    // Check the consistency
    if (gr_check_consistency() != GR_SUCCESS) {
      terminate("GRACKLE: Error in gr_check_consistency.\n");
    }
    
    int grackle_verbose = 0;
    // Enable output
    if(ThisTask == 0) grackle_verbose = 1;
    
    // First, set up the units system.
    // These are conversions from code units to cgs.
    All.GrackleUnits.comoving_coordinates = 0; //All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
    All.GrackleUnits.density_units        = All.UnitDensity_in_cgs;
    All.GrackleUnits.length_units         = All.UnitLength_in_cm;
    All.GrackleUnits.time_units           = All.UnitTime_in_s;
    All.GrackleUnits.velocity_units       = All.UnitVelocity_in_cm_per_s;
    All.GrackleUnits.a_units              = 1.0; // units for the expansion factor; NOTE: Should be 1 always
    set_velocity_units(&All.GrackleUnits);
    
    // Set initial expansion factor (for internal units).
    // Set expansion factor to 1 for non-cosmological simulation.
    double a_value = 1.0; 
    if(All.ComovingIntegrationOn) 
      {a_value = All.TimeBegin;}
    All.GrackleUnits.a_value = a_value;
    
    // Second, create a chemistry object for parameters and rate data.
    chemistry_data *my_grackle_data;
    my_grackle_data = malloc(sizeof(chemistry_data));
    
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
      terminate("\nGRACKLE: Error in set_default_chemistry_parameters.\n");
    }
    
    // Third, set parameter values for chemistry & cooling
    int three_body_rate=0, metal_cooling=0, h2_on_dust=0, photoelectric_heating=0, cmb_temperature_floor=1, UVbackground=1, Compton_xray_heating=1, cie_cooling=0, h2_optical_depth_approximation=0;
    int LWbackground_intensity=0, LWbackground_sawtooth_suppression=0, use_grackle=1, with_radiative_cooling=1, primordial_chemistry=0, dust_chemistry=0, H2_self_shielding=1, self_shielding_method=3; 
    double photoelectric_heating_rate=8.5e-26, Gamma;

    /* optional flags:: */
    three_body_rate = 0; /* Flag to control which three-body H2 formation rate is used.
                          0: Abel, Bryan & Norman (2002), 1: Palla, Salpeter & Stahler (1983), 2: Cohen & Westberg (1983), 3: Flower & Harris (2007), 4: Glover (2008).
                          These are discussed in Turk et. al. (2011). Default: 0. */
    cmb_temperature_floor = 1; // Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate. Default: 1.
    UVbackground = 1; // Flag to enable a UV background. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    Compton_xray_heating = 1; // Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
    cie_cooling = 0; // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0.
    h2_optical_depth_approximation = 0; // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0
    LWbackground_intensity = 0; // Rad_Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field, in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
    LWbackground_sawtooth_suppression = 0; // Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
    H2_self_shielding=0; // Switch to enable approximate H2 self-shielding from both the UV background dissociation rate and the H2 dissociation rate using local Sobolev approximation (see user guide for details)
    self_shielding_method=3; // Switch to enable approximate self-shielding from the UV background (see user guide for details)
#ifdef METALS
    metal_cooling = 1; // Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    dust_chemistry = 1; // Flag to control additional dust cooling and chemistry processes. Default: 0 (no dust). 1: adds the following processes: photo-electric heating, electron recombination on dust, H2 formation on dust.
    h2_on_dust = 0; // Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the dust to gas ratio scales with the metallicity. Default: 0.
    // Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
    // If photoelectric_heating enabled, photoelectric_heating_rate is the heating rate in units of erg cm-3 s-1. Default: 8.5e-26.
    //   (Caution: this tends to heat gas even at extremely high densities to ~3000 K, when it should be entirely self-shielding)
    photoelectric_heating = 1; // photo-electric on [but not adjusted to local background, beware!]
    photoelectric_heating_rate = 8.5e-26; // default rate normalization
#endif
    /* fixed flags:: */
    use_grackle = 1; // Flag to activate the grackle machinery
    with_radiative_cooling = 1; // Flag to include radiative cooling and actually update the thermal energy during the chemistry solver. If off, the chemistry species will still be updated. The most common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
    Gamma = GAMMA; // The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3.
    primordial_chemistry = 0; // Flag to control which primordial chemistry network is used (set by Config file). 0 = fully tabulated cooling
#ifdef GRACKLE_CHEMISTRY
    primordial_chemistry = GRACKLE_CHEMISTRY;
#endif


    my_grackle_data->grackle_data_file = All.GrackleDataFile; // Path to the data file containing the metal cooling and UV background tables
    my_grackle_data->three_body_rate = three_body_rate;
    my_grackle_data->metal_cooling = metal_cooling;
    my_grackle_data->h2_on_dust = h2_on_dust;
    my_grackle_data->photoelectric_heating = photoelectric_heating;
    my_grackle_data->photoelectric_heating_rate = photoelectric_heating_rate;
    my_grackle_data->cmb_temperature_floor = cmb_temperature_floor;
    my_grackle_data->UVbackground = UVbackground;
    my_grackle_data->Compton_xray_heating = Compton_xray_heating;
    my_grackle_data->cie_cooling = cie_cooling;
    my_grackle_data->h2_optical_depth_approximation = h2_optical_depth_approximation;
    my_grackle_data->LWbackground_intensity = LWbackground_intensity;
    my_grackle_data->LWbackground_sawtooth_suppression = LWbackground_sawtooth_suppression;
    my_grackle_data->use_grackle = use_grackle;
    my_grackle_data->with_radiative_cooling = with_radiative_cooling;
    my_grackle_data->Gamma = Gamma;
    my_grackle_data->primordial_chemistry = primordial_chemistry;
    my_grackle_data->dust_chemistry = dust_chemistry;
    my_grackle_data->H2_self_shielding = H2_self_shielding;
    my_grackle_data->self_shielding_method = self_shielding_method;
    
   
    // Create chemistry data storage object to store rates.
    chemistry_data_storage my_grackle_rates;

    // Finally, initialize the chemistry object.
    if (local_initialize_chemistry_data(my_grackle_data, &my_grackle_rates, &All.GrackleUnits) == 0) {
      terminate("Error in initialize_chemistry_data.\n");
    }

    if(ThisTask == 0) {
    printf("\nGRACKLE: Grackle Initialized.\n\n");}
}
#endif  //USE_GRACKLE
