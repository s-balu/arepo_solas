/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/star_formation/sfr_eEOS.c
 * \date        05/2018
 * \brief       Star formation rate routines for the effective multi-phase
 *              model.
 * \details     contains functions:
 *                void cooling_and_starformation(void)
 *                double get_starformation_rate(int i)
 *                void init_clouds(void)
 *                void integrate_sfr(void)
 *                void set_units_sfr(void)
 *                double calc_egyeff(int i, double gasdens, double *ne,
 *                  double *x, double *tsfr, double *factorEVP)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../gravity/forcetree.h"

#if defined(USE_SFR) && defined(JEANS_SF) && !defined(AGORA_SF)

/*! \brief Return the Jeans length of the cell.
 *
 *  \param[in] i the index of the gas cell.
 *
 *  \return Jeans length in code units.
 */
double get_jeans_length(int i)
{
  double sound_speed, jeans_length;

  sound_speed  = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
  
  jeans_length = sqrt(M_PI / All.G / SphP[i].Density) * sound_speed;

  return jeans_length;
}

/*! \brief Return the Jeans Mass of the cell.
 *
 *  \param[in] i the index of the gas cell.
 *
 *  \return Jeans mass in code units.
 */
double get_jeans_mass(int i)
{
  double sound_speed, jeans_mass;

  sound_speed  = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
  
  jeans_mass = pow(M_PI, 2.5) * pow(sound_speed, 3.) / 6. / pow(All.G, 1.5) / sqrt(SphP[i].Density);

  return jeans_mass;
}

/*! \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the
 *  function cool_cell() and the star formation rate is set to 0.
 *
 *  \return void
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin, flag;
  double dt, dtime, ne = 1;
  double unew, du;

  double t_freefall;
#ifdef JEANS_MASS_BASED
  double jeans_mass, jeans_mass_thresh;
#else 
  double cell_size, jeans_length;
#endif

  /* note: assuming FULL ionization */
  double u_to_temp_fac =
      (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      dt    = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      /* apply the temperature floor */

      unew = dmax(All.MinEgySpec, SphP[i].Utherm);

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
      
      cool_cell(i);

      /* check whether conditions for star formation are fulfilled.
       * f=1  normal cooling
       * f=0  star formation
       */

      flag = 1; /* default is normal cooling */

#ifdef JEANS_MASS_BASED
    {
      jeans_mass = get_jeans_mass(i);
      /* enable star formation if jeans mass is smaller than factor x gas cell mass */
      if(jeans_mass < All.JeansMassThreshold * P[i].Mass)
        flag = 0;
    }
#else 
      cell_size = 2.0 * get_cell_radius(i);
      jeans_length = get_jeans_length(i);     
      /* enable star formation if jeans length is smaller than the gas cell, i.e. we are not resolving the jeans length */
      if(jeans_length < cell_size)
        flag = 0;
#endif

      if(P[i].Mass == 0) /* tracer particles don't form stars */
        flag = 1;

      if(flag == 1)
        SphP[i].Sfr = 0;

      /* active star formation */
      if(flag == 0)
        {
          SphP[i].Ne = (HYDROGEN_MASSFRAC + 1) / 2 / HYDROGEN_MASSFRAC; /* note: assuming FULL ionization */

          if(dt > 0)
            {
              if(P[i].TimeBinHydro) /* upon start-up, we need to protect against dt==0 */
                {
                  unew = SphP[i].Utherm;

                  du = unew - SphP[i].Utherm;
                  if(unew < All.MinEgySpec)
                    du = All.MinEgySpec - SphP[i].Utherm;

                  SphP[i].Utherm += du;
                  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
                  if(dtime > 0)
                    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

                  set_pressure_of_cell(i);
                }
            }

            //  Compute the freefall time
            t_freefall=sqrt(3.*M_PI/32/All.G/SphP[i].Density); // freefall time in code units

            SphP[i].Sfr=All.StarFormationEfficiency*P[i].Mass/t_freefall * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
            
            TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
    } /* end of main loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

/*! \brief Return the star formation rate associated with the gas cell i.
 *
 *  \param[in] i the index of the gas cell.
 *
 *  \return star formation rate in solar masses / yr.
 */
double get_starformation_rate(int i)
{
  if(RestartFlag == 3)
    return SphP[i].Sfr;

  double rateOfSF;
  int flag;
  double tsfr;
  double factorEVP, egyeff, ne, x, cloudmass;
  
  double t_freefall;
#ifdef JEANS_MASS_BASED
  double jeans_mass, jeans_mass_thresh;
#else 
  double cell_size, jeans_length;
#endif

  /* note: assuming FULL ionization */
  double u_to_temp_fac =
      (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  flag   = 1; /* default is normal cooling */

#ifdef JEANS_MASS_BASED
    {
      jeans_mass = get_jeans_mass(i);
      /* enable star formation if jeans mass is smaller than factor x gas cell mass */
      if(jeans_mass < All.JeansMassThreshold * P[i].Mass)
        flag = 0;
    }
#else 
      cell_size = 2.0 * get_cell_radius(i);
      jeans_length = get_jeans_length(i);     
      /* enable star formation if jeans length is smaller than the gas cell, i.e. we are not resolving the jeans length */
      if(jeans_length < cell_size)
        flag = 0;
#endif
  
  if(flag == 1)
    return 0;

  t_freefall=sqrt(3.*M_PI/32/All.G/SphP[i].Density); // freefall time in code units

  rateOfSF=All.StarFormationEfficiency*P[i].Mass/t_freefall;

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

/*! \brief Set the appropriate units for the parameters of the multi-phase
 *         model.
 *
 *  \return void
 */
void set_units_sfr(void)
{
  double meanweight;

  All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC); /* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
}

#endif /* #if defined(USE_SFR) && defined(JEANS_SF) && !defined(AGORA_SF) */
