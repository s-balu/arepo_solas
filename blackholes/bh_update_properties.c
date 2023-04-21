#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

static int int_compare(const void *a, const void *b)

void update_bh_accretion_rate(void)
{
  if(NumBh == 0)
    return;

/*calculate bondi accretion rate*/
  int i;
  double density, pressure, sound_speed, velocity_gas_norm;
  double denominator, denominator_inv, BondiRate, EddingtonRate, accretion_rate;

  for(i = 0; i < NumBh; i++)
    {

/*get pressure*/
      if(BhP[i].Density>0)
        {  
          density = BhP[i].Density;
          pressure = GAMMA_MINUS1 * density * BhP[i].InternalEnergyGas;

/*get soundspeed*/
          sound_speed = sqrt(GAMMA * pressure / density);
      
          velocity_gas_norm = sqrt(BhP[i].VelocityGas[0]*BhP[i].VelocityGas[0] + BhP[i].VelocityGas[1]*BhP[i].VelocityGas[1] + BhP[i].VelocityGas[2]*BhP[i].VelocityGas[2]); 

          denominator = (sound_speed*sound_speed + velocity_gas_norm*velocity_gas_norm);
          if(denominator > 0)
            {
              denominator_inv = 1. / sqrt(denominator);
              BondiRate = 4. * M_PI * GRAVITY * GRAVITY * PPB(i).Mass * density *
                   denominator_inv * denominator_inv * denominator_inv;
            }
          else
            terminate("Invalid denominator in Bondi Accretion Rate");
        }
      else
        BondiRate = 0;
  
/*limit by Eddington accretion rate*/
      EddingtonRate = 4. * M_PI * GRAVITY * PPB(i).Mass * PROTONMASS / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *= All.UnitTime_in_s / All.UnitMass_in_g;
  
      accretion_rate = fmin(BondiRate, EddingtonRate);
  
/*efficiency*/
      accretion_rate *= (1. - All.Epsilon_r);

      BhP[i].AccretionRate = accretion_rate;
    }
}

/*update bh-timestep at prior_mesh_construction based on ngb smallest timestep*/
void update_bh_timesteps(void)
{
  if(NumBh == 0)
    return;
  
  int idx, i, binold, bin;
  integertime ti_step;

  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      ti_step = BPP(i).NgbMinStep;
      binold = P[i].TimeBinBh;
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      if(bin < binold || binold == 0) /*need the "or" condition for first loop for these ngb criteria*/
        {
          timebin_move_particle(&TimeBinsBh, i, binold, bin);
          P[i].TimeBinBh = bin;
        }
    }
}

/*call this function after updating the bh-timebin to the ngb condition*/
void update_list_of_active_bh_particles(void)
{
  int i, n;
  TimeBinsBh.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsBh.FirstInTimeBin[n]; i >= 0; i = TimeBinsBh.NextInTimeBin[i])
            {
              if(P[i].Type == 5)
                {
                  if(P[i].Ti_Current != All.Ti_Current)
                    drift_particle(i, All.Ti_Current);

                  TimeBinsBh.ActiveParticleList[TimeBinsBh.NActiveParticles] = i;
                  TimeBinsBh.NActiveParticles++;
                }
            }
        }
    }

    mysort(TimeBinsBh.ActiveParticleList, TimeBinsBh.NActiveParticles, sizeof(int), int_compare);

  n = 1;
  int in;
  long long out;

  in = TimeBinsBh.NActiveParticles;

  sumup_large_ints(n, &in, &out);

  TimeBinsBh.GlobalNActiveParticles = out;
}

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}