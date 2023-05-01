#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

static int int_compare(const void *a, const void *b);

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
              BondiRate = 4. * M_PI * All.G * All.G * PPB(i).Mass * density *
                   denominator_inv * denominator_inv * denominator_inv;
            }
          else
            terminate("Invalid denominator in Bondi Accretion Rate");
        }
      else
        BondiRate = 0;
  
/*limit by Eddington accretion rate*/
      EddingtonRate = 4. * M_PI * GRAVITY * (PPB(i).Mass / All.UnitMass_in_g )* PROTONMASS / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *=  (All.UnitMass_in_g / All.UnitTime_in_s);
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
      bin = get_timestep_bin(ti_step);

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

void perform_end_of_step_bh_physics(void)
{
/*inject feedback to ngb cells*/
    if(All.Time >= All.FeedbackTime)
    {   
      if(All.FeedbackFlag > 0)
        {
          struct pv_update_data pvd;
          if(All.ComovingIntegrationOn)
            {
              pvd.atime    = All.Time;
              pvd.hubble_a = hubble_function(All.Time);
              pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
            }
          else
            pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

          for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              int i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0)
              continue;
              if(SphP[i].ThermalFeed > 0 || SphP[i].KineticFeed > 0)
                {
/*update total energy*/
                  SphP[i].Energy += SphP[i].ThermalFeed + SphP[i].KineticFeed;
                  All.EnergyExchange[1] += SphP[i].ThermalFeed + SphP[i].KineticFeed;
/*update momentum -> include flag for momentum direction*/ 
                  if(P[i].Pos[0] >= 0.5)
                    SphP[i].Momentum[0] += sqrt(2 * P[i].Mass * SphP[i].KineticFeed);
                  if(P[i].Pos[0] < 0.5)
                    SphP[i].Momentum[0] -= sqrt(2 * P[i].Mass * SphP[i].KineticFeed);
/*update velocities*/
                  update_primitive_variables_single(P, SphP, i, &pvd);
/*update internal energy*/
                  update_internal_energy(P, SphP, i, &pvd);
/*update pressure*/
                  set_pressure_of_cell_internal(P, SphP, i);
/*set feed flags to zero*/
                  SphP[i].ThermalFeed = SphP[i].KineticFeed = 0;
                }
            }
#ifdef BURST_MODE
          All.FeedbackFlag = -1;
#endif      
        }
    }
  MPI_Reduce(&All.EnergyExchange, &All.EnergyExchangeTot, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
  mpi_printf("Energy given by BH = %f, Energy taken up by gas particles = %f \n", All.EnergyExchangeTot[0], All.EnergyExchangeTot[1]);

#ifdef BURST_MODE
  if(All.EnergyExchangeTot[0] - All.EnergyExchangeTot[1] > 10)  
    All.FeedbackFlag = 1;
#endif   

/*accrete mass, angular momentum into the bh and drain ngb cells*/
  int i, j, bin;
  double dt;
  
  for(i=0; i<NumBh; i++)
    {
      bin = PPB(i).TimeBinBh;
      dt    = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
      
      PPB(i).Mass += BhP[i].AccretionRate * dt;

      BhP[i].AngularMomentum[0] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[0];
      BhP[i].AngularMomentum[1] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[1];
      BhP[i].AngularMomentum[2] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[2];

      for(j=0; j<NumGas; j++)
        P[j].Mass -= SphP[j].MassDrain;
    }
}

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}