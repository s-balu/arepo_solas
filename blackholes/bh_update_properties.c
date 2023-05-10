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
/*calculate bondi accretion rate*/
  int i;
  double density, pressure, sound_speed, velocity_gas_norm;
  double denominator, denominator_inv, BondiRate, EddingtonRate;
  double accretion_rate, acc_rate_for_print;

  accretion_rate = acc_rate_for_print = 0;

  for(i = 0; i < NumBh; i++)
    {
/*get pressure*/
      if(BhP[i].Density>0)
        {  
          density = BhP[i].Density;
          pressure = GAMMA_MINUS1 * density * BhP[i].InternalEnergyGas;

/*get soundspeed*/
          sound_speed = sqrt(GAMMA * pressure / density);
      
          velocity_gas_norm = sqrt(BhP[i].VelocityGas[0]*BhP[i].VelocityGas[0] + 
          BhP[i].VelocityGas[1]*BhP[i].VelocityGas[1] + BhP[i].VelocityGas[2]*BhP[i].VelocityGas[2]); 

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
      EddingtonRate = 4. * M_PI * GRAVITY * (PPB(i).Mass * All.UnitMass_in_g) * PROTONMASS / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *=  (All.UnitTime_in_s / All.UnitMass_in_g);
      accretion_rate = fmin(BondiRate, EddingtonRate);
  
/*efficiency*/
      accretion_rate *= (1. - All.Epsilon_r);
      BhP[i].AccretionRate = accretion_rate;
    }
  
  MPI_Allreduce(&accretion_rate, &acc_rate_for_print, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
  mpi_printf("BLACK_HOLES: Black hole accretion rate: %e \n", acc_rate_for_print);
}

/*get timestep for bh based on smallest between ngbmin and acc_timestep*/
integertime get_timestep_bh(int p)
{
  MyDouble accretion_timestep;
  integertime acc_timestep;

      accretion_timestep = BhP[p].NgbMass / BhP[p].AccretionRate;
      
      if(acc_timestep > 1) /*if accretion rate is a small number the timestep becomes very large*/
        return BhP[p].NgbMinStep;
     
      acc_timestep = accretion_timestep / All.Timebase_interval;
      acc_timestep *= 0.01;
      
      if(BhP[p].NgbMinStep < acc_timestep)
        return BhP[p].NgbMinStep;
      
      return acc_timestep;
}

/*update bh-timestep at prior_mesh_construction*/
void update_bh_timesteps(void)
{
  int idx, i, binold, bin;
  integertime ti_step;

  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      ti_step = get_timestep_bh(P[i].BhID);
      binold = P[i].TimeBinBh;
      bin = get_timestep_bin(ti_step);

      timebin_move_particle(&TimeBinsBh, i, binold, bin);
      P[i].TimeBinBh = bin;
    }
}

/*call this function as the reconstruct_timebins() bh version*/
void reconstruct_bh_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsBh.TimeBinCount[bin]   = 0;
      TimeBinsBh.FirstInTimeBin[bin] = -1;
      TimeBinsBh.LastInTimeBin[bin]  = -1;
    }
  
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type != 5)
        continue;
      bin = P[i].TimeBinBh;

      if(TimeBinsBh.TimeBinCount[bin] > 0)
        {
          TimeBinsBh.PrevInTimeBin[i]                                  = TimeBinsBh.LastInTimeBin[bin];
          TimeBinsBh.NextInTimeBin[i]                                  = -1;
          TimeBinsBh.NextInTimeBin[TimeBinsBh.LastInTimeBin[bin]]      = i;
          TimeBinsBh.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsBh.FirstInTimeBin[bin] = TimeBinsBh.LastInTimeBin[bin] = i;
          TimeBinsBh.PrevInTimeBin[i] = TimeBinsBh.NextInTimeBin[i] = -1;
        }
      TimeBinsBh.TimeBinCount[bin]++;
    }
  update_list_of_active_bh_particles();
}

/*call this function after updating the bh-timebin to the ngb condition*/
void update_list_of_active_bh_particles(void)
{
  int i, n;
  TimeBinsBh.NActiveParticles = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      //if(TimeBinSynchronized[n]) --> need early bh timestep criteria in the run loop to include this
        //{
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
        //}
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
/*accrete mass, angular momentum onto the bh and drain ngb cells*/
  int i, j, bin;
  double dt;
  double cos_theta, p0, pj;
  
  for(i=0; i<NumBh; i++)
    {
      bin = PPB(i).TimeBinBh;
      dt    = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
      
      PPB(i).Mass += (1-All.Epsilon_f) * BhP[i].AccretionRate * dt;

      BhP[i].AngularMomentum[0] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[0];
      BhP[i].AngularMomentum[1] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[1];
      BhP[i].AngularMomentum[2] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[2];
    
      for(j=0; j<NumGas; j++)
        {
          if(SphP[j].MassDrain > 0)
          {
            if(P[j].Mass - SphP[j].MassDrain < 0.1*P[j].Mass)
              {
                P[j].Mass -= 0.9*P[j].Mass;
                BhP[i].MassToDrain += SphP[j].MassDrain - 0.9*P[j].Mass; 
                /*we're also losing thermal and kinetic energy & momentum*/
                
                /*update total energy*/
                SphP[j].Energy *= 0.1;

                /*update momentum*/
                SphP[j].Momentum[0] *= 0.1;
                SphP[j].Momentum[1] *= 0.1;
                SphP[j].Momentum[2] *= 0.1;
              }
            else
              {
                P[j].Mass -= SphP[j].MassDrain;
                
                /*update total energy*/
                SphP[j].Energy *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);

                /*update momentum*/
                SphP[j].Momentum[0] *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
                SphP[j].Momentum[1] *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
                SphP[j].Momentum[2] *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
              }
          
            SphP[j].MassDrain = 0;
          }
        }
    }

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

          /*momentum kick direction along jet axis*/
          double pos_x_axis[3] = {1, 0, 0};
          double neg_x_axis[3] = {-1, 0, 0};

          for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              int i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0)
              continue;
              if(SphP[i].ThermalFeed > 0 || SphP[i].KineticFeed > 0)
                {
                  /*calculate momentum feed exactly so energy is conserved*/
                  /*-> we need to do this here so that particle properties don't change between loading the buffer and writing it*/
                  p0 = sqrt(pow(SphP[i].Momentum[0], 2) + pow(SphP[i].Momentum[1], 2) + pow(SphP[i].Momentum[2], 2));
              
                  if(p0 < pow(10,-10)) //protect against p0 = 0;
                    cos_theta = 1;
                  else if(SphP[i].PositiveJet)
                    cos_theta = (SphP[i].Momentum[0]*pos_x_axis[0] + SphP[i].Momentum[1]*pos_x_axis[1] + SphP[i].Momentum[2]*pos_x_axis[2]) / 
                    (p0*sqrt(pow(pos_x_axis[0], 2) + pow(pos_x_axis[1], 2) + pow(pos_x_axis[2], 2)));       
                  else
                    cos_theta = (SphP[i].Momentum[0]*neg_x_axis[0] + SphP[i].Momentum[1]*neg_x_axis[1] + SphP[i].Momentum[2]*neg_x_axis[2]) / 
                    (p0*sqrt(pow(neg_x_axis[0], 2) + pow(neg_x_axis[1], 2) + pow(neg_x_axis[2], 2)));
          
                  pj = -p0*cos_theta + sqrt(p0*p0 * cos_theta*cos_theta + 2*P[i].Mass*SphP[i].KineticFeed);

                  if(SphP[i].PositiveJet) 
                    { 
                      SphP[i].MomentumFeed[0] += pos_x_axis[0] * pj / sqrt(pow(pos_x_axis[0], 2) + pow(pos_x_axis[1], 2) +  pow(pos_x_axis[2], 2));
                      SphP[i].MomentumFeed[1] += pos_x_axis[1] * pj / sqrt(pow(pos_x_axis[0], 2) + pow(pos_x_axis[1], 2) +  pow(pos_x_axis[2], 2));
                      SphP[i].MomentumFeed[2] += pos_x_axis[2] * pj / sqrt(pow(pos_x_axis[0], 2) + pow(pos_x_axis[1], 2) +  pow(pos_x_axis[2], 2)); 
                    }
                  else
                    {  
                      SphP[i].MomentumFeed[0] += neg_x_axis[0] * pj / sqrt(pow(neg_x_axis[0], 2) + pow(neg_x_axis[1], 2) +  pow(neg_x_axis[2], 2));
                      SphP[i].MomentumFeed[1] += neg_x_axis[1] * pj / sqrt(pow(neg_x_axis[0], 2) + pow(neg_x_axis[1], 2) +  pow(neg_x_axis[2], 2));
                      SphP[i].MomentumFeed[2] += neg_x_axis[2] * pj / sqrt(pow(neg_x_axis[0], 2) + pow(neg_x_axis[1], 2) +  pow(neg_x_axis[2], 2));
                    }

                  /*update total energy*/
                  SphP[i].Energy += SphP[i].ThermalFeed + SphP[i].KineticFeed;
                  All.EnergyExchange[1] += SphP[i].ThermalFeed + SphP[i].KineticFeed;
                  /*update momentum*/
                    SphP[i].Momentum[0] += SphP[i].MomentumFeed[0];
                    SphP[i].Momentum[1] += SphP[i].MomentumFeed[1];
                    SphP[i].Momentum[2] += SphP[i].MomentumFeed[2];
                  /*update velocities*/
                  update_primitive_variables_single(P, SphP, i, &pvd);
                  /*update internal energy*/
                  update_internal_energy(P, SphP, i, &pvd);
                  /*update pressure*/
                  set_pressure_of_cell_internal(P, SphP, i);
                  /*set feed flags to zero*/
                  SphP[i].ThermalFeed = SphP[i].KineticFeed = 0;
                  SphP[i].MomentumFeed[0] = SphP[i].MomentumFeed[1] = SphP[i].MomentumFeed[2] = 0;
                }
            }
#ifdef BURST_MODE
          All.FeedbackFlag = -1;
#endif      
        }
    }
  MPI_Allreduce(&All.EnergyExchange, &All.EnergyExchangeTot, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
  mpi_printf("BLACK_HOLES: Energy given by BH = %e, Energy taken up by gas particles = %e \n", All.EnergyExchangeTot[0], All.EnergyExchangeTot[1]);

#ifdef BURST_MODE
  if(All.EnergyExchangeTot[0] - All.EnergyExchangeTot[1] > 10)  
    All.FeedbackFlag = 1;
#endif   
}

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}