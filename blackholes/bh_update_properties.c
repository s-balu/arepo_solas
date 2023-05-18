#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

static int int_compare(const void *a, const void *b);

/*sph loop kernel function*/
void kernel(double u, double hinv3, double hinv4, double *wk, double *dwk, int mode)
{
#ifdef CUBIC_SPLINE_KERNEL
#if defined(WENDLAND_C2_KERNEL) || defined(WENDLAND_C4_KERNEL) || defined(WENDLAND_C6_KERNEL)
#error "Only one SPH kernel can be used"
#endif
  if(u < 0.5)
    {
      if(mode >= COMPUTE_WK_AND_DWK)
        *dwk = u * (18.0 * u - 12.0);
      if(mode <= COMPUTE_WK_AND_DWK)
        *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      if(mode >= COMPUTE_WK_AND_DWK)
        *dwk = -6.0 * t2;
      if(mode <= COMPUTE_WK_AND_DWK)
        *wk = 2.0 * t2 * t1;
    }
#endif

#ifdef WENDLAND_C2_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);

  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -12.0 * u * t2;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t2 * t1 * (1.0 + u * 3.0);

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -20.0 * u * t2 * t1;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t4 * (1.0 + u * 4.0);

#endif
#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;
  double t5 = t4 * t1;

  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -14.0 * t4 * (4.0 * u + 1) * u;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t5 * (1.0 + u * (5.0 + 8.0 * u));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t2 * t2 * t2;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -56.0 / 3.0 * u * t4 * t1 * (5.0 * u + 1);
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t6 * (1.0 + u * (6.0 + 35.0 / 3.0 * u));

#endif
#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t4 * t2;
  double t7 = t4 * t2 * t1;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -6.0 * u * t6 * (3.0 + u * (18.0 + 35.0 * u));
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t7 * (1.0 + u * (7.0 + u * (19.0 + 21.0 * u)));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t7 = t4 * t2 * t1;
  double t8 = t4 * t4;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -22.0 * u * (1.0 + u * (7.0 + 16.0 * u)) * t7;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t8 * (1.0 + u * (8.0 + u * (25.0 + 32.0 * u)));

#endif
#endif /* WENDLAND_C6_KERNEL */
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk *= NORM * hinv4;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk *= NORM * hinv3;
}

void update_bh_accretion_rate(void)
{
  int i;
 
  double EddingtonRate;
 /* double acc_rate_for_print;

  acc_rate_for_print = 0;*/

  for(i = 0; i < NumBh; i++)
    {
      EddingtonRate = 4. * M_PI * GRAVITY * (PPB(i).Mass * All.UnitMass_in_g) * PROTONMASS / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *=  (All.UnitTime_in_s / All.UnitMass_in_g);

      BhP[i].AccretionRate  = EddingtonRate; //fix accretion to eddington rate
    }
  
  /*MPI_Allreduce(&BhP[i].AccretionRate, &acc_rate_for_print, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
  mpi_printf("BLACK_HOLES: Black hole accretion rate: %e \n", acc_rate_for_print);*/
}

/*get timestep for bh based on smallest between ngbs*/
integertime get_timestep_bh(int p)
{
  return BhP[p].NgbMinStep;
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

/*call this function after reconstruct_bh_timebins*/
void update_list_of_active_bh_particles(void)
{
  int i, n;
  TimeBinsBh.NActiveParticles = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      //if(TimeBinSynchronized[n]) 
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

/*call this function after updating the bh-timebin to the ngb condition*/
void update_list_of_active_bh_particles_prior_mesh(void)
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
/*accrete mass, angular momentum onto the bh and drain ngb cells*/
  int idx, i;
  double pj;

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

          /*radial momentum kick*/
          double kick_vector[3];

          for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0)
              continue;
              if(SphP[i].KineticFeed > 0)
                {
                  kick_vector[0] = SphP[i].MomentumKickVector[0];
                  kick_vector[1] = SphP[i].MomentumKickVector[1];
                  kick_vector[2] = SphP[i].MomentumKickVector[2];
                  pj = SphP[i].KineticFeed;

                  /*update momentum*/
                  SphP[i].Momentum[0] += kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                  SphP[i].Momentum[1] += kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                  SphP[i].Momentum[2] += kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));  

                  All.EnergyExchange[1] += SphP[i].KineticFeed;     
                 
                  /*update velocities*/
                  update_primitive_variables_single(P, SphP, i, &pvd);  

                  /*update total energy*/
                  SphP[i].Energy = SphP[i].Utherm * P[i].Mass + 0.5 * P[i].Mass * (pow(P[i].Vel[0], 2) + pow(P[i].Vel[1], 2) + pow(P[i].Vel[2], 2));                 
                  /*update internal energy*/
                  update_internal_energy(P, SphP, i, &pvd);
                  /*update pressure*/
                  set_pressure_of_cell_internal(P, SphP, i);
                  /*set feed flag to zero*/
                  SphP[i].KineticFeed = 0;
                }
            }
#ifdef BURST_MODE
          All.FeedbackFlag = -1;
#endif      
        }
    }
  MPI_Allreduce(&All.EnergyExchange, &All.EnergyExchangeTot, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
  mpi_printf("BLACK_HOLES: Momentum given by StarPart = %e, Momentum taken up by gas particles = %e \n", All.EnergyExchangeTot[0], All.EnergyExchangeTot[1]);

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