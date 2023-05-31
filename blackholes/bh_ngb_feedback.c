#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"


static int bh_ngb_feedback_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Bin;
  int IsBh;
  MyDouble BhRho;
  MyDouble NgbMass;
  MyDouble NgbMassFeed;
  MyDouble AccretionRate;
  MyDouble MassToDrain;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0]        = PPB(i).Pos[0];
  in->Pos[1]        = PPB(i).Pos[1];
  in->Pos[2]        = PPB(i).Pos[2];
  in->Hsml          = BhP[i].Hsml;
  in->Bin           = BhP[i].TimeBinBh;
  in->IsBh          = BhP[i].IsBh;
  in->BhRho         = BhP[i].Density;
  in->NgbMass       = BhP[i].NgbMass;
  in->NgbMassFeed   = BhP[i].NgbMassFeed;
  in->AccretionRate = BhP[i].AccretionRate;
  in->MassToDrain   = BhP[i].MassToDrain;
  in->Firstnode     = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.*/

typedef struct
{

} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void*/
 
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      return;
    }
  else /* combine */
    {
      return;
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int idx;

  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        idx = NextParticle++;

        if(idx >= TimeBinsBh.NActiveParticles)
          break;

        int i = TimeBinsBh.ActiveParticleList[idx];
        if(i < 0)
          continue;
        bh_ngb_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        bh_ngb_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void bh_ngb_feedback(void)
{
  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);
}

static int bh_ngb_feedback_evaluate(int target, int mode, int threadid)
{
  int j, n, bin, isbh;
  int numnodes, *firstnode;
  double h, h2, hinv, hinv3, hinv4, wk, dwk;
  double dx, dy, dz, r, r2, u;
  double dt; //dtime;
  MyDouble ngbmass, ngbmass_feed, accretion_rate, mass_to_drain; 
  MyDouble *pos, bh_rho, energy_feed;

  data_in local, *target_data;
  /*data_out out;*/

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos            = target_data->Pos;
  h              = target_data->Hsml;
  isbh           = target_data->IsBh;
  bin            = target_data->Bin;
  bh_rho         = target_data->BhRho;
  ngbmass        = target_data->NgbMass;
  ngbmass_feed   = target_data->NgbMassFeed;
  accretion_rate = target_data->AccretionRate;
  mass_to_drain  = target_data->MassToDrain; 

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;
 
/*bh timestep*/
  dt    = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
  //dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  energyfeed = All.Epsilon_f * All.Epsilon_r * accretion_rate * dt * (CLIGHT * CLIGHT / (All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s));

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);
  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

#ifndef REFLECTIVE_X
      if(dx > boxHalf_X)
        dx -= boxSize_X;
      if(dx < -boxHalf_X)
        dx += boxSize_X;
#endif /* #ifndef REFLECTIVE_X */

#ifndef REFLECTIVE_Y
      if(dy > boxHalf_Y)
        dy -= boxSize_Y;
      if(dy < -boxHalf_Y)
        dy += boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y */

#ifndef REFLECTIVE_Z
      if(dz > boxHalf_Z)
        dz -= boxSize_Z;
      if(dz < -boxHalf_Z)
        dz += boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z */
      r2 = dx * dx + dy * dy + dz * dz;

/*kernel*/
      if(r2 < h2)
        {
          r = sqrt(r2);

          u = r * hinv;

          kernel(u, hinv3, hinv4, &wk, &dwk);

          if(isbh)
            {

            }
          
          if(!isbh)
            {
/*set radial momentum kick*/
              SphP[j].MomentumFeed  += energy_feed / (CLIGHT / All.UnitVelocity_in_cm_per_s) * P[j].Mass / bh_rho * wk;
              All.EnergyExchange[2] += energy_feed / (CLIGHT / All.UnitVelocity_in_cm_per_s) * P[j].Mass / bh_rho * wk;

              SphP[j].MomentumKickVector[0] = -dx;
              SphP[j].MomentumKickVector[1] = -dy;
              SphP[j].MomentumKickVector[2] = -dz;
            }
        }
    }
  return 0;
}
  /* Now collect the result at the right place 
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}*/

