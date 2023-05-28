#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"


static int bh_density_evaluate(int target, int mode, int threadid);
static int bh_density_isactive(int n);

static MyFloat *BhNumNgb; //*BhDhsmlDensityFactor;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;
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
  in->Vel[0]        = PPB(i).Vel[0];
  in->Vel[1]        = PPB(i).Vel[1];
  in->Vel[2]        = PPB(i).Vel[2];
  in->Hsml          = BhP[i].Hsml;
  in->Firstnode     = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble DhsmlDensity;
  MyDouble Ngb;
  MyDouble Rho;
  MyDouble Mass;
  integertime NgbMinStep;
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
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      /*BhDhsmlDensityFactor[i]          = out->DhsmlDensity;*/
      BhNumNgb[i]                      = out->Ngb;
      BhP[i].Density                   = out->Rho;
      BhP[i].NgbMass                   = out->Mass;
      BhP[i].NgbMinStep                = out->NgbMinStep;
    }
  else /* combine */
    {
      /*BhDhsmlDensityFactor[i]          += out->DhsmlDensity;*/
      BhNumNgb[i]                      += out->Ngb;
      BhP[i].Density                   += out->Rho;
      BhP[i].NgbMass                   += out->Mass;
      if(out->NgbMinStep < BhP[i].NgbMinStep)
        BhP[i].NgbMinStep               = out->NgbMinStep;
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
        if(bh_density_isactive(P[i].BhID))
          bh_density_evaluate(P[i].BhID, MODE_LOCAL_PARTICLES, threadid);
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

        bh_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Main function of SPH density calculation.
 *
 *  This function computes the local density for each active SPH particle and
 *  the number of weighted neighbors in the current smoothing radius. If a
 *  particle with its smoothing region is fully inside the local domain, it is
 *  not exported to the other processors. The function also detects particles
 *  that have a number of neighbors outside the allowed tolerance range. For
 *  these particles, the smoothing length is adjusted accordingly, and the
 *  computation is called again.
 *
 *  \return void
 */
void bh_density(void)
{
  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;
  double desnumngb, t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  BhNumNgb             = (MyFloat *)mymalloc("BhNumNgb", NumBh * sizeof(MyFloat));
  /*BhDhsmlDensityFactor = (MyFloat *)mymalloc("BhDhsmlDensityFactor", NumBh * sizeof(MyFloat));*/
  Left               = (MyFloat *)mymalloc("Left", NumBh * sizeof(MyFloat));
  Right              = (MyFloat *)mymalloc("Right", NumBh * sizeof(MyFloat));

  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(bh_density_isactive(P[i].BhID))
        {
          Left[P[i].BhID] = Right[P[i].BhID] = 0;
        }
    }

  generic_set_MaxNexport();

  desnumngb = All.BhDesNumNgb;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

      for(idx=0, npleft=0; idx<TimeBinsBh.NActiveParticles; idx++)
        {
          i = TimeBinsBh.ActiveParticleList[idx];
          if(i < 0)
            continue;
          /*if(bh_density_isactive(P[i].BhID))
            {
              if(BPP(i).Density > 0)
                {
                  BhDhsmlDensityFactor[P[i].BhID] *= BPP(i).Hsml / (NUMDIMS * BPP(i).Density);
                  if(BhDhsmlDensityFactor[P[i].BhID] > -0.9)*/ /* note: this would be -1 if only a single particle at zero lag is found */ /*
                      BhDhsmlDensityFactor[P[i].BhID] = 1 / (1 + BhDhsmlDensityFactor[P[i].BhID]);
                  else
                      BhDhsmlDensityFactor[P[i].BhID] = 1;
                }
            }
          */

          if(BhNumNgb[P[i].BhID] < (desnumngb - All.BhDesNumNgbDev) || BhNumNgb[P[i].BhID] > (desnumngb + All.BhDesNumNgbDev))
          {
                  /* need to redo this particle */
            npleft++;

            if(Left[P[i].BhID] > 0 && Right[P[i].BhID] > 0)
              {
                if((Right[P[i].BhID] - Left[P[i].BhID]) < 1.0e-3 * Left[P[i].BhID])
                  {
                        /* this one should be ok */
                    npleft--;
                    P[i].TimeBinBh = -P[i].TimeBinBh - 1; /* Mark as inactive */
                    continue;
                }
              } 

            if(BPP(i).NgbMass < (desnumngb - All.BhDesNumNgbDev))
              Left[P[i].BhID] = dmax(BPP(i).Hsml, Left[P[i].BhID]);
            else
              {
                if(Right[P[i].BhID] != 0)
                  {
                    if(BPP(i).Hsml < Right[P[i].BhID])
                        Right[P[i].BhID] = BPP(i).Hsml;
                  }
                    else
                        Right[P[i].BhID] = BPP(i).Hsml;
              }

            if(iter >= MAXITER - 10)
              {
                printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g NgbsMass=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i, ThisTask,
                    (int)P[i].ID, BPP(i).Hsml, Left[P[i].BhID], Right[P[i].BhID], (float)BPP(i).NgbMass, Right[P[i].BhID] - Left[P[i].BhID], P[i].Pos[0],
                    P[i].Pos[1], P[i].Pos[2]);
                myflush(stdout);
              }

            if(Right[P[i].BhID] > 0 && Left[P[i].BhID] > 0)
                BPP(i).Hsml = pow(0.5 * (pow(Left[P[i].BhID], 3) + pow(Right[P[i].BhID], 3)), 1.0 / 3);
            else
              {
                if(Right[P[i].BhID] == 0 && Left[P[i].BhID] == 0)
                    terminate("should not occur");

                if(Right[P[i].BhID] == 0 && Left[P[i].BhID] > 0)
                  {
                    BPP(i).Hsml *= 1.26;
                  }

                if(Right[P[i].BhID] > 0 && Left[P[i].BhID] == 0)
                  {
                    BPP(i).Hsml /= 1.26;
                  }
              }
          }
        else
             P[i].TimeBinBh = -P[i].TimeBinBh - 1; /* Mark as inactive */ 
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("BH_DENSITY: ngb iteration %3d: need to repeat for %12lld particles. NgbsMass=%g. (took %g sec)\n", iter, ntot, (float)BPP(i).NgbMass,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in density()\n");
        }
    }
  while(ntot > 0);


  myfree(Right);
  myfree(Left);
  /*myfree(BhDhsmlDensityFactor);*/
  myfree(BhNumNgb);

  /* mark as active again */
for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinBh < 0)
        P[i].TimeBinBh = -P[i].TimeBinBh - 1;
    }
  /* collect some timing information */
  CPU_Step[CPU_INIT] += measure_time();
}

/*! \brief Inner function of the SPH density calculation
 *
 *  This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 *
 *  \param[in] target Index of particle in local data/import buffer.
 *  \param[in] mode Mode in which function is called (local or impored data).
 *  \param[in] threadid ID of local thread.
 *
 *  \return 0
 */
static int bh_density_evaluate(int target, int mode, int threadid)
{
  int j, n;
  int numngb, numnodes, *firstnode;
  double h, h2, hinv, hinv3, hinv4;
  double rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat weighted_numngb;
  /*MyFloat dhsmlrho;*/
  MyDouble *pos;
  MyDouble mass;
  integertime ngb_min_step;
  int bin = TIMEBINS;
   
  data_in local, *target_data;
  data_out out;

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

  pos = target_data->Pos;
  h   = target_data->Hsml;

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;

  numngb = 0;
  rho = weighted_numngb =/*= dhsmlrho = */0;
  mass = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

/*compute the min hydro step for neighbors*/     
      if(bin > P[j].TimeBinHydro)
        bin = P[j].TimeBinHydro;

/*compute the bh-ngb-mass*/
      mass += P[j].Mass;

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

/* now find the closest image in the given box size */
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
          numngb++;

          r = sqrt(r2);

          u = r * hinv;
          
          kernel(u, hinv3, hinv4, &wk, &dwk);

          mass_j = P[j].Mass;

/*compute bh density*/
          rho += FLT(mass_j * wk);

          weighted_numngb += FLT(NORM_COEFF * wk / hinv3); /* 4.0/3 * PI = 4.188790204786 */

          //dhsmlrho += FLT(-mass_j * (NUMDIMS * hinv * wk + u * dwk));
        }
    }

/*compute bh timestep based on min ngb timestep*/
  if(bin == 0)
    ngb_min_step = 0;
  else
    ngb_min_step   = (((integertime)1) << bin);
  
  //out.DhsmlDensity            = dhsmlrho;
  out.Ngb                     = weighted_numngb;
  out.Rho                     = rho;
  out.Mass                    = mass;
  out.NgbMinStep              = ngb_min_step;


  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/* \brief Determines if a BhP is active in current timestep.
 *
 *  If the BhP is not active in a timestep, its value in TimeBinBh is
 *  negative.
 *
 *  \param[in] n Index of BhP in Particle array
 *
 *  \return 1: cell active; 0: BhP not active.
 */
int bh_density_isactive(int n)
{
  if(PPB(n).TimeBinBh < 0)
    return 0;

  return 1;
}
