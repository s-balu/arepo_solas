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
 * \file        src/init/density.c
 * \date        05/2018
 * \brief       SPH density computation and smoothing length determination.
 * \details     This file contains the "first SPH loop", where the SPH
 *              densities and smoothing lengths are calculated.
 *              In Arepo, this is used in setup_smoothinglengths() (init.c) to
 *              get an initial guess for MaxDelaunayRadius.
 *              Note that the SPH density is NOT used in the subsequent
 *              hydrodynamics calculation, but the density is either set by the
 *              initial conditions explicitly (DENSITY_AS_MASS_IN_INPUT) or
 *              calculated by the mass given in the initial conditions divided
 *              by the volume of the cell calculated by the Voronoi
 *              tessellation algorithm.
 *              contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void density(void)
 *                static int density_evaluate(int target, int mode, int
 *                  threadid)
 *                int density_isactive(int n)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

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

static MyFloat *BhNumNgb, *BhDhsmlDensityFactor;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
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
  in->Pos[0] = P[BhP[i].PID].Pos[0];
  in->Pos[1] = P[BhP[i].PID].Pos[1];
  in->Pos[2] = P[BhP[i].PID].Pos[2];
  in->Hsml   = BhP[i].Hsml;

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyFloat Rho;
  MyFloat DhsmlDensity;
  MyFloat Ngb;
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
      BhNumNgb[i] = out->Ngb;
/*      if(P[i].Type == 5)
        {
*/
      BhP[i].Density       = out->Rho;
      BhDhsmlDensityFactor[i] = out->DhsmlDensity;
/*        }
*/
    }
  else /* combine */
    {
      BhNumNgb[i] += out->Ngb;
/*      if(P[i].Type == 5)
        {
*/
      BhP[i].Density += out->Rho;
      BhDhsmlDensityFactor[i] += out->DhsmlDensity;
/*        }
*/
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

       /*if(idx >= TimeBinsHydro.NActiveParticles)
         break;

        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(density_isactive(i))*/
        
        if(idx >= NumBh)
          break;
        if(bh_density_isactive(idx))
          bh_density_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

static MyFloat *BhNumNgb, *BhDhsmlDensityFactor;

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
  int idx, npleft, iter = 0;
  long long ntot;
  double desnumngb, t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  BhNumNgb             = (MyFloat *)mymalloc("BhNumNgb", NumBh * sizeof(MyFloat));
  BhDhsmlDensityFactor = (MyFloat *)mymalloc("BhDhsmlDensityFactor", NumBh * sizeof(MyFloat));
  Left               = (MyFloat *)mymalloc("Left", NumBh * sizeof(MyFloat));
  Right              = (MyFloat *)mymalloc("Right", NumBh * sizeof(MyFloat));

/*  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(density_isactive(i))
        {
          Left[i] = Right[i] = 0;
        }
    }
*/

  for(idx=0; idx<NumBh; idx++)
  {
      BhP[idx].mark=1;
  }
  for(idx=0; idx<NumBh; idx++)
  {
      if(bh_density_isactive(idx))
      {
          Left[idx] = Right[idx] = 0;
      } 
  }

  generic_set_MaxNexport();

  desnumngb = All.DesNumNgb;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(NumBh, kernel_local, kernel_imported);

      /* do final operations on results */
/*      for(idx = 0, npleft = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(density_isactive(i))
            {
              if(P[i].Type == 5)
                {
                  if(BhP[i].Density > 0)
                    {
                      DhsmlDensityFactor[i] *= BhP[i].Hsml / (NUMDIMS * BhP[i].Density);
                      if(DhsmlDensityFactor[i] > -0.9)  note: this would be -1 if only a single particle at zero lag is found 
                        DhsmlDensityFactor[i] = 1 / (1 + DhsmlDensityFactor[i]);
                      else
                        DhsmlDensityFactor[i] = 1;
                    }
                }
*/

      for(idx=0, npleft=0; idx<NumBh; idx++)
        {
          if(bh_density_isactive(idx))
            {
          
              if(BhP[idx].Density > 0)
                {
                  BhDhsmlDensityFactor[idx] *= BhP[idx].Hsml / (NUMDIMS * BhP[idx].Density);
                  if(BhDhsmlDensityFactor[idx] > -0.9) /* note: this would be -1 if only a single particle at zero lag is found */
                      BhDhsmlDensityFactor[idx] = 1 / (1 + BhDhsmlDensityFactor[idx]);
                  else
                      BhDhsmlDensityFactor[idx] = 1;
                }
            } 
        

          if(BhNumNgb[idx] < (desnumngb - All.MaxNumNgbDeviation) || BhNumNgb[idx] > (desnumngb + All.MaxNumNgbDeviation))
          {
                  /* need to redo this particle */
            npleft++;

            if(Left[idx] > 0 && Right[idx] > 0)
              {
                if((Right[idx] - Left[idx]) < 1.0e-3 * Left[idx])
                  {
                        /* this one should be ok */
                    npleft--;
                    BhP[idx].mark = -1 - BhP[idx].mark; /* Mark as inactive */
                    continue;
                }
              } 

            if(BhNumNgb[idx] < (desnumngb - All.MaxNumNgbDeviation))
              Left[idx] = dmax(BhP[idx].Hsml, Left[idx]);
            else
              {
                if(Right[idx] != 0)
                  {
                    if(BhP[idx].Hsml < Right[idx])
                        Right[idx] = BhP[idx].Hsml;
                  }
                    else
                        Right[idx] = BhP[idx].Hsml;
              }

            if(iter >= MAXITER - 10)
              {
                printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", idx, ThisTask,
                    (int)BhP[idx].PID, BhP[idx].Hsml, Left[idx], Right[idx], (float)BhNumNgb[idx], Right[idx] - Left[idx], P[BhP[idx].PID].Pos[0],
                    P[BhP[idx].PID].Pos[1], P[BhP[idx].PID].Pos[2]);
                myflush(stdout);
              }

            if(Right[idx] > 0 && Left[idx] > 0)
                BhP[idx].Hsml = pow(0.5 * (pow(Left[idx], 3) + pow(Right[idx], 3)), 1.0 / 3);
            else
              {
                if(Right[idx] == 0 && Left[idx] == 0)
                    terminate("should not occur");

                if(Right[idx] == 0 && Left[idx] > 0)
                  {
                    BhP[idx].Hsml *= 1.26;
                  }

                if(Right[idx] > 0 && Left[idx] == 0)
                  {
                    BhP[idx].Hsml /= 1.26;
                  }
              }
          }
        else
            BhP[idx].mark = -1 - BhP[idx].mark; /* Mark as inactive */
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("DENSITY: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in density()\n");
        }
    }
  while(ntot > 0);


  myfree(Right);
  myfree(Left);
  myfree(BhDhsmlDensityFactor);
  myfree(BhNumNgb);

  /* mark as active again */
/*  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinHydro < 0)
        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;
    }
*/
    for(idx=0; idx<NumBh; idx++)
    {
        if(BhP[idx].mark < 0)
        {
            BhP[idx].mark = 1;
        }
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
  MyFloat rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat weighted_numngb;
  MyFloat dhsmlrho;
  MyDouble *pos;

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
  rho = weighted_numngb = dhsmlrho = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

/*  now find the closest image in the given box size  */
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

      if(r2 < h2)
        {
          numngb++;

          r = sqrt(r2);

          u = r * hinv;

          if(u < 0.5)
            {
              wk  = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
            }
          else
            {
              wk  = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
              dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
            }

          mass_j = P[j].Mass;

          rho += FLT(mass_j * wk);

          weighted_numngb += FLT(NORM_COEFF * wk / hinv3); /* 4.0/3 * PI = 4.188790204786 */

          dhsmlrho += FLT(-mass_j * (NUMDIMS * hinv * wk + u * dwk));
        }
    }

  out.Rho          = rho;
  out.Ngb          = weighted_numngb;
  out.DhsmlDensity = dhsmlrho;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/* \brief Determines if a cell is active in current timestep.
 *
 *  If the cell is not active in a timestep, its value in TimeBinHydro is
 *  negative.
 *
 *  \param[in] n Index of cell in P and SphP arrays.
 *
 *  \return 1: cell active; 0: cell not active or not a cell.
 */
int bh_density_isactive(int n)
{
  if(NumBh==0)
    return 0;
  if(BhP[n].mark <= 0)
    return 0;

  return 1;
}
