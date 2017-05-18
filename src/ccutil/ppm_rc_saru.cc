    //!-------------------------------------------------------------------------
    //! Copyright (c) 2016 MOSAIC Group (MPI-CBG Dresden)
    //!
    //!
    //! This file is part of the PPM_RC program.
    //!
    //! PPM_RC is free software: you can redistribute it and/or modify
    //! it under the terms of the GNU Lesser General Public License
    //! as published by the Free Software Foundation, either
    //! version 3 of the License, or (at your option) any later
    //! version.
    //!
    //! PPM_RC is distributed in the hope that it will be useful,
    //! but WITHOUT ANY WARRANTY; without even the implied warranty of
    //! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    //! GNU General Public License for more details.
    //!
    //! You should have received a copy of the GNU General Public License
    //! and the GNU Lesser General Public License along with PPM_RC. If not,
    //! see <http://www.gnu.org/licenses/>.
    //!-------------------------------------------------------------------------
    //!  MOSAIC Group
    //!  Max Planck Institute of Molecular Cell Biology and Genetics
    //!  Pfotenhauerstr. 108, 01307 Dresden, Germany
    //!
    //!  Author           - y.afshar           Nov   2015
    //!-------------------------------------------------------------------------
#ifndef PPM_RC_SARUPRNG_H
#define PPM_RC_SARUPRNG_H

#include "saruprng.h"

// 2^31-2
static const unsigned int TWO_31=2147483646;

Saru saru1;
Saru saru2;

extern "C" {
  int SaruInitialize(unsigned int *iseed)
  {
    unsigned int seed=*iseed;

    /* creating a new independent stream, seeded using
    current generator's state and input seed.*/

    unsigned int state, wstate(12345678);
    const unsigned int churned1=0xDEADBEEF^(0x1fc4ce47*(seed^(seed>>13)));
    const unsigned int churned2=0x1234567+(0x82948463*(churned1^(churned1>>20)));
    const unsigned int churned3=0x87654321^(0x87655677*(churned2^(churned2>>16)));

    state=churned2+0x12345678+(churned3^wstate);

    unsigned int add=(state+churned1)>>1;

    if (wstate-0x8009d14b<0xda879add-add)
      wstate+=add;
    else
      wstate+=add-0xda879add;

    saru1.setstate(state,wstate);

    return 0;
  }

  /* Mixing two seeds to assure of a unique seed. */
  unsigned int Saru_SEEDPRNG(unsigned int *iseed1, unsigned int *iseed2) {
    unsigned int seed1=*iseed1;
    unsigned int seed2=*iseed2;

    seed2+=seed1<<16;
    seed1+=seed2<<11;
    seed2+=((signed int)seed1)>>7;
    seed1^=((signed int)seed2)>>3;
    seed2*=0xA5366B4D;
    seed2^=seed2>>10;
    seed2^=((signed int)seed2)>>19;
    seed1+=seed2^0x6d2d4e11;

    seed1 = 0x79dedea3*(seed1^(((signed int)seed1)>>14));
    seed2 = (seed1 + seed2) ^ (((signed int)seed1)>>8);
    seed1 = seed1 + (seed2*(seed2^0xdddf97f5));
    seed2 = 0xABCB96F7 + (seed2>>1);

    seed1 = 0x4beb5d59*seed1 + 0x2600e1f7;                                // LCG
    seed2 = seed2+0x8009d14b + ((((signed int)seed2)>>31)&0xda879add);    // OWS

    unsigned int MTseed=(seed1 ^ (seed1>>26))+seed2;
    return (MTseed^(MTseed>>20))*0x6957f5a7;
  }

  int Saru_SEED1(unsigned int *iseed1)
  {
    Saru saru(*iseed1);
    saru2=saru;
    return 0;
  }

  int Saru_SEED2(unsigned int *iseed1, unsigned int *iseed2)
  {
    Saru saru(*iseed1,*iseed2);
    saru2=saru;
    return 0;
  }

  int Saru_SEED3(unsigned int *iseed1, unsigned int *iseed2, unsigned int *iseed3)
  {
    Saru saru(*iseed1,*iseed2,*iseed3);
    saru2=saru;
    return 0;
  }

  // Advance state by 1, and output a 32 bit integer pseudo-random value.
  // with variate in [1, high] for high < 2^32
  unsigned int SaruGetIntegerVariate1(unsigned int *high)
  {
    if (*high <= 1) return 1;

    unsigned int usedhigh = *high-1;

    usedhigh |= usedhigh >> 1;
    usedhigh |= usedhigh >> 2;
    usedhigh |= usedhigh >> 4;
    usedhigh |= usedhigh >> 8;
    usedhigh |= usedhigh >> 16;

    // Draw numbers until one is found in [0,n]
    unsigned int i = saru1.u32() & usedhigh;
    while (i > *high-1) {
      i = saru1.u32() & usedhigh;
    }

    // For FORTRAN Index
    return ++i;
  }

  // Advance state by 1, and output a 32 bit integer pseudo-random value.
  // with variate in [1, n] for n < 2^31-2
  unsigned int SaruGetIntegerVariate2()
  {
    // Draw numbers until one is found in [0,n]
    unsigned int i = saru1.u32() & TWO_31;
    while (i > TWO_31) {
      i = saru1.u32() & TWO_31;
    }

    // For FORTRAN Index
    return ++i;
  }

  // Advance state by 1, and output a single precision [0..1) floating point
  float SaruGetRealVariateS1()
  {
    return saru1.f();
  }

  float SaruGetRealVariateS2(float *low, float *high)
  {
    return saru1.f(*low,*high);
  }

  // Advance state by 1, and output a double precision [0..1) floating point
  double SaruGetRealVariateD1()
  {
    return saru1.d();
  }

  // Advance state by 1, and output a double precision [0..1) floating point
  double SaruGetRealVariateD2(double *low, double *high)
  {
    return saru1.d(*low,*high);
  }

  // This stream will be seeded at each step
  float Saru_PRNG()
  {
    return saru2.f();;
  }

  double Saru_PRNGD()
  {
    return saru2.d();;
  }
}

#endif /* PPM_RC_SARUPRNG_H */