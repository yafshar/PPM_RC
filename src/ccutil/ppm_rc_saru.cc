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

#include "saruprng.h"

// 2^31-1
static const unsigned int TWO_31=2147483646;

Saru s;
Saru n;

extern "C" {
  int SaruInitialize1(unsigned int *iseed)
  {
    s.advance(*iseed);
    return 0;
  }

  int SaruInitialize2()
  {
    n.fork<123>();
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
    unsigned int i = s.u32() & usedhigh;
    while (i > *high-1) {
      i = s.u32() & usedhigh;
    }
    return ++i;
  }

  // Advance state by 1, and output a 32 bit integer pseudo-random value.
  // with variate in [1, n] for n < 2^31-2
  unsigned int SaruGetIntegerVariate2()
  {
    // Draw numbers until one is found in [0,n]
    unsigned int i = s.u32() & TWO_31;
    while (i > TWO_31) {
      i = s.u32() & TWO_31;
    }
    return ++i;
  }

  // Advance state by 1, and output a single precision [0..1) floating point
  float SaruGetRealVariateS1()
  {
    return n.f();
  }

  float SaruGetRealVariateS2(float *low, float *high)
  {
    return n.f(*low,*high);
  }

  // Advance state by 1, and output a double precision [0..1) floating point
  double SaruGetRealVariateD1()
  {
    return n.d();
  }

  // Advance state by 1, and output a double precision [0..1) floating point
  double SaruGetRealVariateD2(double *low, double *high)
  {
    return n.d(*low,*high);
  }
}



