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
    //!  Author           - y.afshar           March 2016
    //!-------------------------------------------------------------------------

#include <stdint.h>

static const uint32_t c1 = 0xcc9e2d51;
static const uint32_t c2 = 0x1b873593;

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here
inline __attribute__((always_inline)) uint32_t getblock32 (const uint32_t * p, int i)
{
  return p[i];
}

//-----------------------------------------------------------------------------
extern "C" {

  uint32_t HashFunc (const void *key, uint32_t *seed, uint32_t *tablenrow)
  {
    const uint8_t * data = (const uint8_t*)key;

    uint32_t h1 = *seed;

    //----------
    // body
    const uint32_t *blocks = (const uint32_t *)(data + 8);

    // 64 bits or 2 blocks of 32
    for (int i = -2; i; i++) {

      uint32_t k1 = getblock32(blocks,i);

      k1 *= c1;
      k1  = (k1 << 15) | (k1 >> 17);
      k1 *= c2;

      h1 ^= k1;
      h1  = (h1 << 13) | (h1 >> 19);
      h1  = h1*5+0xe6546b64;
    }

    //----------
    // finalization
    h1 ^= 8;
    // Finalization mix - force all bits of a hash block to avalanche
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    h1 &= *tablenrow-1;

    return h1;
  }

  uint32_t HashFunc_XY (uint32_t *X, uint32_t *Y, uint32_t *seed, uint32_t *tablenrow, int64_t *XY)
  {
    //this is fine as *X and *Y are positive numbers and less than 2^16 (<=65535)
    *XY = *Y<<16 | *X;

    const uint8_t * data = (const uint8_t*)XY;

    uint32_t h1 = *seed;

    //----------
    // body
    const uint32_t *blocks = (const uint32_t *)(data + 8);

    // 64 bits or 2 blocks of 32
    for (int i = -2; i; i++) {

      uint32_t k1 = getblock32(blocks,i);

      k1 *= c1;
      k1  = (k1 << 15) | (k1 >> 17);
      k1 *= c2;

      h1 ^= k1;
      h1  = (h1 << 13) | (h1 >> 19);
      h1  = h1*5+0xe6546b64;
    }

    //----------
    // finalization
    h1 ^= 8;
    // Finalization mix - force all bits of a hash block to avalanche
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    h1 &= *tablenrow-1;

    return h1;
  }

  uint32_t HashFunc_XYLabel (uint32_t *X, uint32_t *Y, uint32_t *Label, uint32_t *seed, uint32_t *tablenrow, int64_t *XYLabel)
  {
    // Label is less than 2^31 and this shifting is safe here
    *XYLabel  = 0;
    *XYLabel |= *Label;
    //this is fine as *X and *Y are positive numbers and less than 2^16 (<=65535)
    *XYLabel  = (*XYLabel << 32) | *Y << 16 | *X;

    const uint8_t * data = (const uint8_t*)XYLabel;

    uint32_t h1 = *seed;

    //----------
    // body
    const uint32_t *blocks = (const uint32_t *)(data + 8);

    // (8 Bytes) 64 bits or 2 blocks of (4 Bytes) 32 bits
    for (int i = -2; i; i++) {
      uint32_t k1 = getblock32(blocks,i);

      k1 *= c1;
      k1  = (k1 << 15) | (k1 >> 17);
      k1 *= c2;

      h1 ^= k1;
      h1  = (h1 << 13) | (h1 >> 19);
      h1  = h1*5+0xe6546b64;
    }

    //----------
    // finalization
    h1 ^= 8;
    // Finalization mix - force all bits of a hash block to avalanche
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    h1 &= *tablenrow-1;

    return h1;
  }

  uint32_t HashFunc_XYZ (uint32_t *X, uint32_t *Y, uint32_t *Z, uint32_t *seed, uint32_t *tablenrow, int64_t *XYZ)
  {
    *XYZ  = 0;
    *XYZ |= *Z;
    //this is fine as *X and *Y are positive numbers and less than 2^11 (<=2047)
    *XYZ  = *XYZ << 22 | *Y << 11 | *X;

    const uint8_t * data = (const uint8_t*)XYZ;

    uint32_t h1 = *seed;

    //----------
    // body
    const uint32_t *blocks = (const uint32_t *)(data + 8);

    // (8 Bytes) 64 bits or 2 blocks of (4 Bytes) 32 bits
    for (int i = -2; i; i++) {

      uint32_t k1 = getblock32(blocks,i);

      k1 *= c1;
      k1  = (k1 << 15) | (k1 >> 17);
      k1 *= c2;

      h1 ^= k1;
      h1  = (h1 << 13) | (h1 >> 19);
      h1  = h1*5+0xe6546b64;
    }

    //----------
    // finalization
    h1 ^= 8;
    // Finalization mix - force all bits of a hash block to avalanche
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    h1 &= *tablenrow-1;

    return h1;
  }

  uint32_t HashFunc_XYZLabel (uint32_t *X, uint32_t *Y, uint32_t *Z, uint32_t *Label, uint32_t *seed, uint32_t *tablenrow, int64_t *XYZLabel)
  {
    // Label is less than 2^31 and this shifting is safe here
    *XYZLabel  =0;
    *XYZLabel |= *Label;
    *XYZLabel  = *XYZLabel << 11 | *Z;
    //this is fine as *X and *Y are positive numbers and less than 2^11 (<=2047)
    *XYZLabel  = *XYZLabel << 22 | *Y << 11 | *X;

    const uint8_t * data = (const uint8_t*)XYZLabel;

    uint32_t h1 = *seed;

    //----------
    // body
    const uint32_t *blocks = (const uint32_t *)(data + 8);

    //(8 Bytes) 64 bits or 2 blocks of (4 Bytes) 32 bits
    for (int i = -2; i; i++) {

      uint32_t k1 = getblock32(blocks,i);

      k1 *= c1;
      k1  = (k1 << 15) | (k1 >> 17);
      k1 *= c2;

      h1 ^= k1;
      h1  = (h1 << 13) | (h1 >> 19);
      h1  = h1*5+0xe6546b64;
    }

    //----------
    // finalization
    h1 ^= 8;
    // Finalization mix - force all bits of a hash block to avalanche
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    h1 &= *tablenrow-1;

    return h1;
  }

}