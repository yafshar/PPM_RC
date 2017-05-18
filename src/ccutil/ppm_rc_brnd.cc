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
    //!  Author           - y.afshar           Dec   2015
    //!-------------------------------------------------------------------------
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

boost::mt19937 NumberGeneratorBoost;

typedef boost::random::discrete_distribution<> DiscrDistrType;

DiscrDistrType *ImageDiscrDistr=NULL;
DiscrDistrType *PartDiscrDistr=NULL;

template <typename T>
class ArrayWrapper
{
  class iterator : public std::iterator<std::input_iterator_tag, T>
  {
  public:
    iterator(const T *aPointer) : iPosition(aPointer) {}

    bool operator==(const iterator& rhs) {return iPosition == rhs.iPosition;}
    bool operator!=(const iterator& rhs) {return iPosition != rhs.iPosition;}

    void operator++() {++iPosition;}

    T operator*() {
      return *iPosition;
    }

    ~iterator() {
      // nothing to do
    }

  private:
    const T *iPosition;
  };

public:
  ArrayWrapper(const T *aInputArray, long aNumOfElements) : iArray(aInputArray), iNumOfElements(aNumOfElements) {}

  iterator begin() {
    return iterator(iArray);
  }
  iterator end() {
    return iterator(iArray + iNumOfElements);
  }

private:
  const T *iArray;
  long iNumOfElements;
};

extern "C" {
  int MTInitialize(unsigned int *iseed)
  {
    //seed used on construction or seeding
    NumberGeneratorBoost.seed(*iseed);
    return 0;
  }

  int GenerateFloatImageDiscrDistr2D(float *image, unsigned int *width, unsigned int *length)
  {
    long int l=(long int) *width * (long int) *length;
    ArrayWrapper<float> aw(image, l);

    // Sets the parameters of the distribution.
    // param_type --- Constructs a discrete_distribution from an iterator range (vVec.begin(), vVec.end())
    // the values of the range represent weights for the possible values of the distribution
    // ImageDiscrDistr.param(DiscrDistrType::param_type(aw.begin(), aw.end()));
    ImageDiscrDistr= new DiscrDistrType(aw.begin(), aw.end());

    return 0;
  }

  int GenerateDoubleImageDiscrDistr2D(double *image, unsigned int *width, unsigned int *length)
  {
    long int l=(long int) *width * (long int) *length;
    ArrayWrapper<double> aw(image, l);

    // Sets the parameters of the distribution.
    // param_type --- Constructs a discrete_distribution from an iterator range (vVec.begin(), vVec.end())
    // the values of the range represent weights for the possible values of the distribution
    // ImageDiscrDistr.param(DiscrDistrType::param_type(aw.begin(), aw.end()));
    ImageDiscrDistr= new DiscrDistrType(aw.begin(), aw.end());

    return 0;
  }

  int GenerateFloatImageDiscrDistr3D(float *image, unsigned int *width, unsigned int *length, unsigned int *depth)
  {
    long int l=(long int) *width * (long int) *length * (long int) *depth;
    ArrayWrapper<float> aw(image, l);

    // Sets the parameters of the distribution.
    // param_type --- Constructs a discrete_distribution from an iterator range (vVec.begin(), vVec.end())
    // the values of the range represent weights for the possible values of the distribution
    // ImageDiscrDistr.param(DiscrDistrType::param_type(aw.begin(), aw.end()));
    ImageDiscrDistr= new DiscrDistrType(aw.begin(), aw.end());

    return 0;
  }

  int GenerateDoubleImageDiscrDistr3D(double *image, unsigned int *width, unsigned int *length, unsigned int *depth)
  {
    long int l=(long int) *width * (long int) *length * (long int) *depth;
    ArrayWrapper<double> aw(image, l);

    // Sets the parameters of the distribution.
    // param_type --- Constructs a discrete_distribution from an iterator range (vVec.begin(), vVec.end())
    // the values of the range represent weights for the possible values of the distribution
    // ImageDiscrDistr.param(DiscrDistrType::param_type(aw.begin(), aw.end()));
    ImageDiscrDistr= new DiscrDistrType(aw.begin(), aw.end());

    return 0;
  }

  long int GetImageDistrIndex()
  {
    return (*ImageDiscrDistr)(NumberGeneratorBoost)+1;
  }

  int DestroyImageDiscrDistr()
  {
    if (ImageDiscrDistr != NULL) {
      delete ImageDiscrDistr;
      ImageDiscrDistr=NULL;
    }
    return 0;
  }

  int GenerateFloatParticlesFwdProposalsDiscrDistr(float *part, unsigned int *size)
  {
    if (PartDiscrDistr != NULL) {
      delete PartDiscrDistr;
      PartDiscrDistr=NULL;
    }
    ArrayWrapper<float> aw(part, *size);
    PartDiscrDistr = new DiscrDistrType(aw.begin(),aw.end());
    return 0;
  }

  int GenerateDoubleParticlesFwdProposalsDiscrDistr(double *part, unsigned int *size)
  {
    if (PartDiscrDistr != NULL) {
      delete PartDiscrDistr;
      PartDiscrDistr=NULL;
    }
    ArrayWrapper<double> aw(part, *size);
    PartDiscrDistr = new DiscrDistrType(aw.begin(),aw.end());
    return 0;
  }

  unsigned int GetPartDistrIndex()
  {
    return (*PartDiscrDistr)(NumberGeneratorBoost)+1;
  }

  int DestroyParticlesDiscrDistr()
  {
    if (PartDiscrDistr != NULL) {
      delete PartDiscrDistr;
      PartDiscrDistr=NULL;
    }
    return 0;
  }
}

