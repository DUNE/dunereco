////////////////////////////////////////////////////////////////////////
/// \file    SparsePixelMap.cxx
/// \brief   Sparse pixel map for CVN
/// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "dune/CVN/func/SparsePixelMap.h"
#include "canvas/Utilities/Exception.h"

namespace cvn {

  SparsePixelMap::SparsePixelMap(unsigned int dim, unsigned int views)
    : fDim(dim), fViews(views)
  {
    fCoordinates.resize(fViews);
    fValues.resize(fViews);
  }

  void SparsePixelMap::AddHit(unsigned int view, std::vector<unsigned int> coordinates, float value) {

    if (coordinates.size() != fDim)
      throw art::Exception(art::errors::LogicError)
        << "Coordinate vector with size " << coordinates.size()
        << " does not match sparse pixel map dimension " << fDim;

    fCoordinates[view].push_back(coordinates);
    fValues[view].push_back(value);
  }

  std::vector<unsigned int> SparsePixelMap::GetNPixels() const {
    
    std::vector<unsigned int> ret(fViews);
    for (size_t it = 0; it < fViews; ++it) {
      ret[it] = fValues[it].size();
    }
    return ret;
  }

  std::vector<std::vector<unsigned int>> SparsePixelMap::GetCoordinatesFlat() const {

    std::vector<std::vector<unsigned int>> ret;
    for (auto it : fCoordinates) {
      ret.insert(ret.end(), it.begin(), it.end());
    }
    return ret;
  }

  std::vector<float> SparsePixelMap::GetValuesFlat() const {

    std::vector<float> ret;
    for (auto it : fValues) {
      ret.insert(ret.end(), it.begin(), it.end());
    }
    return ret;
  }

} // namespace cvn
