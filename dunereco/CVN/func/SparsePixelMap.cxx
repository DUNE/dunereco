////////////////////////////////////////////////////////////////////////
/// \file    SparsePixelMap.cxx
/// \brief   Sparse pixel map for CVN
/// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "dune/CVN/func/SparsePixelMap.h"
#include "canvas/Utilities/Exception.h"

namespace cvn {

  SparsePixelMap::SparsePixelMap(unsigned int dim, unsigned int views, bool usePixelTruth)
    : fDim(dim), fViews(views), fUsePixelTruth(usePixelTruth)
  {
    fCoordinates.resize(fViews);
    fValues.resize(fViews);
    if (fUsePixelTruth) {
      fPixelPDG.resize(fViews);
      fPixelTrackID.resize(fViews);
    }
  }

  /// Default AddHit implementation, which just adds pixel value and coordinates
  void SparsePixelMap::AddHit(unsigned int view, std::vector<unsigned int> coordinates, float value) {

    if (coordinates.size() != fDim) {
      throw art::Exception(art::errors::LogicError)
        << "Coordinate vector with size " << coordinates.size()
        << " does not match sparse pixel map dimension " << fDim;
    }

    if (fUsePixelTruth) {
      throw art::Exception(art::errors::LogicError)
        << "Pixel truth is enabled for this SparsePixelMap, so you must include pixel PDG and "
        << "track ID when calling AddHit.";
    }

    fCoordinates[view].push_back(coordinates);
    fValues[view].push_back(value);
  }

  /// AddHit function that includes per-pixel truth labelling for segmentation
  void SparsePixelMap::AddHit(unsigned int view, std::vector<unsigned int> coordinates,
    float value, int pdg, int trueID) {

    if (coordinates.size() != fDim) {
      throw art::Exception(art::errors::LogicError)
        << "Coordinate vector with size " << coordinates.size()
        << " does not match sparse pixel map dimension " << fDim;
    }

    if (!fUsePixelTruth) {
      throw art::Exception(art::errors::LogicError)
        << "Pixel truth is disabled for this SparsePixelMap, but AddHit call includes "
        << "pixel PDG and track ID.";
    }

    fCoordinates[view].push_back(coordinates);
    fValues[view].push_back(value);
    fPixelPDG[view].push_back(pdg);
    fPixelTrackID[view].push_back(trueID);

  }

  std::vector<unsigned int> SparsePixelMap::GetNPixels() const {
    
    std::vector<unsigned int> ret(fViews);
    for (size_t it = 0; it < fViews; ++it) {
      ret[it] = fValues[it].size();
    }
    return ret;
  }

  // Return flat coordinates vector
  std::vector<std::vector<unsigned int>> SparsePixelMap::GetCoordinatesFlat() const {
    return FlattenVector<std::vector<unsigned int>>(fCoordinates);
  }

  // Return flat pixel values vector
  std::vector<float> SparsePixelMap::GetValuesFlat() const {
    return FlattenVector<float>(fValues);
  }

  // Return flat true PDG vector
  std::vector<int> SparsePixelMap::GetPixelPDGFlat() const {
    return FlattenVector<int>(fPixelPDG);
  }

  // Return flat true G4 track ID vector
  std::vector<int> SparsePixelMap::GetPixelTrackIDFlat() const {
    return FlattenVector<int>(fPixelTrackID);
  }

} // namespace cvn
