////////////////////////////////////////////////////////////////////////
/// \file    SparsePixelMap.cxx
/// \brief   Sparse pixel map for CVN
/// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include "dune/CVN/func/SparsePixelMap.h"
#include "canvas/Utilities/Exception.h"


namespace cvn {

  SparsePixelMap::SparsePixelMap(unsigned int dim, unsigned int views, bool usePixelTruth)
    : fDim(dim), fViews(views), fUsePixelTruth(usePixelTruth)
  {
    fCoordinates.resize(fViews);
    fFeatures.resize(fViews);
    if (fUsePixelTruth) {
      fPixelPDGs.resize(fViews);
      fPixelTrackIDs.resize(fViews);
      fPixelEnergies.resize(fViews);
      fProcesses.resize(fViews);
// *******************************************
    }
  }

  /// Default AddHit implementation, which just adds pixel value and coordinates
  void SparsePixelMap::AddHit(unsigned int view, std::vector<float> coordinates,
    std::vector<float> features) {

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
    fFeatures[view].push_back(features);
  }

  /// AddHit function that includes per-pixel truth labelling for segmentation
  void SparsePixelMap::AddHit(unsigned int view, std::vector<float> coordinates,
    std::vector<float> features, std::vector<int> pdgs, std::vector<int> tracks,
    std::vector<float> energies, std::vector<std::string> processes) {

    if (coordinates.size() != fDim) {
      throw art::Exception(art::errors::LogicError)
        << "Coordinate vector with size " << coordinates.size()
        << " does not match sparse pixel map dimension " << fDim;
    }

    if (!fUsePixelTruth) {
      throw art::Exception(art::errors::LogicError)
        << "Pixel truth is disabled for this SparsePixelMap, but AddHit call includes "
        << "pixel PDG, track ID and Energy";
    }

    fCoordinates[view].push_back(coordinates);
    fFeatures[view].push_back(features);
    fPixelPDGs[view].push_back(pdgs);
    fPixelTrackIDs[view].push_back(tracks);
    fPixelEnergies[view].push_back(energies);
    fProcesses[view].push_back(processes);

  }

  std::vector<unsigned int> SparsePixelMap::GetNPixels() const {
    
    std::vector<unsigned int> ret(fViews);
    for (size_t it = 0; it < fViews; ++it) {
      ret[it] = fFeatures[it].size();
    }
    return ret;
  }

} // namespace cvn
