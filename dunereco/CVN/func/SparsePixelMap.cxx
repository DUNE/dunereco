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
    fValues.resize(fViews);
    if (fUsePixelTruth) {
      fPixelPDGs.resize(fViews);
      fPixelTrackIDs.resize(fViews);
      fPixelEnergy.resize(fViews);
// *******************************************
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
    float value, std::vector<int> pdgs, std::vector<int> Tracks, std::vector<float> Energy) {

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
    fValues[view].push_back(value);
    fPixelPDGs[view].push_back(pdgs);
    fPixelTrackIDs[view].push_back(Tracks);
    fPixelEnergy[view].push_back(Energy);
    //std::cout<< "*Carlos*  view " <<  view << "trackIDs size " << TrackIDs.size() << std::endl;
    //for (auto k : TrackIDs){
      // fPixelTrackIDs[view].push_back(k);

    //}
    
    //std::cout<< "Tracks IDs = " << TrackIDs.at(TrackIDs.size()-1) << std::endl;
    //std::cout<< "TrueID = " << trueID<< std::endl;

     
  }

  std::vector<unsigned int> SparsePixelMap::GetNPixels() const {
    
    std::vector<unsigned int> ret(fViews);
    for (size_t it = 0; it < fViews; ++it) {
      ret[it] = fValues[it].size();
    }
    return ret;
  }

  // Return flat coordinates vector
 // std::vector<std::vector<unsigned int>> SparsePixelMap::GetCoordinatesFlat() const {
 //   return FlattenVector<std::vector<unsigned int>>(fCoordinates);
 // }

  // Return flat pixel values vector
  //std::vector<float> SparsePixelMap::GetValuesFlat() const {
  //  return FlattenVector<float>(fValues);
 // }

  // Return flat particle PDG vector
  //std::vector<std::vector<int>> SparsePixelMap::GetPixelPDGsFlat() const {
  //  return FlattenVector<std::vector<int>>(fPixelPDGs);
  //}

  // Return flat true G4 track ID vector
  //std::vector<int> SparsePixelMap::GetPixelTrackIDsFlat() const {
   // return FlattenVector<int>(fPixelTrackID);
  //}


} // namespace cvn
