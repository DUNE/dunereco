////////////////////////////////////////////////////////////////////////
/// \file    SparsePixelMap.h
/// \brief   Sparse pixel map for CVN
/// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>

namespace cvn
{

  class SparsePixelMap
  {
  public:
    SparsePixelMap(unsigned int dim, unsigned int views, bool usePixelTruth=false);
    SparsePixelMap() {};
    ~SparsePixelMap() {};

    void AddHit(unsigned int view, std::vector<float> coordinates, std::vector<float> features);
    void AddHit(unsigned int view, std::vector<float> coordinates, std::vector<float> features, std::vector<int> pdgs, 
      std::vector<int> tracks, std::vector<float> energies, std::vector<std::string> processes);
    unsigned int GetDim() const { return fDim; };
    unsigned int GetViews() const { return fViews; }
    std::vector<unsigned int> GetNPixels() const;
    unsigned int GetNPixels(size_t view) const { return fFeatures[view].size(); };

    std::vector<std::vector<std::vector<float>>> GetCoordinates() const { return fCoordinates; };
    std::vector<std::vector<float>> GetCoordinates(size_t view) const { return fCoordinates[view]; };

    std::vector<std::vector<std::vector<float>>> GetFeatures() const { return fFeatures; };
    std::vector<std::vector<float>> GetFeatures(size_t view) const { return fFeatures[view]; };

    std::vector<std::vector<std::vector<int>>> GetPixelPDGs() const { return fPixelPDGs; };
    std::vector<std::vector<int>> GetPixelPDGs(size_t view) const { return fPixelPDGs[view]; };

    std::vector<std::vector<std::vector<int>>> GetPixelTrackIDs() const {return fPixelTrackIDs; };
    std::vector<std::vector<int>> GetPixelTrackIDs(size_t view) const {return fPixelTrackIDs[view]; }; 

    std::vector<std::vector<std::vector<float>>> GetPixelEnergies() const {return fPixelEnergies; };
    std::vector<std::vector<float>> GetPixelEnergies(size_t view) const {return fPixelEnergies[view]; }; 

    std::vector<std::vector<std::vector<std::string>>> GetProcesses() const { return fProcesses; };
    std::vector<std::vector<std::string>> GetProcesses(size_t view) const { return fProcesses[view]; };



  private:

    template <class T>
    std::vector<T> FlattenVector(std::vector<std::vector<T>> vec) const {
      std::vector<T> ret;
      for (auto it : vec) {
        ret.insert(ret.end(), it.begin(), it.end());
      }
      return ret;
    }

    unsigned int fDim; ///< Dimensionality of each pixel map
    unsigned int fViews; ///< Number of views
    bool fUsePixelTruth; ///< Whether to use a per-pixel ground truth for pixel segmentation
    std::vector<std::vector<std::vector<float>>> fCoordinates; ///< Coordinates of non-zero pixels
    std::vector<std::vector<std::vector<float>>> fFeatures; ///< Features of non-zero pixels
    std::vector<std::vector<std::vector<int>>> fPixelPDGs; ///< True particle PDG responsible for pixel
    std::vector<std::vector<std::vector<int>>> fPixelTrackIDs; ///< G4 track IDs responsible for pixelel
    std::vector<std::vector<std::vector<float>>> fPixelEnergies;
    std::vector<std::vector<std::vector<std::string>>> fProcesses; // physics process that created the particle

  }; // class SparsePixelMap
} // namespace cvn
