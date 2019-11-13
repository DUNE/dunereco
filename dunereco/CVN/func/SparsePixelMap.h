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

    void AddHit(unsigned int view, std::vector<unsigned int> coordinates, float value);
    void AddHit(unsigned int view, std::vector<unsigned int> coordinates, float value, std::vector<int> pdgs, 
      std::vector<int> TrackIDs, std::vector<float> Energy, std::vector<int> pdgParent, std::vector<std::string> process);
    unsigned int GetDim() const { return fDim; };
    unsigned int GetViews() const { return fViews; }
    std::vector<unsigned int> GetNPixels() const;
    unsigned int GetNPixels(size_t view) const { return fValues[view].size(); };

    std::vector<std::vector<std::vector<unsigned int>>> GetCoordinates() const { return fCoordinates; };
    std::vector<std::vector<unsigned int>> GetCoordinates(size_t view) const { return fCoordinates[view]; };


    std::vector<std::vector<float>> GetValues() const { return fValues; };
    std::vector<float> GetValues(size_t view) const { return fValues[view]; };


    std::vector<std::vector<std::vector<int>>> GetPixelPDGs() const { return fPixelPDGs; };
    std::vector<std::vector<int>> GetPixelPDGs(size_t view) const { return fPixelPDGs[view]; };


    std::vector<std::vector<std::vector<int>>> GetPixelTrackIDs() const {return fPixelTrackIDs; };
    std::vector<std::vector<int>> GetPixelTrackIDs(size_t view) const {return fPixelTrackIDs[view]; }; 


    std::vector<std::vector<std::vector<float>>> GetPixelEnergy() const {return fPixelEnergy; };
    std::vector<std::vector<float>> GetPixelEnergy(size_t view) const {return fPixelEnergy[view]; }; 


    // Particle Parents Information  
    std::vector<std::vector<std::vector<int>>> GetParentPDGs() const { return fParentPDG; };
    std::vector<std::vector<int>> GetParentPDGs(size_t view) const { return fParentPDG[view]; };

    std::vector<std::vector<std::vector<std::string>>> GetProcess() const { return fProcess; };
    std::vector<std::vector<std::string>> GetProcess(size_t view) const { return fProcess[view]; };



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
    std::vector<std::vector<std::vector<unsigned int>>> fCoordinates; ///< Coordinates of non-zero pixels
    std::vector<std::vector<float>> fValues; ///< Values of non-zero pixels
    std::vector<std::vector<std::vector<int>>> fPixelPDGs; ///< True particle PDG responsible for pixel
    std::vector<std::vector<std::vector<int>>> fPixelTrackIDs; ///< G4 track IDs responsible for pixelel
    std::vector<std::vector<std::vector<float>>> fPixelEnergy;

    std::vector<std::vector<std::vector<int>>> fParentPDG; // particle's parent
    std::vector<std::vector<std::vector<std::string>>> fProcess; // physics process that created the particle
 


  }; // class SparsePixelMap
} // namespace cvn
