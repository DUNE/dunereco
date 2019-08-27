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
    void AddHit(unsigned int view, std::vector<unsigned int> coordinates, float value, int pdg, int trueID);
    unsigned int GetDim() const { return fDim; };
    unsigned int GetViews() const { return fViews; }
    std::vector<unsigned int> GetNPixels() const;
    unsigned int GetNPixels(size_t view) const { return fValues[view].size(); };

    std::vector<std::vector<std::vector<unsigned int>>> GetCoordinates() const { return fCoordinates; };
    std::vector<std::vector<unsigned int>> GetCoordinates(size_t view) const { return fCoordinates[view]; };
    std::vector<std::vector<unsigned int>> GetCoordinatesFlat() const;

    std::vector<std::vector<float>> GetValues() const { return fValues; };
    std::vector<float> GetValues(size_t view) const { return fValues[view]; };
    std::vector<float> GetValuesFlat() const;

    std::vector<std::vector<int>> GetPixelPDG() const { return fPixelPDG; };
    std::vector<int> GetPixelPDG(size_t view) const { return fPixelPDG[view]; };
    std::vector<int> GetPixelPDGFlat() const;

    std::vector<std::vector<int>> GetPixelTrackID() const { return fPixelTrackID; };
    std::vector<int> GetPixelTrackID(size_t view) const { return fPixelTrackID[view]; };
    std::vector<int> GetPixelTrackIDFlat() const;

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
    std::vector<std::vector<int>> fPixelPDG; ///< True particle PDG responsible for pixel
    std::vector<std::vector<int>> fPixelTrackID; ///< G4 track ID responsible for pixel

  }; // class SparsePixelMap
} // namespace cvn
