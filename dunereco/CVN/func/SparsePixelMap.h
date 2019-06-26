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
    SparsePixelMap(unsigned int dim, unsigned int views);
    SparsePixelMap() {};
    ~SparsePixelMap() {};

    void AddHit(unsigned int view, std::vector<unsigned int> coordinates, float value);
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

  private:
    unsigned int fDim; ///< Dimensionality of each pixel map
    unsigned int fViews; ///< Number of views
    std::vector<std::vector<std::vector<unsigned int>>> fCoordinates; ///< Coordinates of non-zero pixels
    std::vector<std::vector<float>> fValues; ///< Values of non-zero pixels

  }; // class SparsePixelMap
} // namespace cvn
