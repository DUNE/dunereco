////////////////////////////////////////////////////////////////////////
/// \file    CTPResult.h
/// \brief   Class storing the result from the convolutional track PID
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CTPRESULT_H
#define CTPRESULT_H

#include <vector>

namespace ctp
{

  /// Class containing some utility functions for all things CVN
  class CTPResult
  {
  public:
    CTPResult();
    CTPResult(const std::vector<float> &vals);
    ~CTPResult();

    // Function to calculate the PID for a given track
    std::vector<float> GetResults() const;

    bool IsValid() const;
    
    // Individual scores
    float GetMuonScore() const {return fMuonScore;};
    float GetPionScore() const {return fPionScore;};
    float GetProtonScore() const {return fProtonScore;};

  private:

    float fMuonScore;
    float fPionScore;
    float fProtonScore;

  };

}

#endif  // CTPRESULT_H
