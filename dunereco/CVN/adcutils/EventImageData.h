//////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       EventImageData
//
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan (Dorota.Stefan@cern.ch), May 2016
//              Stripped from SPMultiTpcDump by Leigh Whitehead (leigh.howard.whitehead@cern.ch)
//
// Produces "global image" of events in multi-tpc detector, could be FD workspace or ProtoDUNE (SP).
// Projections of ADC are saved to 1-byte / pixel (TH1C) format, int-size (TH2I) projection of PDG
// truth and float-size (TH2F) projection of true charge deposits can be saved as well. We use this
// to dump deconvoluted ADC images for preparation of various classifiers and CNN based pattern reco.
//
// Moved from larreco/RecoAlg/ImagePatternAlgs since it became specific to DUNE geometries.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENTIMAGEDATA_HH
#define EVENTIMAGEDATA_HH

#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// C++ Includes
#include <vector>

namespace nnet
{

  class EventImageData // full image for one plane
  {
  public:
    EventImageData(size_t w, size_t d, bool saveDep);

    void addTpc(const TrainingDataAlg & dataAlg, size_t gw, bool flipw, size_t gd, bool flipd);
    bool findCrop(size_t max_area_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const;

    const std::vector< std::vector<float> > & adcData(void) const { return fAdc; }
    const std::vector<float> & wireAdc(size_t widx) const { return fAdc[widx]; }

    const std::vector< std::vector<float> > & depData(void) const { return fDeposit; }
    const std::vector<float> & wireDep(size_t widx) const { return fDeposit[widx]; }

    const std::vector< std::vector<int> > & pdgData(void) const { return fPdg; }
    const std::vector<int> & wirePdg(size_t widx) const { return fPdg[widx]; }

    void setProjXY(const TrainingDataAlg & dataAlg, float x, float y, size_t gw, bool flipw, size_t gd, bool flipd);
    float getProjX(void) const { return fProjX; }
    float getProjY(void) const { return fProjY; }

    int getVtxX(void) const { return fVtxX; }
    int getVtxY(void) const { return fVtxY; }

  private:
    std::vector< std::vector<float> > fAdc, fDeposit;
    std::vector< std::vector<int> > fPdg;
    int fVtxX, fVtxY;
    float fProjX, fProjY;
    bool fSaveDep;
  };

}

#endif

