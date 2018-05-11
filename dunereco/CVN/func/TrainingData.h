////////////////////////////////////////////////////////////////////////
// \file    TrainingData.h
/// \brief   The TrainingData objects contains a PixelMap and the
///          output class type, and any other bit that goes into the ANN
// \author   radovic -- a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef CVN_TRAININGDATA_H
#define CVN_TRAININGDATA_H

#include "dune/CVN/func/InteractionType.h"

namespace cvn
{


  /// \brief   The TrainingData objects contains a PixelMap and the
  ///          output class type, and any other bit that goes into the ANN

  class TrainingData
  {

  public:
    TrainingData(){};
    TrainingData(const InteractionType& interaction,
                 float nuEnergy, float lepEnergy,
                 float nueEnergy, float numuEnergy,
                 float weight,
                 const PixelMap& pMap);

    unsigned int NOutput() const {return (unsigned int)kNIntType;};

    void FillOutputVector(float* output) const;

    InteractionType  fInt;    ///< Class of the event
    float    fNuEnergy;       ///< True energy of neutrino event
    float    fLepEnergy;      ///< True energy of outgoing lepton
    float    fRecoNueEnergy;  ///< Reconstructed energy under nue hypothesis
    float    fRecoNumuEnergy; ///< Reconstructed energy under nue hypothesis
    float    fEventWeight;    ///< The event weight (norm * oscProb)
    PixelMap fPMap;           ///< PixelMap for the event
  };

} // end namespace

#endif // CVN_TRAININGDATA_H
//////////////////////////////////////////////////////////////////////////////
