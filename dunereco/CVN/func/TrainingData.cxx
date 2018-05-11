////////////////////////////////////////////////////////////////////////
// \file    TrainingData.h
/// \brief   The TrainingData objects contains a PixelMap and the
///          output class type, and any other bit that goes into the ANN
// \author   radovic -- a.radovic@gmail.com

//#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dune/CVN/func/TrainingData.h"
//#include "MCCheater/BackTracker.h"

namespace cvn
{

  TrainingData::TrainingData(const InteractionType& interaction,
                             float nuEnergy, float lepEnergy,
                             float nueEnergy, float numuEnergy,
                             float weight,
                             const PixelMap& pMap):
  fInt(interaction),
  fNuEnergy(nuEnergy),
  fLepEnergy(lepEnergy),
  fRecoNueEnergy(nueEnergy),
  fRecoNumuEnergy(numuEnergy),
  fEventWeight(weight),
  fPMap(pMap)
  {  }


  void TrainingData::FillOutputVector(float* output) const
  {
    for(unsigned int i = 0; i < kNIntType; ++i)
      output[i] = 0;

    output[fInt] = 1;



  }

} // end namespace cvn
////////////////////////////////////////////////////////////////////////
