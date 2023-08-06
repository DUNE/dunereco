////////////////////////////////////////////////////////////////////////
// \file    TrainingData.h
/// \brief   The TrainingData objects contains a PixelMap and the
///          output class type, and any other bit that goes into the ANN
// \author   radovic -- a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef CVN_TRAININGDATA_H
#define CVN_TRAININGDATA_H

#include "dunereco/CVN/func/InteractionType.h"
#include "dunereco/CVN/func/AssignLabels.h"

namespace cvn
{


  /// \brief   The TrainingData objects contains a PixelMap and the
  ///          output class type, and any other bit that goes into the ANN

  class TrainingData
  {

  public:
    TrainingData(){};
    TrainingData(const InteractionType& interaction,
                 float nuEnergy, float lepEnergy, float lepAngle,
                 float nueEnergy, float numuEnergy,
                 float nutauEnergy, float weight,
                 const PixelMap& pMap);

    unsigned int NOutput() const {return (unsigned int)kNIntType;};

    void FillOutputVector(float* output) const;

    // Set topology information separately to save having a large number of 
    // arguments in the constructor.
    void SetTopologyInformation(int pdg, int nproton, int npion,
                                int npizero, int nneutron, int toptype,
                                int toptypealt);

    InteractionType  fInt;     ///< Class of the event
    float    fNuEnergy;        ///< True energy of neutrino event
    float    fLepEnergy;       ///< True energy of outgoing lepton
    float    fLepAngle;        ///< True lepton angle wrt neutrino direction 
    float    fRecoNueEnergy;   ///< Reconstructed energy under nue hypothesis
    float    fRecoNumuEnergy;  ///< Reconstructed energy under numu hypothesis
    float    fRecoNutauEnergy; ///< Reconstructed energy under nutau hypothesis
    float    fEventWeight;     ///< The event weight (norm * oscProb)
 
    // If we are using topology information, store it here
    bool fUseTopology;
    int  fNuPDG;
    int  fNProton;
    int  fNPion;
    int  fNPizero;
    int  fNNeutron;
    int  fTopologyType;
    int  fTopologyTypeAlt;

    PixelMap fPMap;           ///< PixelMap for the event
  };

} // end namespace

#endif // CVN_TRAININGDATA_H
//////////////////////////////////////////////////////////////////////////////
