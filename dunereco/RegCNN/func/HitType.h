////////////////////////////////////////////////////////////////////////
// \file    HitType.h
///\brief   Defines an enumeration for cellhit classification
///
// \author psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef REGCNN_HITTYPE_H
#define REGCNN_HITTYPE_H


namespace cnn
{

  typedef enum HType
  {
    kElectronHit,
    kMuonHit,
    kProtonHit,
    kNeutronHit,
    kPionHit,
    kPiZeroHit,
    kGammaHit,
    kOtherPDGhit,
    kUnknownHit,
    kEmptyHit
  } HitType;

}

#endif // REGCNN_HITTYPE_H
