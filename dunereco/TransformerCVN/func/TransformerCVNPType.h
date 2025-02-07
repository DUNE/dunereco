////////////////////////////////////////////////////////////////////////
// \file    ProngType.h
///\brief   Defines an enumeration for prong classification
////////////////////////////////////////////////////////////////////////
#ifndef TRANSFORMERCVN_PRONGTYPE_H
#define TRANSFORMERCVN_PRONGTYPE_H

namespace cnn
{

  typedef enum PType
  {
    kElectron,
    kMuon,
    kProton,
    kNeutron,
    kPion,
    kPiZero,
    kGamma,
    kOtherPDG,
    kUnknown,
    kEmpty,
    kEM,
    kHadron
  } ProngType;

  int GetPDGByPType(PType ptype);

}

#endif // TRANSFORMERCVN_PRONGTYPE_H
