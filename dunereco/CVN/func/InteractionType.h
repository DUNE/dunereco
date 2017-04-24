////////////////////////////////////////////////////////////////////////
// \file    Interaction.h
///\brief   Defines an enumeration for interaction type
///
// \author   radovic -- a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef CVN_INTERACTION_H
#define CVN_INTERACTION_H

#include "dune/CVN/func/PixelMap.h"

namespace cvn
{

  typedef enum Interaction
  {
    kNumuQE,           ///< Numu CC QE interaction
    kNumuRes,          ///< Numu CC Resonant interaction
    kNumuDIS,          ///< Numu CC DIS interaction
    kNumuOther,        ///< Numu CC, other than above
    kNueQE,            ///< Nue CC QE interaction
    kNueRes,           ///< Nue CC Resonant interaction
    kNueDIS,           ///< Nue CC DIS interaction
    kNueOther,         ///< Nue CC, other than above
    kNutauQE,          ///< Nutau CC QE interaction
    kNutauRes,         ///< Nutau CC Resonant interaction
    kNutauDIS,         ///< Nutau CC DIS interaction
    kNutauOther,       ///< Nutau CC, other than above
    kNuElectronElastic,///< NC Nu On E Scattering
    kNC,               ///< NC interaction
    kCosmic,           ///< Cosmic ray background
    kOther,            ///< Something else.  Tau?  Hopefully we don't use this
    kNIntType          ///< Number of interaction types, used like a vector size
  } InteractionType;


}

#endif // CVN_INTERACTIONTYPE_H
