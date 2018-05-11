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

  /// Enumeration to describe the order of the TF network output
  typedef enum TFResult
  {
    kTFNumuQE,           ///< Numu CC QE interaction
    kTFNumuRes,          ///< Numu CC Resonant interaction
    kTFNumuDIS,          ///< Numu CC DIS interaction
    kTFNumuOther,        ///< Numu CC, other than above
    kTFNueQE,            ///< Nue CC QE interaction
    kTFNueRes,           ///< Nue CC Resonant interaction
    kTFNueDIS,           ///< Nue CC DIS interaction
    kTFNueOther,         ///< Nue CC, other than above
    kTFNutauQE,          ///< Nutau CC QE interaction
    kTFNutauRes,         ///< Nutau CC Resonant interaction
    kTFNutauDIS,         ///< Nutau CC DIS interaction
    kTFNutauOther,       ///< Nutau CC, other than above
    kTFNC                ///< NC interaction
  } TFResultType;

  // Enumeration to describe the flavour-reduced TF network output
  typedef enum TFFlavour
  {
    kFlavNumuCC,
    kFlavNueCC,
    kFlavNutauCC,
    kFlavNC
  } TFFlavourType; 

}

#endif // CVN_INTERACTIONTYPE_H
