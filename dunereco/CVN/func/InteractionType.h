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

  // It might be a good idea to consider topology information. 
  // We should use base two here so we can & and | results.
  typedef enum Topology
  {
    kTopUnset   = 0x00000,
    // Flavours
    kTopNumu    = 0x00001,
    kTopNue     = 0x00002,
    kTopNutau   = 0x00004,
    kTopNC      = 0x00008,
    // Topologies - protons
    kTop0proton = 0x00010,
    kTop1proton = 0x00020,
    kTop2proton = 0x00040,
    kTopNproton = 0x00080,
    // Topologies - pions
    kTop0pion   = 0x00100,
    kTop1pion   = 0x00200,
    kTop2pion   = 0x00400,
    kTopNpion   = 0x00800,
    // Topologies - pizeros
    kTop0pizero = 0x01000,
    kTop1pizero = 0x02000,
    kTop2pizero = 0x04000,
    kTopNpizero = 0x08000,
    // Topologies - neutrons
    kTop0neutron = 0x10000,
    kTop1neutron = 0x20000,
    kTop2neutron = 0x40000,
    kTopNneutron = 0x80000
    
  } TopologyType;

}

#endif // CVN_INTERACTIONTYPE_H
