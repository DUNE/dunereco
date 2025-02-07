#include "dunereco/TransformerCVN/func/TransformerCVNPType.h"

namespace cnn
{
  int GetPDGByPType(PType ptype)
  {

    switch (ptype) {
    case  kElectron:
      return 11;
    case  kMuon:
      return  13;
    case  kProton:
      return  2212;
    case  kNeutron:
      return  2112;
    case  kPion:
      return  211;
    case  kPiZero:
      return  111;
    case  kGamma:
      return  22;
      
      // In the case of composite CVN scores like EMID and HadronID,
      // represent the PDG code of those types as the product of the two 
      // pdg codes. 
      // Eg. pdgEM = 22 * 11 = 242

      // This way we can make comparisons with %
      //     pdgEM % 11 == 0
    case kEM:
      // pdg 22*11 = 242 describes a bound state of an up and charm quark with spin 1/2
      // This breaks quantum mechanics, so we're safe
      return 22 * 11;
    case kHadron:
      // pdg 2212*211 = 466732 describes a tb's baryon in the 4th radial eigen state
      // and an unphysical orbital momentum state.
      // Pretty sure we won't have to deal with those
      return 2212 * 211;
    default:
      return 0;
    }

  }// GetPDGByPType
}
