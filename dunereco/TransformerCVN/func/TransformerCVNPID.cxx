#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "dunereco/TransformerCVN/func/TransformerCVNPID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cnn
{

  TransformerCVNPID::TransformerCVNPID(const int pdg, const float val){
    fPdg = pdg;
    fVal = val;
  }

  TransformerCVNPID::TransformerCVNPID():
  fPdg(),fVal()
  {}

}
