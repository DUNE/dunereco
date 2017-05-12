////////////////////////////////////////////////////////////////////////
/// \file    Result.h
/// \brief   Result for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include  <iostream>
#include  <ostream>

#include "dune/CVN/func/Result.h"

namespace cvn
{

  Result::Result(const float* output, unsigned int& nOutputs):
  fOutput(nOutputs)
  {
    for(size_t i = 0; i < nOutputs; ++i) fOutput[i] = output[i];
  }

  Result::Result():
  fOutput()
  {}


}
