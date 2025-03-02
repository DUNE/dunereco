////////////////////////////////////////////////////////////////////////
/// \file    TransformerCVNPID.h
/// \brief   TransformerCVNPID for TransformerCVN prong predictions
/// \author  Alejandro Yankelevich - ayankele@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef TRANSFORMERCVN_PID_H
#define TRANSFORMERCVN_PID_H

#include <vector>

namespace cnn
{
  /// TransformerCVNPID, output of CNN neural net for PID
  class TransformerCVNPID
  {
  public:
    TransformerCVNPID(const int pdg, const float val);
    TransformerCVNPID();

    int Pdg()      const   {return fPdg; }
    float Value()  const   {return fVal;}

  protected:
    int    fPdg;  ///< pdg code 
    float  fVal;  ///< network output pid value

  };
}

#endif  // TRANSFORMERCVN_PID_H

