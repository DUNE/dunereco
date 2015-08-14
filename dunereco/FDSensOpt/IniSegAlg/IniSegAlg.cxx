#include "IniSegAlg.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT libraries
#include "TMath.h"

// C/C++ libraries
#include <memory>

dunefd::Hit2D::Hit2D(TVector2 point2d, size_t key) :
fPoint(point2d),
fKey(key)
{
}

dunefd::IniSegAlg::IniSegAlg() 
{
}

