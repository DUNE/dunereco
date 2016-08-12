/////////////////////////////////////////////////////////////////
//
//  Disambiguation algorithm based on Tom Junk's idea of 3-view matching
//  Start from code written by talion@gmail.com
//  trj@fnal.gov
//  tjyang@fnal.gov
//  iamjaejunkim@gmail.com
//
////////////////////////////////////////////////////////////////////
#ifndef TimeBasedDisambig_H
#define TimeBasedDisambig_H
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <stdint.h>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Event.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "APAGeometryAlg.h"

#include "TMatrixD.h"
#include "TVectorD.h"


//hit positions need to be saved as a histogram to apply position-correction function but since
//the defined hit object does not have elements to save position informations, induction plane
//hits with their corresponding wireid and returned positions need to be saved in a separate set
//of variables thus defined a struct of a vector to save all predefined variables in the hit object
//with returned and simulated positions and wireid

//define a struct of a vector to save hit information including simulated hit information
struct HitPos{
  unsigned int Channel;
  unsigned int StartTick;
  unsigned int EndTick;
  double PeakTime;
  double SigmaPeakTime;
  double RMS;
  double PeakAmplitude;
  double SigmaPeakAmplitude;
  double SummedADC;
  double Integral;
  double SigmaIntegral;
  unsigned int Multiplicity;
  unsigned int LocalIndex;
  double GoodnessOfFit;
  int DegreesOfFreedom;
  geo::View_t View;
  geo::SigType_t SignalType;
  geo::WireID WireID;
  double FinalYPos;
  double FinalZPos;
  unsigned int TPC;
  unsigned int DisambigWireID;
  //double SimYPos;
  //double SimZPos;
  //double SimWireID;
};


namespace dune{

  //--------------------------------------------------------------- 
  class TimeBasedDisambig {
  public:
    
    
    TimeBasedDisambig(fhicl::ParameterSet const& pset);
    
    void               reconfigure(fhicl::ParameterSet const& p);

    void               RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits );
    //void               RunDisambig();
                                                                  ///< Run disambiguation as currently configured


    std::vector< std::pair<art::Ptr<recob::Hit>, geo::WireID> > fDisambigHits;
                                                                   ///< The final list of hits to pass back to be made

    



    private:

    // other classes we will use
    apa::APAGeometryAlg                           fAPAGeo;
    double fTimeCut;
    double fDistanceCut;
    double fDistanceCutClu;
  }; // class TimeBasedDisambig

} // namespace dune

#endif // ifndef TimeBasedDisambig_H
