/////////////////////////////////////////////////////////////////
//
//  Disambiguation algorithm for the protoDUNE single phase detector
//  using two views and the third for confirmation if present.  Does not allow for hits to be on
//  wires on the outside sides of the APA's.
//  Start from code written by talion@gmail.com
//  trj@fnal.gov
//  tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////
#ifndef DisambigAlgProtoDUNESP_H
#define DisambigAlgProtoDUNESP_H
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
#include "larreco/RecoAlg/DBScanAlg.h"

#include "TMatrixD.h"
#include "TVectorD.h"



namespace dune{

  //--------------------------------------------------------------- 
  class DisambigAlgProtoDUNESP {
  public:
    
    
    DisambigAlgProtoDUNESP(fhicl::ParameterSet const& pset);
    
    void               reconfigure(fhicl::ParameterSet const& p);

    void               RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits );
                                                                  ///< Run disambiguation as currently configured


    std::vector< std::pair<art::Ptr<recob::Hit>, geo::WireID> > fDisambigHits;
                                                                   ///< The final list of hits to pass back to be made

    private:

    bool notOuterWire(geo::WireID wireid);

    int longTPC(int apa);

    // other classes we will use
    apa::APAGeometryAlg  fAPAGeo;
    double fTimeCut;
    double fDistanceCut;
    double fDcut2;
  }; // class DisambigAlgProtoDUNESP

} // namespace dune

#endif // ifndef DisambigAlgProtoDUNESP_H
