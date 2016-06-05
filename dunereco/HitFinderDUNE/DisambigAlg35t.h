/////////////////////////////////////////////////////////////////
//
//  Disambiguation algorithm based on Tom Junk's idea of 3-view matching
//  Start from code written by talion@gmail.com
//  trj@fnal.gov
//  tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////
#ifndef DisambigAlg35t_H
#define DisambigAlg35t_H
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
  class DisambigAlg35t {
  public:
    
    
    DisambigAlg35t(fhicl::ParameterSet const& pset);
    
    void               reconfigure(fhicl::ParameterSet const& p);

    void               RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits );
                                                                  ///< Run disambiguation as currently configured


    std::vector< std::pair<art::Ptr<recob::Hit>, geo::WireID> > fDisambigHits;
                                                                   ///< The final list of hits to pass back to be made

    private:

    // other classes we will use
    apa::APAGeometryAlg                           fAPAGeo;
    double fTimeCut;
    double fDistanceCut;
    double fDistanceCutClu;
    double fTimeWiggle;
    int    fColChanWiggle;
    bool   fDoCleanUpHits;
    cluster::DBScanAlg fDBScan; ///< object that implements the DB scan algorithm
  }; // class DisambigAlg35t

} // namespace dune

#endif // ifndef DisambigAlg35t_H
