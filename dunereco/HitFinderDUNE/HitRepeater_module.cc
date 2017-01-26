#ifndef HITREPEATER__H
#define HITREPEATER__H

////////////////////////////////////////////////////////////////////////
//
// HitRepeater class
// 
// pplonski@ire.pw.edu.pl
//  
// This class repeats hits on wires if channel has more than one wire.
// Based on Robert Sulej's idea - we don't need a disambiguation step 
// before reconstruction; disambiguation step can be done during 
// 3D tracks reconstruction.
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <utility> 
#include <memory>  
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "DisambigAlg35t.h"
#include "TimeBasedDisambig.h"

// ROOT Includes
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVector2.h"
#include "TVector3.h"

namespace dune {
  class HitRepeater : public art::EDProducer {
  public:
    explicit HitRepeater(fhicl::ParameterSet const& pset);
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt);
    void beginJob() {}
    void endJob() {}

  private:
    art::ServiceHandle<geo::Geometry> geom;
    std::string fChanHitLabel;
  };

  // implementation

  HitRepeater::HitRepeater(fhicl::ParameterSet const& pset) {
    this->reconfigure(pset);
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
  }

  void HitRepeater::reconfigure(fhicl::ParameterSet const& p) {
    fChanHitLabel =  p.get< std::string >("ChanHitLabel");
  }

  void HitRepeater::produce(art::Event& evt) {
    art::Handle< std::vector<recob::Hit> > ChannelHits;
    evt.getByLabel(fChanHitLabel, ChannelHits);

    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> ChannelHitRawDigits(ChannelHits, evt, fChanHitLabel);
    art::FindOneP<recob::Wire>   ChannelHitWires    (ChannelHits, evt, fChanHitLabel);

    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this, evt,
       ChannelHitWires.isValid(), // doWireAssns
       ChannelHitRawDigits.isValid() // doRawDigitAssns
      );
    
    // Make hits collection
    std::vector< art::Ptr<recob::Hit> >  ChHits;
    art::fill_ptr_vector(ChHits, ChannelHits);
      
    for( size_t h = 0; h < ChHits.size(); h++ ) {
      std::vector<geo::WireID> cwids = geom->ChannelToWire(ChHits[h]->Channel());
      for(size_t w = 0; w < cwids.size(); w++) {
        art::Ptr<recob::Hit>  hit = ChHits[h];
	      geo::WireID           wid = cwids[w];
        recob::HitCreator repeated_hit(*hit, wid);
        art::Ptr<recob::Hit>::key_type hit_index = hit.key();
	      art::Ptr<recob::Wire> wire = ChannelHitWires.at(hit_index);
	      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(hit_index);

	      hcol.emplace_back(repeated_hit.move(), wire, rawdigits);
      }
    }

    // put the hit collection and associations into the event
    hcol.put_into(evt);    
  }

  DEFINE_ART_MODULE(HitRepeater)

} // end of dune namespace
#endif 

