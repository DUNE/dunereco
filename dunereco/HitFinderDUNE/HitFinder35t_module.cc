#ifndef HITFINDER35T_H
#define HITFINDER35T_H

////////////////////////////////////////////////////////////////////////
//
// HitFinder35t class
// 
// tjyang@fnal.gov
//  
// Based on Tom Junk's idea of 3-view matching
// HitFinder for 35t, input is undisambiguated hits, out is disambiguated hits
// Start from code written by talion@gmail.com
//
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <utility> // std::move()
#include <memory> // std::unique_ptr<>

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

namespace dune{
  class HitFinder35t : public art::EDProducer {
    
  public:
    
    explicit HitFinder35t(fhicl::ParameterSet const& pset); 
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                
    
    
  private:
    
    DisambigAlg35t    fDisambigAlg;
    TimeBasedDisambig fTimeBasedDisambigAlg;

    art::ServiceHandle<geo::Geometry> fGeom;
    
    std::string fChanHitLabel;
    std::string fAlg;    // which algorithm to use
    
  protected: 
    
    
  }; // class HitFinder35t
  
  
  //-------------------------------------------------
  //-------------------------------------------------
  HitFinder35t::HitFinder35t(fhicl::ParameterSet const& pset) :
    fDisambigAlg(pset.get< fhicl::ParameterSet >("DisambigAlg")),
    fTimeBasedDisambigAlg(pset.get< fhicl::ParameterSet >("TimeBasedDisambigAlg"))
  {
    this->reconfigure(pset);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
  }
  
  
  //-------------------------------------------------
  //-------------------------------------------------
  void HitFinder35t::reconfigure(fhicl::ParameterSet const& p)
  {
    
    fChanHitLabel =  p.get< std::string >("ChanHitLabel");
    fAlg = p.get < std::string >("Algorithm");
    
  }  
  
  //-------------------------------------------------
  //-------------------------------------------------
  void HitFinder35t::beginJob()
  {
    
  }
  
  //-------------------------------------------------
  //-------------------------------------------------
  void HitFinder35t::endJob()
  {
    
  }
  

  //-------------------------------------------------
  void HitFinder35t::produce(art::Event& evt)
  {
    
    art::Handle< std::vector<recob::Hit> > ChannelHits;
    evt.getByLabel(fChanHitLabel, ChannelHits);
    
    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> ChannelHitRawDigits
      (ChannelHits, evt, fChanHitLabel);
    art::FindOneP<recob::Wire> ChannelHitWires
      (ChannelHits, evt, fChanHitLabel);

    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this, evt,
      /* doWireAssns */ ChannelHitWires.isValid(),
      /* doRawDigitAssns */ ChannelHitRawDigits.isValid()
      );
    
    // Make unambiguous collection hits
    std::vector< art::Ptr<recob::Hit> >  ChHits;
    art::fill_ptr_vector(ChHits, ChannelHits);
    for( size_t h = 0; h < ChHits.size(); h++ ){
      if( ChHits[h]->View() != geo::kZ ) continue;
      
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(h);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(h);
    
      // just copy it
      hcol.emplace_back(*ChHits[h], wire, rawdigits);
      
    } // for
    
    // Run alg on all APAs

    if (fAlg == "TripletMatch")
      {
        fDisambigAlg.RunDisambig(ChHits);

	for( size_t t=0; t < fDisambigAlg.fDisambigHits.size(); t++ ){
	  art::Ptr<recob::Hit>  hit = fDisambigAlg.fDisambigHits[t].first;
	  geo::WireID           wid = fDisambigAlg.fDisambigHits[t].second;
      
	  // create a new hit copy of the original one, but with new wire ID
	  recob::HitCreator disambiguous_hit(*hit, wid);
      
	  // get the objects associated with the original hit;
	  // since hit comes from ChannelHits, its key is the index in that collection
	  // and also the index for the query of associated objects
	  art::Ptr<recob::Hit>::key_type hit_index = hit.key();
	  art::Ptr<recob::Wire> wire = ChannelHitWires.at(hit_index);
	  art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(hit_index);
      
	  hcol.emplace_back(disambiguous_hit.move(), wire, rawdigits);
	} // for

      }
    else if (fAlg == "TimeBased")
      {
        fTimeBasedDisambigAlg.RunDisambig(ChHits);

	for( size_t t=0; t < fTimeBasedDisambigAlg.fDisambigHits.size(); t++ ){
	  art::Ptr<recob::Hit>  hit = fTimeBasedDisambigAlg.fDisambigHits[t].first;
	  geo::WireID           wid = fTimeBasedDisambigAlg.fDisambigHits[t].second;
      
	  // create a new hit copy of the original one, but with new wire ID
	  recob::HitCreator disambiguous_hit(*hit, wid);
      
	  // get the objects associated with the original hit;
	  // since hit comes from ChannelHits, its key is the index in that collection
	  // and also the index for the query of associated objects
	  art::Ptr<recob::Hit>::key_type hit_index = hit.key();
	  art::Ptr<recob::Wire> wire = ChannelHitWires.at(hit_index);
	  art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(hit_index);
      
	  hcol.emplace_back(disambiguous_hit.move(), wire, rawdigits);
	} // for

      }
    else
      {
        throw cet::exception("HitFinder35t") << "Disambiguation algorithm name: " << fAlg << " is not supported.\n" ;
      }
    
    
    // put the hit collection and associations into the event
    hcol.put_into(evt);
    
  }
  
  
  DEFINE_ART_MODULE(HitFinder35t)
  
} // end of dune namespace
#endif // HITFINDER35T_H
