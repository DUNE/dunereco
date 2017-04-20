#ifndef HITFINDERPROTODUNESP_H
#define HITFINDERPROTODUNESP_H

////////////////////////////////////////////////////////////////////////
//
// HitFinderProtoDUNESP class
// 
// trj@fnal.gov
//  
// Runs the disambiguation algorithm for the single-phase ProtoDUNE detector
// 
// 
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
#include "DisambigAlgProtoDUNESP.h"

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
  class HitFinderProtoDUNESP : public art::EDProducer {
    
  public:
    
    explicit HitFinderProtoDUNESP(fhicl::ParameterSet const& pset); 
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                
    
    
  private:
    
    DisambigAlgProtoDUNESP    fDisambigAlg;

    art::ServiceHandle<geo::Geometry> fGeom;
    
    std::string fChanHitLabel;
    std::string fAlg;    // which algorithm to use
    
  protected: 
    
    
  }; // class HitFinderProtoDUNESP
  
  
  //-------------------------------------------------
  //-------------------------------------------------
  HitFinderProtoDUNESP::HitFinderProtoDUNESP(fhicl::ParameterSet const& pset) :
    fDisambigAlg(pset.get< fhicl::ParameterSet >("DisambigAlg"))
  {
    this->reconfigure(pset);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
  }
  
  
  //-------------------------------------------------
  //-------------------------------------------------
  void HitFinderProtoDUNESP::reconfigure(fhicl::ParameterSet const& p)
  {
    
    fChanHitLabel =  p.get< std::string >("ChanHitLabel");
    fAlg = p.get < std::string >("Algorithm");  // switch on which algorithm to run
    
  }  
  
  //-------------------------------------------------
  //-------------------------------------------------
  void HitFinderProtoDUNESP::beginJob()
  {
    
  }
  
  //-------------------------------------------------
  //-------------------------------------------------
  void HitFinderProtoDUNESP::endJob()
  {
    
  }
  

  //-------------------------------------------------
  void HitFinderProtoDUNESP::produce(art::Event& evt)
  {
    
    //art::Handle< std::vector<recob::Hit> > ChannelHits;
    //evt.getByLabel(fChanHitLabel, ChannelHits);
    
    auto ChannelHits = evt.getValidHandle<std::vector<recob::Hit>>(fChanHitLabel);

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
    
    // Run alg on all APAs -- check fAlg if we have more than one algorithm defined (future expansion)

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

    
  // put the hit collection and associations into the event
  hcol.put_into(evt);

  }
    
  
DEFINE_ART_MODULE(HitFinderProtoDUNESP)
  
} // end of dune namespace
#endif // HITFINDERPROTODUNESP_H
