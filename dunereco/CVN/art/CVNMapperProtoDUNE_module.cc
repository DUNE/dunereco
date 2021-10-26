////////////////////////////////////////////////////////////////////////
// \file    CVNMapperProtoDUNE_module.cc
// \brief   Producer module for creating CVN PixelMaps for individual
//          particles in ProtoDUNE
// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "dune/CVN/art/PixelMapProducer.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/TrainingData.h"

#include "dune/CVN/func/CVNProtoDUNEUtils.h"

namespace cvn {

  class CVNMapperProtoDUNE : public art::EDProducer {
  public:
    explicit CVNMapperProtoDUNE(fhicl::ParameterSet const& pset);
    ~CVNMapperProtoDUNE();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    /// Module label for input hits
    std::string    fHitsModuleLabel;

    /// Module label for input particles
    std::string    fParticleModuleLabel;

    /// Module label for input tracks
    std::string    fTrackLabel;

    /// Module label for input showers
    std::string    fShowerLabel;

    /// Instance lablel for cluster pixelmaps
    std::string    fClusterPMLabel;

    /// Minimum number of hits for cluster to be converted to pixel map
    unsigned short fMinClusterHits;

    /// Width of pixel map in tdcs
    unsigned short fTdcWidth;

    /// Length of pixel map in wires
    unsigned short fWireLength;

    /// Length of pixel map in wires
    double fTimeResolution;

    /// Maximum gap in wires at front of cluster to prevent pruning of upstream
    /// hits
    //unsigned int fMaxWireGap;

    /// Use unwrapped pixel maps?
    // 0 means no unwrap, 1 means unwrap in wire, 2 means unwrap in wire and time
    unsigned short fUnwrappedPixelMap;

    /// Track length cut
    float fTrackLengthCut;

    /// Use the whole event instead of each track and shower
    bool fUseWholeEvent;

    /// For protoDUNE vertex finding, we only want the beam slice
    bool fUseBeamSliceOnly;

    /// PixelMapProducer does the work for us
    PixelMapProducer fProducer;

  };



  //.......................................................................
  CVNMapperProtoDUNE::CVNMapperProtoDUNE(fhicl::ParameterSet const& pset): EDProducer{pset},
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fParticleModuleLabel  (pset.get<std::string>    ("ParticleModuleLabel")),
  fTrackLabel  (pset.get<std::string>    ("TrackLabel")),
  fShowerLabel  (pset.get<std::string>    ("ShowerLabel")),
  fClusterPMLabel(pset.get<std::string>    ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short> ("MinClusterHits")),
  fTdcWidth     (pset.get<unsigned short> ("TdcWidth")),
  fWireLength   (pset.get<unsigned short> ("WireLength")),
  fTimeResolution   (pset.get<unsigned short> ("TimeResolution")),
  fUnwrappedPixelMap(pset.get<unsigned short> ("UnwrappedPixelMap")),
  fTrackLengthCut(pset.get<unsigned short> ("TrackLengthCut")),
  fUseWholeEvent(pset.get<bool> ("UseWholeEvent")),
  fUseBeamSliceOnly(pset.get<bool> ("UseBeamSliceOnly")),
  fProducer      (fWireLength, fTdcWidth, fTimeResolution)
  {

    produces< std::vector<cvn::PixelMap>   >(fClusterPMLabel);

  }

  //......................................................................
  CVNMapperProtoDUNE::~CVNMapperProtoDUNE()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNMapperProtoDUNE::beginJob()
  {  }

  //......................................................................
  void CVNMapperProtoDUNE::endJob()
  {
  }

  //......................................................................
  void CVNMapperProtoDUNE::produce(art::Event& evt)
  {

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::PixelMap> > pmCol(new std::vector<cvn::PixelMap>);

    // Use unwrapped pixel maps if requested
    // For protoDUNE unwrapped if > 0
    fProducer.SetUnwrapped(fUnwrappedPixelMap);
    fProducer.SetProtoDUNE();

    // Use the whole event just like we would for the FD
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    if(fUseWholeEvent){
      std::vector< art::Ptr< recob::Hit > > hitlist;
      auto hitListHandle = evt.getHandle< std::vector< recob::Hit > >(fHitsModuleLabel);
      if (hitListHandle)
        art::fill_ptr_vector(hitlist, hitListHandle);
      unsigned short nhits = hitlist.size();
  
      std::cout << "nhits: " << nhits << std::endl; // REMOVE 
 
      if(nhits>fMinClusterHits){
        PixelMap pm = fProducer.CreateMap(detProp, hitlist);
        pmCol->push_back(pm);
      }
    }
    else if(fUseBeamSliceOnly){
      // We want to make a pixel map for just APA3. We can pad out to 500 pixels in the wire number
      // The best way to do this is to get the Pandora beam slice
      cvn::CVNProtoDUNEUtils pfpUtil;
      const unsigned int beamSlice = pfpUtil.GetBeamSlice(evt,fParticleModuleLabel);
      if(beamSlice < 500){
        const std::vector<const recob::Hit*> sliceHits = pfpUtil.GetRecoSliceHits(beamSlice,evt,fParticleModuleLabel);

        PixelMap pm = fProducer.CreateMap(detProp, sliceHits);
        pmCol->push_back(pm);
      }
    }
    else{
      // Get the list of tracks and showers and their associated hits
      auto allTracks=evt.getValidHandle<std::vector<recob::Track> >(fTrackLabel);
      const art::FindManyP<recob::Hit> findTrackHits(allTracks,evt,fTrackLabel);
  
      auto allShowers=evt.getValidHandle<std::vector<recob::Shower> >(fShowerLabel);
      const art::FindManyP<recob::Hit> findShowerHits(allShowers,evt,fShowerLabel);
  
      std::cout << "Event contains " << allTracks->size() << " tracks and " << allShowers->size() << " showers" << std::endl;
  
      // Iterate over the tracks and showers making a PixelMap for each 
      for(unsigned int t = 0; t < allTracks->size(); ++t){
  
        std::vector<art::Ptr<recob::Hit> > trackHits = findTrackHits.at(t);
        if(trackHits.size()>fMinClusterHits && (*allTracks)[t].Length() > fTrackLengthCut){
  //        std::cout << "Making PixelMap number " << pmCol->size() + 1 << " with " << trackHits.size() << " hits" << std::endl; 
  
          PixelMap pm = fProducer.CreateMap(detProp, trackHits);
          pmCol->push_back(pm);
        }
  
      }
  
      for(unsigned int s = 0; s < allShowers->size(); ++s){
  
        std::vector<art::Ptr<recob::Hit> > showerHits = findShowerHits.at(s);
        if(showerHits.size()>fMinClusterHits){
          PixelMap pm = fProducer.CreateMap(detProp, showerHits);
          pmCol->push_back(pm);
        }
      }
    }

    std::cout << "Made " << pmCol->size() << " pixel maps for this event" << std::endl;
    
    evt.put(std::move(pmCol), fClusterPMLabel);

  }

  //----------------------------------------------------------------------



DEFINE_ART_MODULE(cvn::CVNMapperProtoDUNE)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
