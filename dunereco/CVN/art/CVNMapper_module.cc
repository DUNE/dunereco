////////////////////////////////////////////////////////////////////////
// \file    CVNMapper_module.cc
// \brief   Producer module for creating CVN PixelMap objects
// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include "TTree.h"
#include "TH2F.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/CVN/art/PixelMapProducer.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/TrainingData.h"




namespace cvn {

  class CVNMapper : public art::EDProducer {
  public:
    explicit CVNMapper(fhicl::ParameterSet const& pset);
    ~CVNMapper();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    /// Module lablel for input clusters
    std::string    fHitsModuleLabel;

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
    unsigned int fMaxWireGap;

    /// Use unwrapped pixel maps?
    bool fUnwrappedPixelMap;

    /// PixelMapProducer does the work for us
    PixelMapProducer fProducer;

  };



  //.......................................................................
  CVNMapper::CVNMapper(fhicl::ParameterSet const& pset):
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>    ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short> ("MinClusterHits")),
  fTdcWidth     (pset.get<unsigned short> ("TdcWidth")),
  fWireLength   (pset.get<unsigned short> ("WireLength")),
  fTimeResolution   (pset.get<unsigned short> ("TimeResolution")),
  fUnwrappedPixelMap(pset.get<bool> ("UnwrappedPixelMap")),
  fProducer      (fWireLength, fTdcWidth, fTimeResolution)
  {

    produces< std::vector<cvn::PixelMap>   >(fClusterPMLabel);

  }

  //......................................................................
  CVNMapper::~CVNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNMapper::beginJob()
  {  }

  //......................................................................
  void CVNMapper::endJob()
  {
  }

  //......................................................................
  void CVNMapper::produce(art::Event& evt)
  {
    // Use unwrapped pixel maps if requested
    fProducer.SetUnwrapped(fUnwrappedPixelMap);

    art::Handle< std::vector< recob::Hit > > hitListHandle;
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    unsigned short nhits = hitlist.size();

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::PixelMap> >
      pmCol(new std::vector<cvn::PixelMap>);

    if(nhits>fMinClusterHits){
      //std::cout<<"nhits: "<<nhits<<std::endl;
      PixelMap pm = fProducer.CreateMap(hitlist);
      pmCol->push_back(pm);
    }
    //pm.Print();
    //Boundary bound = pm.Bound();
    //}
    evt.put(std::move(pmCol), fClusterPMLabel);
    //std::cout<<"Map Complete!"<<std::endl;
  }

  //----------------------------------------------------------------------



DEFINE_ART_MODULE(cvn::CVNMapper)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







