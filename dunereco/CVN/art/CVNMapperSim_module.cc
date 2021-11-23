////////////////////////////////////////////////////////////////////////
// \file    CVNMapperSim_module.cc
// \brief   Producer module for creating CVN PixelMap objects
// \author  Alexander Radovic - a.radovic@gmail.com
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

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"

#include "dune/CVN/art/PixelMapSimProducer.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/TrainingData.h"




namespace cvn {

  class CVNMapperSim : public art::EDProducer {
  public:
    explicit CVNMapperSim(fhicl::ParameterSet const& pset);
    ~CVNMapperSim();

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
    //unsigned int fMaxWireGap;

    /// Use unwrapped pixel maps?
    // 0 means no unwrap, 1 means unwrap in wire, 2 means unwrap in wire and time
    unsigned short fUnwrappedPixelMap;
    
    /// ADC threshold for calculating charge from wires directly
    double fThreshold;

    /// PixelMapProducer does the work for us
    PixelMapSimProducer fProducer;

  };



  //.......................................................................
  CVNMapperSim::CVNMapperSim(fhicl::ParameterSet const& pset): EDProducer{pset},
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>    ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short> ("MinClusterHits")),
  fTdcWidth     (pset.get<unsigned short> ("TdcWidth")),
  fWireLength   (pset.get<unsigned short> ("WireLength")),
  fTimeResolution   (pset.get<unsigned short> ("TimeResolution")),
  fUnwrappedPixelMap(pset.get<unsigned short> ("UnwrappedPixelMap")),
  fThreshold        (pset.get<double>("Threshold")),
  fProducer      (fWireLength, fTdcWidth, fTimeResolution, fThreshold)
  {

    produces< std::vector<cvn::PixelMap>   >(fClusterPMLabel);

  }

  //......................................................................
  CVNMapperSim::~CVNMapperSim()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNMapperSim::beginJob()
  {  }

  //......................................................................
  void CVNMapperSim::endJob()
  {
  }

  //......................................................................
  void CVNMapperSim::produce(art::Event& evt)
  {
    // Use unwrapped pixel maps if requested
    // 0 means no unwrap, 1 means unwrap in wire, 2 means unwrap in wire and time
    fProducer.SetUnwrapped(fUnwrappedPixelMap);

    std::vector< art::Ptr< sim::SimChannel > > hitlist;
    auto hitListHandle = evt.getHandle< std::vector< sim::SimChannel > >(fHitsModuleLabel);
    if (hitListHandle)
      art::fill_ptr_vector(hitlist, hitListHandle);
    // unsigned short nhits = hitlist.size();

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::PixelMap> >
      pmCol(new std::vector<cvn::PixelMap>);

    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    PixelMap pm = fProducer.CreateMap(detProp, hitlist);
    std::cout << "NROI : " << fProducer.NROI() << std::endl;
    if(fProducer.NROI() > fMinClusterHits)
      pmCol->push_back(pm);
    //pm.Print();
    //Boundary bound = pm.Bound();
    //}
    evt.put(std::move(pmCol), fClusterPMLabel);
    //std::cout<<"Map Complete!"<<std::endl;
  }

  //----------------------------------------------------------------------



DEFINE_ART_MODULE(cvn::CVNMapperSim)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
