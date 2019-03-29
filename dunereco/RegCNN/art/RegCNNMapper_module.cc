////////////////////////////////////////////////////////////////////////
// \file    RegCNNMapper_module.cc
// \brief   Producer module for creating RegCNN PixelMap objects modified from CVNMapper_module.cc
// \author  Ilsoo Seong - iseong@uci.edu
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
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/RegCNN/art/RegPixelMapProducer.h"
#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/func/RegCNNResult.h"


namespace cnn {

  class RegCNNMapper : public art::EDProducer {
  public:
    explicit RegCNNMapper(fhicl::ParameterSet const& pset);
    ~RegCNNMapper();

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

    /// select which global wire method
    int fGlobalWireMethod;
    /// select how to choose center of pixel map
    int fUseRecoVertex;
    std::string fRegCNNResultLabel;
    std::string fRegCNNModuleLabel;

    /// PixelMapProducer does the work for us
    RegPixelMapProducer fProducer;

  };



  //.......................................................................
  RegCNNMapper::RegCNNMapper(fhicl::ParameterSet const& pset):
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel   (pset.get<std::string>       ("ClusterPMLabel")),
  fMinClusterHits   (pset.get<unsigned short>    ("MinClusterHits")),
  fTdcWidth         (pset.get<unsigned short>     ("TdcWidth")),
  fWireLength       (pset.get<unsigned short>     ("WireLength")),
  fTimeResolution   (pset.get<unsigned short>     ("TimeResolution")),
  fGlobalWireMethod (pset.get<int>                ("GlobalWireMethod")),
  fUseRecoVertex    (pset.get<int>                ("UseRecoVertex")),
  fRegCNNResultLabel (pset.get<std::string>       ("RegCNNResultLabel")),
  fRegCNNModuleLabel (pset.get<std::string>       ("RegCNNModuleLabel")),
  fProducer      (fWireLength, fTdcWidth, fTimeResolution, fGlobalWireMethod)
  {

    produces< std::vector<cnn::RegPixelMap>   >(fClusterPMLabel);

  }

  //......................................................................
  RegCNNMapper::~RegCNNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void RegCNNMapper::beginJob()
  {  }

  //......................................................................
  void RegCNNMapper::endJob()
  {
  }

  //......................................................................
  void RegCNNMapper::produce(art::Event& evt)
  {

    art::Handle< std::vector< recob::Hit > > hitListHandle;
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    art::FindManyP<recob::Wire> fmwire(hitListHandle, evt, fHitsModuleLabel);

    unsigned short nhits = hitlist.size();

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cnn::RegPixelMap> >
      pmCol(new std::vector<cnn::RegPixelMap>);

    if(nhits>fMinClusterHits){
      RegPixelMap pm;
      if (!(fUseRecoVertex)){
          // create pixel map based on mean of wire and ticks
          pm = fProducer.CreateMap(hitlist, fmwire);
      } else{
        // create pixel map based on the reconstructed vertex
        // Get RegCNN Results
        art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle;
        evt.getByLabel(fRegCNNModuleLabel, fRegCNNResultLabel, cnnresultListHandle);
        float vtx[3] = {-99999, -99999, -99999};
        if (!cnnresultListHandle.failedToGet())
        {
            if (!cnnresultListHandle->empty())
            {
                const std::vector<float>& v = (*cnnresultListHandle)[0].fOutput;
                for (unsigned int ii = 0; ii < 3; ii++){
                     vtx[ii] = v[ii]; 
                }
            }
        }
        pm = fProducer.CreateMap(hitlist, fmwire, vtx);
      }
      // skip if PixelMap is empty
      if (pm.fInPM) pmCol->push_back(pm);
      //pm.Print();
    }
    //Boundary bound = pm.Bound();
    //}
    evt.put(std::move(pmCol), fClusterPMLabel);
    //std::cout<<"Map Complete!"<<std::endl;
  }

  //----------------------------------------------------------------------



DEFINE_ART_MODULE(cnn::RegCNNMapper)
} // end namespace cnn
////////////////////////////////////////////////////////////////////////







