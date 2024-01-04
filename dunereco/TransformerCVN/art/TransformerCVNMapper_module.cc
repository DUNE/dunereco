////////////////////////////////////////////////////////////////////////
// \file    TransformerCVNMapper_module.cc
// \brief   Producer module for creating TransformerCVN PixelMap objects modified from TransformerCVNMapper_module.cc
// \author  Alejandro Yankelevich - ayankele@uci.edu
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
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dunereco/TransformerCVN/art/TransformerPixelMapProducer.h"
#include "dunereco/TransformerCVN/func/TransformerPixelMap.h"
//#include "dunereco/TransformerCVN/func/TransformerPixelMap3D.h"
#include "dunereco/RegCNN/func/RegCNNResult.h"

namespace lar_pandora{class LArPandoraHelper;}

namespace cnn {

  class TransformerCVNMapper : public art::EDProducer {
  public:
    explicit TransformerCVNMapper(fhicl::ParameterSet const& pset);
    ~TransformerCVNMapper();

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

    /// parameters for 2D pixel maps
    /// Width of pixel map in tdcs
    unsigned short fTdcWidth;
    /// Length of pixel map in wires
    unsigned short fWireLength;
    // resolution of ticks
    double fTimeResolution;
    // resolution of wires
    double fWireResolution;

    /// Maximum gap in wires at front of cluster to prevent pruning of upstream
    /// hits
    //unsigned int fMaxWireGap;

    /// Use unwrapped pixel maps?
    //dla bool fUnwrappedPixelMap;

    /// select which global wire method
    int fGlobalWireMethod;
    /// select how to choose center of pixel map
    int fUseRecoVertex;
    /// select how to tag prong
    int fProngTagMethod;
    // Use Prong only pixel maps
    bool fProngOnly;
    // generate pixel map by raw charge or by reco. hit
    bool fByHit;

    std::string fShowerModuleLabel;
    std::string fTrackModuleLabel;
    std::string fVertexModuleLabel;
    std::string fPFParticleModuleLabel;
    std::string fPandoraNuVertexModuleLabel;

    std::string fRegCNNResultLabel;
    std::string fRegCNNModuleLabel;

    TransformerPixelMapProducer fProducer;
  };
  
   //.......................................................................
  TransformerCVNMapper::TransformerCVNMapper(fhicl::ParameterSet const& pset):
  EDProducer(pset),
  fHitsModuleLabel  (pset.get<std::string>         ("HitsModuleLabel")),
  fClusterPMLabel   (pset.get<std::string>         ("ClusterPMLabel")),
  fMinClusterHits   (pset.get<unsigned short>      ("MinClusterHits")),
  fTdcWidth         (pset.get<unsigned short>      ("TdcWidth")),
  fWireLength       (pset.get<unsigned short>      ("WireLength")),
  fTimeResolution   (pset.get<unsigned short>      ("TimeResolution")),
  fWireResolution   (pset.get<unsigned short>      ("WireResolution")),
  fGlobalWireMethod (pset.get<int>                 ("GlobalWireMethod")),
  fUseRecoVertex    (pset.get<int>                 ("UseRecoVertex")),
  fProngTagMethod   (pset.get<int>                 ("ProngTagMethod")),
  fProngOnly        (pset.get<bool>                ("ProngOnly")),
  fByHit            (pset.get<bool>                ("ByHit")),
  fShowerModuleLabel(pset.get<std::string>         ("ShowerModuleLabel")),
  fTrackModuleLabel (pset.get<std::string>         ("TrackModuleLabel")),
  fVertexModuleLabel(pset.get<std::string>         ("VertexModuleLabel")),
  fPFParticleModuleLabel(pset.get<std::string>     ("PFParticleModuleLabel")),
  fPandoraNuVertexModuleLabel(pset.get<std::string>("PandoraNuVertexModuleLabel")),
  fRegCNNResultLabel (pset.get<std::string>        ("RegCNNResultLabel")),
  fRegCNNModuleLabel (pset.get<std::string>        ("RegCNNModuleLabel")),
  fProducer(TransformerPixelMapProducer(fWireLength, fWireResolution, fTdcWidth, fTimeResolution, fGlobalWireMethod, fProngOnly, fByHit))
  { 
      produces< std::vector<cnn::TransformerPixelMap> >(fClusterPMLabel);
  }
  
  //......................................................................
  TransformerCVNMapper::~TransformerCVNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }
  
  //......................................................................
  void TransformerCVNMapper::beginJob()
  {  
  }

  //......................................................................
  void TransformerCVNMapper::endJob()
  {
  }

  void TransformerCVNMapper::produce(art::Event& evt)
  {
    std::cout<<"ievet : "<<evt.id()<<std::endl;

    std::vector< art::Ptr< recob::Hit > > hitlist;
    auto hitListHandle = evt.getHandle< std::vector< recob::Hit > >(fHitsModuleLabel);
    if (hitListHandle)
      art::fill_ptr_vector(hitlist, hitListHandle);

    art::FindManyP<recob::Wire> fmwire(hitListHandle, evt, fHitsModuleLabel);
    art::FindManyP<recob::Shower> fmshwhit(hitListHandle, evt, fShowerModuleLabel);
    art::FindManyP<recob::Track> fmtrkhit(hitListHandle, evt, fTrackModuleLabel);

    unsigned short nhits = hitlist.size();

    // Get the vertex out of the event record
    std::vector<art::Ptr<recob::Vertex> > vertex_list;
    auto vertexHandle = evt.getHandle<std::vector<recob::Vertex> >(fVertexModuleLabel);
    if (vertexHandle)
        art::fill_ptr_vector(vertex_list, vertexHandle);
    art::FindMany<recob::PFParticle> fmPFParticle(vertexHandle, evt, fPFParticleModuleLabel);

    art::FindManyP<recob::SpacePoint> fmSPFromHits(hitListHandle, evt, fPFParticleModuleLabel);

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cnn::TransformerPixelMap> >
      pmCol(new std::vector<cnn::TransformerPixelMap>);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    if(nhits>fMinClusterHits){
            TransformerPixelMap pm;
            if (fUseRecoVertex==0){
                // create pixel map based on mean of wire and ticks
                pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire);
            } 
            else if (fUseRecoVertex==1) {
                // create pixel map based on the reconstructed vertex
                // Get RegCNN Results
	        art::InputTag itag1(fRegCNNModuleLabel, fRegCNNResultLabel);
	        auto cnnresultListHandle = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag1);
                std::vector<float> vtx(3, -99999);
                if (!cnnresultListHandle.failedToGet())
                {
                    if (!cnnresultListHandle->empty())
                    {
                        const std::vector<float>& v = (*cnnresultListHandle)[0].fOutput;
                        for (unsigned int ii = 0; ii < 3; ii++){
                            vtx[ii] = v[ii];
                            std::cout<<"vertex "<<ii<<" : "<<vtx[ii]<<std::endl;
                        }
                    }
                }
                if (fProngTagMethod==1) {
                    pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, fmtrkhit, vtx);
                } else {
                    pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, fmshwhit, vtx);
                }
            }
            else if (fUseRecoVertex==2) {
                // create pixel map based on pandora vertex
                lar_pandora::PFParticleVector particleVector;
                lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
                lar_pandora::VertexVector vertexVector;
                lar_pandora::PFParticlesToVertices particlesToVertices;
                lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);

                int n_prims = 0;
                int n_pand_particles = particleVector.size();
                //std::cout<<"n_pand_particles "<<n_pand_particles<<std::endl;
                double xyz_temp[3] = {0.0, 0.0, 0.0};
                std::vector<float> vtx(3, -99999);
                for (int ipfp = 0; ipfp< n_pand_particles; ipfp++) {
                    const art::Ptr<recob::PFParticle> particle = particleVector.at(ipfp);
                    if (!particle->IsPrimary()) continue;
                    n_prims++;
                    // Particles <-> Vertices
                    lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
                    if (particlesToVertices.end() != vIter) {
                        const lar_pandora::VertexVector& vertexVector = vIter->second;
                        if (!vertexVector.empty()) {
                            if (vertexVector.size()!=1)
                                std::cout<<" Warning: Found particle with more than one associated vertex "<<std::endl;

                            const art::Ptr<recob::Vertex> vertex_pfp = *(vertexVector.begin());
                            vertex_pfp->XYZ(xyz_temp);
                            for (unsigned int ii=0; ii<3; ++ii) {
                                vtx[ii] = xyz_temp[ii];
                            }
                            //std::cout<<"Pandora vertex "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<std::endl;
                        }
                    }
                } 
                if (fProngTagMethod==1) {
                    pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, fmtrkhit, vtx);
                } else {
                    std::cout<<"FIXME: CreateMap from fmshwhit and pandora vtx"<<std::endl;
                    pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, fmshwhit, vtx);
                }
            }
            else {
                // create pixel map based on the reconstructed vertex on wire and tick coordinate
                // Get RegCNN Results
	        art::InputTag itag2(fRegCNNModuleLabel, fRegCNNResultLabel);
                auto cnnresultListHandle = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag2);
                std::vector<float> vtx(6, -99999);
                if (!cnnresultListHandle.failedToGet())
                {
                    if (!cnnresultListHandle->empty())
                    {
                        const std::vector<float>& v = (*cnnresultListHandle)[0].fOutput;
                        for (unsigned int ii = 0; ii < 6; ii++){
                            vtx[ii] = v[ii];
                        }
                    }
                }
                if (fProngTagMethod==1) {
                    pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, fmtrkhit, vtx);
                } else {
                    pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, fmshwhit, vtx);
                }
            }
            // skip if PixelMap is empty
            if (pm.fInPM) pmCol->push_back(pm);
            pm.Print();
    }
    
    evt.put(std::move(pmCol), fClusterPMLabel);
    
    std::cout<<"Map Complete!"<<std::endl;
}

  //----------------------------------------------------------------------

DEFINE_ART_MODULE(cnn::TransformerCVNMapper)
} // end namespace cnn
////////////////////////////////////////////////////////////////////////

