////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       RegCNNVertexAna
// Module Type: analyzer
// File:        RegCNNVertexAna_module.cc
// Author:      
////////////////////////////////////////////////////////////////////////////////////////////////

#include <string> 
#include <vector>
#include "TTree.h"
#include "TBranch.h"

// Framework includes: 
#include "art/Framework/Core/ModuleMacros.h"  
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/RegCNN/func/RegCNNResult.h"

const int kMax = 1000;
namespace myana {
  class RegCNNVertexAna;
}

class myana::RegCNNVertexAna : public art::EDAnalyzer
{
  public:
    RegCNNVertexAna(fhicl::ParameterSet const& pset);
    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& p);
    void reset();
  
  private:
    TTree* fTree; 
    int ievt;
    double trueEnergy;
    float  regcnn_energy;
    float  CNN_vertex[3];
    float  CNN_vertex1ststep[3];
    float  CNN_vertex2ndstep[3];
    float  CNN_vertex2[3];
    
    int    InDet;
    int    FidCut;
    int    nhits;
    int    nu_truth_N;
    int    nupdg_truth[kMax];
    int    numode_truth[kMax];    
    int    nuccnc_truth[kMax];     
    double nueng_truth[kMax];
    double nuvtxx_truth[kMax];
    double nuvtxy_truth[kMax];
    double nuvtxz_truth[kMax];

    double ErecoNue;
    double RecoLepEnNue; 
    double RecoHadEnNue; 
    int    RecoMethodNue; 


    int pandora_nu_vtx;
    double pandora_nu_vtx_x;
    double pandora_nu_vtx_y;
    double pandora_nu_vtx_z;


    std::string fMCGenModuleLabel;   
    std::string fHitsModuleLabel;
    std::string fEnergyRecoNueLabel;

    std::string fPandoraNuVertexModuleLabel; 

    std::string fRegCNNResultLabel;
    std::string fRegCNNModuleLabel;

    std::string fRegCNNVtxResultLabel;
    std::string fRegCNNVtx1ststepResultLabel;
    std::string fRegCNNVtx2ndstepResultLabel;
    std::string fRegCNNVtxResult2Label;

    std::string fRegCNNVtxModuleLabel;
    std::string fRegCNNVtx1ststepModuleLabel;
    std::string fRegCNNVtx2ndstepModuleLabel;
    std::string fRegCNNVtxModule2Label;


    art::ServiceHandle<art::TFileService> tfs;

};

myana::RegCNNVertexAna::RegCNNVertexAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
  fTree = tfs->make<TTree>("anatree", "anatree");
  fTree->Branch("ievent",            &ievt,              "ievent/I");
  fTree->Branch("InDet",             &InDet,             "InDet/I");
  fTree->Branch("FidCut",            &FidCut,            "FidCut/I");
  fTree->Branch("TrueEnergy",        &trueEnergy,        "TrueEnergy/D");
  fTree->Branch("CNNEnergy",         &regcnn_energy,     "CNNEnergy/F");
  fTree->Branch("CNNVertex",         CNN_vertex,         "CNNVertex[3]/F");
  fTree->Branch("CNNVertex1ststep",  CNN_vertex1ststep,  "CNNVertex1ststep[3]/F");
  fTree->Branch("CNNVertex2ndstep",  CNN_vertex2ndstep,  "CNNVertex2ndstep[3]/F");
  fTree->Branch("CNNVertex2",        CNN_vertex2,        "CNNVertex2[3]/F");

  fTree->Branch("NuTruthN",          &nu_truth_N,    "NuTruthN/I");  
  fTree->Branch("NuEngTruth",        nueng_truth,    "NuEngTruth[NuTruthN]/D");   
  fTree->Branch("NuPDGTruth",        nupdg_truth,    "NuPDGTruth[NuTruthN]/I");   
  fTree->Branch("NuModeTruth",       numode_truth,   "NuModeTruth[NuTruthN]/I"); 
  fTree->Branch("NuCCNCTruth",       nuccnc_truth,   "NuCCNCTruth[NuTruthN]/I");   

  fTree->Branch("NuVtxXTruth",       nuvtxx_truth,   "NuVtxXTruth[NuTruthN]/D"); 
  fTree->Branch("NuVtxYTruth",       nuvtxy_truth,   "NuVtxYTruth[NuTruthN]/D"); 
  fTree->Branch("NuVtxZTruth",       nuvtxz_truth,   "NuVtxZTruth[NuTruthN]/D"); 

  fTree->Branch("PandNuVtx",         &pandora_nu_vtx,  "PandNuVtx/I"); 
  fTree->Branch("PandNuVertexX",     &pandora_nu_vtx_x,  "PandNuVertexX/D"); 
  fTree->Branch("PandNuVertexY",     &pandora_nu_vtx_y,  "PandNuVertexY/D"); 
  fTree->Branch("PandNuVertexZ",     &pandora_nu_vtx_z,  "PandNuVertexZ/D"); 

  fTree->Branch("NHits",             &nhits,         "NHits/I");

  fTree->Branch("ErecoNue",          &ErecoNue,      "ErecoNue/D");
  fTree->Branch("RecoLepEnNue",      &RecoLepEnNue,  "RecoLepEnNue/D");
  fTree->Branch("RecoHadEnNue",      &RecoHadEnNue,  "RecoHadEnNue/D");
  fTree->Branch("RecoMethodNue",     &RecoMethodNue, "RecoMethodNue/I");


}

void myana::RegCNNVertexAna::reconfigure(fhicl::ParameterSet const& pset)
{
  fMCGenModuleLabel  = pset.get<std::string>("MCGenModuleLabel");
  fHitsModuleLabel   = pset.get<std::string>("HitsModuleLabel");
  fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");

  fPandoraNuVertexModuleLabel = pset.get<std::string>("PandoraNuVertexModuleLabel");

  fRegCNNResultLabel = pset.get<std::string>("RegCNNEngResultLabel");
  fRegCNNModuleLabel = pset.get<std::string>("RegCNNEngModuleLabel");

  fRegCNNVtxResultLabel = pset.get<std::string>("RegCNNVtxResultLabel");
  fRegCNNVtx1ststepResultLabel = pset.get<std::string>("RegCNNVtx1ststepResultLabel");
  fRegCNNVtx2ndstepResultLabel = pset.get<std::string>("RegCNNVtx2ndstepResultLabel");
  fRegCNNVtxResult2Label = pset.get<std::string>("RegCNNVtxResult2Label");

  fRegCNNVtxModuleLabel = pset.get<std::string>("RegCNNVtxModuleLabel");
  fRegCNNVtx1ststepModuleLabel = pset.get<std::string>("RegCNNVtx1ststepModuleLabel");
  fRegCNNVtx2ndstepModuleLabel = pset.get<std::string>("RegCNNVtx2ndstepModuleLabel");
  fRegCNNVtxModule2Label = pset.get<std::string>("RegCNNVtxModule2Label");

}

void myana::RegCNNVertexAna::analyze(art::Event const& evt)
{
  this->reset();
  ievt = evt.id().event();
  bool isMC = !evt.isRealData(); 

  // * MC truth information 

  std::vector<art::Ptr<simb::MCTruth> > mclist;      
  if (isMC){    
    auto mctruthListHandle = evt.getHandle< std::vector<simb::MCTruth> >(fMCGenModuleLabel);
    if (mctruthListHandle)
      art::fill_ptr_vector(mclist, mctruthListHandle);  
  }  

  // Get the hits out of the event record
  std::vector<art::Ptr<recob::Hit> > hits;
  auto hitHandle = evt.getHandle<std::vector<recob::Hit> >(fHitsModuleLabel);
  if (hitHandle)
    art::fill_ptr_vector(hits, hitHandle);

  // Get DUNE energy Reco
  auto engrecoHandle = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNueLabel);

  // Get RegCNN Results
  art::InputTag itag(fRegCNNModuleLabel, fRegCNNResultLabel);
  auto cnnresultListHandle = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag);

  art::InputTag itag2(fRegCNNVtxModuleLabel, fRegCNNVtxResultLabel);
  auto cnnvtxresultListHandle = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag2);

  art::InputTag itag3(fRegCNNVtx1ststepModuleLabel, fRegCNNVtx1ststepResultLabel);
  auto cnnvtxresultListHandle1st = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag3);

  art::InputTag itag4(fRegCNNVtx2ndstepModuleLabel, fRegCNNVtx2ndstepResultLabel);
  auto cnnvtxresultListHandle2nd = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag4);

  art::InputTag itag5(fRegCNNVtxModule2Label, fRegCNNVtxResult2Label);
  auto cnnvtxresultListHandle2 = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag5);

  // Get Truth information
  if (mclist.size()>0)
  {
        int neutrino_i = 0;   
        for(size_t iList = 0; (iList < mclist.size()) && (neutrino_i < kMax) ; ++iList)
        {
          if (mclist[iList]->NeutrinoSet())
          {
            nueng_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().Momentum().E();
            nupdg_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().PdgCode(); 
            nuccnc_truth[neutrino_i] = mclist[iList]->GetNeutrino().CCNC(); 
            numode_truth[neutrino_i] = mclist[iList]->GetNeutrino().Mode();

            nuvtxx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vx();
            nuvtxy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vy();
            nuvtxz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vz();

            neutrino_i++; 
          }        
        }   
        nu_truth_N = neutrino_i;  
  }
  // Get Hit information
  nhits = hits.size();
  InDet = 0;
  //cut with true vertex in fiducial volume
  FidCut = 0;
  if (nu_truth_N>0){
    if(fabs(nuvtxx_truth[0]) < 310. && fabs(nuvtxy_truth[0]) < 550. && nuvtxz_truth[0] > 50. && nuvtxz_truth[0] < 1250.) FidCut = 1;
  }

  // Pandora Nu Vertex
  lar_pandora::PFParticleVector particleVector;
  lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
  lar_pandora::VertexVector vertexVector;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);

  pandora_nu_vtx = 0, pandora_nu_vtx_x = -10000, pandora_nu_vtx_y = -10000, pandora_nu_vtx_z = -10000;

  double xyz_temp[3] = {0.0, 0.0, 0.0} ;
  for (unsigned int ipfp = 0; ipfp < particleVector.size(); ipfp++){
    const art::Ptr<recob::PFParticle> particle = particleVector.at(ipfp);
    if (!particle->IsPrimary()) continue;

   // Particles <-> Vertices
   lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
   if (particlesToVertices.end() != vIter)
   {
       const lar_pandora::VertexVector &vertexVector = vIter->second;
       if (!vertexVector.empty())
       {
           if (vertexVector.size() !=1)
               std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

           const art::Ptr<recob::Vertex> vertex_pfp = *(vertexVector.begin());
           vertex_pfp->XYZ(xyz_temp);

	   pandora_nu_vtx = 1;
           pandora_nu_vtx_x = xyz_temp[0];
           pandora_nu_vtx_y = xyz_temp[1];
           pandora_nu_vtx_z = xyz_temp[2];
       }
   } // end of if
  } // end of loop particleVector



  // Get RecoE from DUNE
  if (!engrecoHandle.failedToGet())
  {
      ErecoNue          = engrecoHandle->fNuLorentzVector.E();
      RecoLepEnNue      = engrecoHandle->fLepLorentzVector.E();
      RecoHadEnNue      = engrecoHandle->fHadLorentzVector.E();
      RecoMethodNue     = engrecoHandle->recoMethodUsed;
      std::cout<< ErecoNue << std::endl;
  }

  // Get RegCNN Results
  if (!cnnresultListHandle.failedToGet())
  {
	if (!cnnresultListHandle->empty())
	{
  	  const std::vector<float>& v = (*cnnresultListHandle)[0].fOutput;
          for (unsigned int ii = 0; ii < 1; ii++){
	      regcnn_energy = v[ii];
          }
        }
  }


  // Get RegCNN Results
  if (!cnnvtxresultListHandle.failedToGet())
  {
	if (!cnnvtxresultListHandle->empty())
	{
  	  const std::vector<float>& v = (*cnnvtxresultListHandle)[0].fOutput;
          for (unsigned int ii = 0; ii < 3; ii++){
	      CNN_vertex[ii] = v[ii];
          }
        }
  }

  // Get RegCNN Results
  if (!cnnvtxresultListHandle2.failedToGet())
  {
	if (!cnnvtxresultListHandle2->empty())
	{
  	  const std::vector<float>& v = (*cnnvtxresultListHandle2)[0].fOutput;
          for (unsigned int ii = 0; ii < 3; ii++){
	      CNN_vertex2[ii] = v[ii];
          }
        }
  }

  // Get RegCNN Results
  if (!cnnvtxresultListHandle1st.failedToGet())
  {
	if (!cnnvtxresultListHandle1st->empty())
	{
  	  const std::vector<float>& v = (*cnnvtxresultListHandle1st)[0].fOutput;
          //std::cout << v[0] << std::endl;
          for (unsigned int ii = 0; ii < 3; ii++){
	      CNN_vertex1ststep[ii] = v[ii];
          }
        }
  }
  // Get RegCNN Results
  if (!cnnvtxresultListHandle2nd.failedToGet())
  {
	if (!cnnvtxresultListHandle2nd->empty())
	{
  	  const std::vector<float>& v = (*cnnvtxresultListHandle2nd)[0].fOutput;
          //std::cout << v[0] << std::endl;
          for (unsigned int ii = 0; ii < 3; ii++){
	      CNN_vertex2ndstep[ii] = v[ii];
          }
        }
  }



  // fill entry
  fTree->Fill();
}

void myana::RegCNNVertexAna::reset()
{
    ievt = -9999;
    trueEnergy = -99999;
    regcnn_energy = -99999;
    ErecoNue = -99999;
    RecoLepEnNue = -99999;
    RecoHadEnNue = -99999;
    RecoMethodNue = -99999;
   
    for (int ii = 0; ii < 3; ii++){ 
        CNN_vertex[ii] = -99999;
        CNN_vertex1ststep[ii] = -99999;
        CNN_vertex2ndstep[ii] = -99999;
        CNN_vertex2[ii] = -99999;
    }
    nu_truth_N = 0;
    for (int ii = 0; ii < kMax; ++ii)
    {
        nupdg_truth[ii] = -99999;
        numode_truth[ii] = -99999;
        nuccnc_truth[ii] = -99999;
        nueng_truth[ii] = -99999;

        nuvtxx_truth[ii] = -99999;
        nuvtxy_truth[ii] = -99999;
        nuvtxz_truth[ii] = -99999;

    }

}

DEFINE_ART_MODULE(myana::RegCNNVertexAna)
