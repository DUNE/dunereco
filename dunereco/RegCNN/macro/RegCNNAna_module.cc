////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       RegCNNAna
// Module Type: analyzer
// File:        RegCNNAna_module.cc
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
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/RegCNN/func/RegCNNResult.h"

const int kMax = 1000;
namespace myana {
  class RegCNNAna;
}

class myana::RegCNNAna : public art::EDAnalyzer
{
  public:
    RegCNNAna(fhicl::ParameterSet const& pset);
    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& p);
    void reset();
  
  private:
    TTree* fTree; 
    int ievt;
    double trueEnergy;
    float  regcnn_energy;
    float  regcnn_vertex[3];
    
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

    std::string fMCGenModuleLabel;   
    std::string fRegCNNModuleLabel;   
    std::string fHitsModuleLabel;
    std::string fEnergyRecoNueLabel;
    std::string fRegCNNResultLabel;   

    art::ServiceHandle<art::TFileService> tfs;

};

myana::RegCNNAna::RegCNNAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
  fTree = tfs->make<TTree>("anatree", "anatree");
  fTree->Branch("ievent",            &ievt,          "ievent/I");
  fTree->Branch("InDet",             &InDet,         "InDet/I");
  fTree->Branch("FidCut",            &FidCut,         "FidCut/I");
  fTree->Branch("TrueEnergy",        &trueEnergy,    "TrueEnergy/D");
  fTree->Branch("CNNEnergy",         &regcnn_energy, "CNNEnergy/F");
  fTree->Branch("CNNVertex",         regcnn_vertex, "CNNVertex[3]/F");

  fTree->Branch("NuTruthN",          &nu_truth_N,    "NuTruthN/I");  
  fTree->Branch("NuEngTruth",        nueng_truth,    "NuEngTruth[NuTruthN]/D");   
  fTree->Branch("NuPDGTruth",        nupdg_truth,    "NuPDGTruth[NuTruthN]/I");   
  fTree->Branch("NuModeTruth",       numode_truth,   "NuModeTruth[NuTruthN]/I"); 
  fTree->Branch("NuCCNCTruth",       nuccnc_truth,   "NuCCNCTruth[NuTruthN]/I");   

  fTree->Branch("NuVtxXTruth",       nuvtxx_truth,   "NuVtxXTruth[NuTruthN]/I"); 
  fTree->Branch("NuVtxYTruth",       nuvtxy_truth,   "NuVtxYTruth[NuTruthN]/I"); 
  fTree->Branch("NuVtxZTruth",       nuvtxz_truth,   "NuVtxZTruth[NuTruthN]/I"); 

  fTree->Branch("NHits",             &nhits,         "NHits/I");

  fTree->Branch("ErecoNue",          &ErecoNue,      "ErecoNue/D");
  fTree->Branch("RecoLepEnNue",      &RecoLepEnNue,  "RecoLepEnNue/D");
  fTree->Branch("RecoHadEnNue",      &RecoHadEnNue,  "RecoHadEnNue/D");
  fTree->Branch("RecoMethodNue",     &RecoMethodNue, "RecoMethodNue/I");


}

void myana::RegCNNAna::reconfigure(fhicl::ParameterSet const& pset)
{
  fMCGenModuleLabel  = pset.get<std::string>("MCGenModuleLabel");
  fRegCNNModuleLabel = pset.get<std::string>("RegCNNModuleLabel");
  fHitsModuleLabel   = pset.get<std::string>("HitsModuleLabel");
  fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
  fRegCNNResultLabel = pset.get<std::string>("RegCNNResultLabel");
}

void myana::RegCNNAna::analyze(art::Event const& evt)
{
  this->reset();
  ievt = evt.id().event();
  bool isMC = !evt.isRealData(); 

  // * MC truth information 
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;             
  std::vector<art::Ptr<simb::MCTruth> > mclist;      
  if (isMC){    
    if (evt.getByLabel(fMCGenModuleLabel,mctruthListHandle))         
      art::fill_ptr_vector(mclist, mctruthListHandle);  
  }  

  // Get the hits out of the event record
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Get DUNE energy Reco
  art::Handle<dune::EnergyRecoOutput> engrecoHandle;
  evt.getByLabel(fEnergyRecoNueLabel,engrecoHandle);


  // Get RegCNN Results
  art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle;
  evt.getByLabel(fRegCNNModuleLabel, fRegCNNResultLabel, cnnresultListHandle);
  //std::vector<art::Ptr<cnn::Result> > cnnlist;
  //if (evt.getByLabel(fRegCNNModuleLabel, cnnresultListHandle))
  //  art::fill_ptr_vector(cnnlist, cnnresultListHandle);


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
  if (hits.size()>0)
  {
      bool fid_flag = false;
      for(size_t iHit = 0; iHit < hits.size(); ++iHit)
      {
        art::Ptr<recob::Hit> hit = hits.at(iHit);
    	float peakT = hit->PeakTime();
    	//unsigned int channel = hit->Channel();
	unsigned int plane = hit->WireID().Plane;
        unsigned int wire = hit->WireID().Wire;
        unsigned int tpc  = hit->WireID().TPC;
        // for dunefd 1x2x6
    	if (peakT > 4482.) fid_flag = true; 
        if (plane == 2){
            if (tpc < 4 && wire < 5) fid_flag = true;
            if (tpc > 19 && wire > 474) fid_flag = true;
        }
	if (fid_flag) break; 
      }
      if (!fid_flag) InDet = 1; 
  }

  //cut with true vertex in fiducial volume
  FidCut = 0;
  if (nu_truth_N>0){
    //if(nuccnc_truth[0] == 0 && fabs(nuvtxx_truth[0]) < 310. && fabs(nuvtxy_truth[0]) < 550. && fabs(nuvtxz_truth[0]) > 50. && fabs(nuvtxz_truth[0]) < 1250.) nutrue_fid = 1; 
    if(fabs(nuvtxx_truth[0]) < 310. && fabs(nuvtxy_truth[0]) < 550. && nuvtxz_truth[0] > 50. && nuvtxz_truth[0] < 1250.) FidCut = 1;
  }


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
	  regcnn_energy = v[0];
          //std::cout << v[0] << std::endl;
          for (unsigned int ii = 0; ii < 3; ii++){
	      regcnn_vertex[ii] = v[ii];
          }
        }
  }
  // fill entry
  fTree->Fill();
}

void myana::RegCNNAna::reset()
{
    ievt = -9999;
    trueEnergy = -99999;
    regcnn_energy = -99999;
    ErecoNue = -99999;
    RecoLepEnNue = -99999;
    RecoHadEnNue = -99999;
    RecoMethodNue = -99999;
   
    for (int ii = 0; ii < 3; ii++){ 
        regcnn_vertex[ii] = -99999;
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

DEFINE_ART_MODULE(myana::RegCNNAna)
