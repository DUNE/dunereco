////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       RegCNNAna
// Module Type: analyzer
// File:        RegCNNAna_module.cc for analyzing CNN and Standard Reco. results
// Author:      Ilsoo Seong - iseong@uci.edu
//              Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////////////////////////////

//STL
#include <limits>
#include <string> 
#include <vector>
#include <algorithm>
//ROOT
#include "TTree.h"
#include "TBranch.h"
//ART
#include "art/Framework/Core/ModuleMacros.h"  
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/EDAnalyzer.h"
//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
//DUNE
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/RegCNN/func/RegCNNResult.h"
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"

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

        art::Ptr<recob::Track> GetLongestTrack(const art::Event& event);
        art::Ptr<recob::Shower> GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                       detinfo::DetectorPropertiesData const& detProp,
                                                       const art::Event& event);

        /**
         * @brief  Converts deposited charge into energy by converting to number of electrons and correcting for average recombination
         *
         * @param  charge the deposited charge
         *
         * @return the reconstructed deposited energy
         */
        double CalculateEnergyFromCharge(const double charge);
        std::vector<std::pair<const simb::MCParticle*, double> >  get_sortedMCParticle(
                std::unordered_map<const simb::MCParticle*, double> mcEMap);

        TTree* fTree; 
        int ievt;
        float  regcnn_energy;
        float  regcnn_vertex[3];
        float  regcnn_dir[3];

        int    InDet;
        int    FidCut;
        int    mutrk_contain;
        int    nhits;
        int    nu_truth_N;
        int    nupdg_truth[kMax];
        int    numode_truth[kMax];    
        int    nuccnc_truth[kMax];     
        double nueng_truth[kMax];
        double nuvtxx_truth[kMax];
        double nuvtxy_truth[kMax];
        double nuvtxz_truth[kMax];
        int    lepPDG_truth[kMax];
        double lepEng_truth[kMax];

        double ErecoNu;
        double RecoLepEnNu; 
        double RecoHadEnNu; 
        int    RecoMethodNu; 
        int    RecoTrackMomMethod;

	    int n_track_pad;
        int track_id[kMax];
	    double all_track_px[kMax];
	    double all_track_py[kMax];
	    double all_track_pz[kMax];
	    double all_track_true_px[kMax];
	    double all_track_true_py[kMax];
	    double all_track_true_pz[kMax];

        //uncorrected hadronic energy
        double fUncorrectedHadEn;
        double RecoMuTrackLength;
        double fUncorrectedMuMomMCS;

        calo::CalorimetryAlg fCalorimetryAlg;         ///< the calorimetry algorithm
        std::string fMCGenModuleLabel;   
        std::string fRegCNNModuleLabel;   
        std::string fRegCNNResultLabel;   
        std::string fRegCNNDirModuleLabel;   
        std::string fRegCNNDirResultLabel;
        std::string fHitsModuleLabel;
        std::string fEnergyRecoNuLabel;
        std::string fTrackToHitLabel;
        std::string fShowerToHitLabel;
        std::string fParticleLabel;
        std::string fTrackLabel;
        std::string fTrackLabelDir;
        std::string fShowerLabel;
        double fRecombFactor;                          ///< the average reccombination factor

        art::ServiceHandle<art::TFileService> tfs;
        art::ServiceHandle<cheat::BackTrackerService> backtracker;
        art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;

};

myana::RegCNNAna::RegCNNAna(fhicl::ParameterSet const& pset) : 
    EDAnalyzer      (pset),
    fCalorimetryAlg (pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
    this->reconfigure(pset);

    fTree = tfs->make<TTree>("anatree", "anatree");
    fTree->Branch("ievent",              &ievt,           "ievent/I");
    fTree->Branch("InDet",               &InDet,          "InDet/I");
    fTree->Branch("FidCut",              &FidCut,         "FidCut/I");
    fTree->Branch("MuTrackCont",         &mutrk_contain,  "MuTrackCont/I");
    fTree->Branch("CNNEnergy",           &regcnn_energy,  "CNNEnergy/F");
    fTree->Branch("CNNVertex",           regcnn_vertex,   "CNNVertex[3]/F");
    fTree->Branch("CNNDir",              regcnn_dir,      "CNNDir[3]/F");

    // MC truth
    fTree->Branch("NuTruthN",            &nu_truth_N,     "NuTruthN/I");  
    fTree->Branch("NuEngTruth",          nueng_truth,     "NuEngTruth[NuTruthN]/D");   
    fTree->Branch("NuPDGTruth",          nupdg_truth,     "NuPDGTruth[NuTruthN]/I");   
    fTree->Branch("NuModeTruth",         numode_truth,    "NuModeTruth[NuTruthN]/I"); 
    fTree->Branch("NuCCNCTruth",         nuccnc_truth,    "NuCCNCTruth[NuTruthN]/I");   
    fTree->Branch("NuVtxXTruth",         nuvtxx_truth,    "NuVtxXTruth[NuTruthN]/D"); 
    fTree->Branch("NuVtxYTruth",         nuvtxy_truth,    "NuVtxYTruth[NuTruthN]/D"); 
    fTree->Branch("NuVtxZTruth",         nuvtxz_truth,    "NuVtxZTruth[NuTruthN]/D"); 

    // lepton truth
    fTree->Branch("LepPDGTruth",         lepPDG_truth,    "LepPDGTruth[NuTruthN]/I");
    fTree->Branch("LepEngTruth",         lepEng_truth,    "LepEngTruth[NuTruthN]/D");

    // Reco. from DUNE, energy
    fTree->Branch("NHits",               &nhits,          "NHits/I");
    fTree->Branch("ErecoNu",             &ErecoNu,        "ErecoNu/D");
    fTree->Branch("RecoLepEnNu",         &RecoLepEnNu,    "RecoLepEnNu/D");
    fTree->Branch("RecoHadEnNu",         &RecoHadEnNu,    "RecoHadEnNu/D");
    fTree->Branch("RecoMethodNu",        &RecoMethodNu,   "RecoMethodNu/I");
    fTree->Branch("RecoTrackMomMethod",  &RecoTrackMomMethod, "RecoTrackMomMethod/I");
    fTree->Branch("RecoMuTrackLength",   &RecoMuTrackLength,  "RecoMuTrackLength/D");
    fTree->Branch("UncorrectedHadEn",    &fUncorrectedHadEn,    "UncorrectedHadEn/D");
    fTree->Branch("UncorrectedMuMomMCS", &fUncorrectedMuMomMCS, "UncorrectedMuMomMCS/D");
    
    // direction
    fTree->Branch("n_track_pad", &n_track_pad, "n_track_pad/I");
    fTree->Branch("track_id", track_id, "track_id[n_track_pad]/I");
    fTree->Branch("all_track_px", all_track_px, "all_track_px[n_track_pad]/D");
    fTree->Branch("all_track_py", all_track_py, "all_track_py[n_track_pad]/D");
    fTree->Branch("all_track_pz", all_track_pz, "all_track_pz[n_track_pad]/D");
    fTree->Branch("all_track_true_px", all_track_true_px, "all_track_true_px[n_track_pad]/D");
    fTree->Branch("all_track_true_py", all_track_true_py, "all_track_true_py[n_track_pad]/D");
    fTree->Branch("all_track_true_pz", all_track_true_pz, "all_track_true_pz[n_track_pad]/D");
}

void myana::RegCNNAna::reconfigure(fhicl::ParameterSet const& pset)
{
    //fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"));
    fMCGenModuleLabel     = pset.get<std::string>("MCGenModuleLabel");
    fRegCNNModuleLabel    = pset.get<std::string>("RegCNNModuleLabel");
    fRegCNNResultLabel    = pset.get<std::string>("RegCNNResultLabel");
    fRegCNNDirModuleLabel = pset.get<std::string>("RegCNNDirModuleLabel");
    fRegCNNDirResultLabel = pset.get<std::string>("RegCNNDirResultLabel");
    fHitsModuleLabel      = pset.get<std::string>("HitsModuleLabel");
    fEnergyRecoNuLabel    = pset.get<std::string>("EnergyRecoNuLabel");
    fTrackToHitLabel      = pset.get<std::string>("TrackToHitLabel");
    fShowerToHitLabel     = pset.get<std::string>("ShowerToHitLabel");
    fParticleLabel        = pset.get<std::string>("ParticleLabel");
    fTrackLabel           = pset.get<std::string>("TrackLabel");
    fTrackLabelDir        = pset.get<std::string>("TrackLabelDir");
    fShowerLabel          = pset.get<std::string>("ShowerLabel");
    fRecombFactor         = pset.get<double>("RecombFactor");
}

void myana::RegCNNAna::analyze(art::Event const& evt)
{
    this->reset();
    ievt = evt.id().event();
    bool isMC = !evt.isRealData(); 

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    art::Ptr<recob::Shower> highestChargeShower(this->GetHighestChargeShower(clockData, detProp, evt));
    // Get the longest track of current event and then the associated hits of that track
    // Get lepton charge from track hits and event charge from event hits
    // Get hadronic charge by the subtraction of lepton charge from event charge
    // Get hadronic energy from hadronic charge
    // Compare hadronic energy to true energy, get calibrated parameters: grad and interception
    art::Ptr<recob::Track> longestTrack(this->GetLongestTrack(evt));
    if (!longestTrack.isAvailable() || longestTrack.isNull()) {
        mf::LogWarning("RegCNNAna")<<" Cannot access the muon track which is needed for this energy reconstruction method.\n"
            <<"Don't set value"<<std::endl;
    } else {
        const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(longestTrack, evt, fTrackToHitLabel), 2));
        const double leptonObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, muonHits));
        const std::vector<art::Ptr<recob::Hit> > eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel), 2));
        const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits));
        const double hadronicObservedCharge(eventObservedCharge-leptonObservedCharge);
        fUncorrectedHadEn = this->CalculateEnergyFromCharge(hadronicObservedCharge);
        RecoMuTrackLength = longestTrack->Length();
        trkf::TrackMomentumCalculator TrackMomCalc;
        fUncorrectedMuMomMCS = TrackMomCalc.GetMomentumMultiScatterChi2(longestTrack);
        std::cout<<"track length :"<<RecoMuTrackLength<<", uncorrectedMuMomMCS :"<<fUncorrectedMuMomMCS<<std::endl;
    }

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

    // Get the pandoratrack out of the event record
    art::Handle<std::vector<recob::Track> > trackHandleDir;
    std::vector<art::Ptr<recob::Track> > trackDir;
    if (evt.getByLabel(fTrackLabelDir, trackHandleDir))
        art::fill_ptr_vector(trackDir, trackHandleDir);
    art::FindManyP<recob::Hit> fmtrkDir(trackHandleDir, evt, fTrackLabelDir);

    if (trackHandleDir.isValid()) {
        n_track_pad = trackDir.size();
        for (int itrack= 0; itrack< n_track_pad; ++itrack) {
            art::Ptr<recob::Track> ptrack = trackDir.at(itrack);
            track_id[itrack] = ptrack->ID();
            all_track_px[itrack] = ptrack->VertexDirection().X();
            all_track_py[itrack] = ptrack->VertexDirection().Y();
            all_track_pz[itrack] = ptrack->VertexDirection().Z();

            if (fmtrkDir.isValid()) {
                std::vector<art::Ptr<recob::Hit> > vhit = fmtrkDir.at(itrack);
                std::unordered_map<const simb::MCParticle*, double> mcEMap;
                for (size_t ihit= 0; ihit< vhit.size(); ++ihit) {
                    std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackIDEs(clockData, vhit[ihit]);
                    for (size_t e= 0; e< trackIDs.size(); ++e) {
                        mcEMap[particleinventory->TrackIdToParticle_P(trackIDs[e].trackID)] += trackIDs[e].energy;
                    }
                }
                std::vector<std::pair<const simb::MCParticle*, double> > trkTrue = get_sortedMCParticle(mcEMap);
                if (trkTrue.size()> 0) {
                    TVector3 v3_true(trkTrue[0].first->Momentum().Vect());
                    all_track_true_px[itrack] = v3_true.X();
                    all_track_true_py[itrack] = v3_true.Y();
                    all_track_true_pz[itrack] = v3_true.Z();
                }
            }
        }
    }

    // Get DUNE energy Reco
    art::Handle<dune::EnergyRecoOutput> engrecoHandle;
    evt.getByLabel(fEnergyRecoNuLabel,engrecoHandle);

    // Get RegCNN Results
    art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle;
    evt.getByLabel(fRegCNNModuleLabel, fRegCNNResultLabel, cnnresultListHandle);

    art::Handle<std::vector<cnn::RegCNNResult> > RegCnnDirResultListHandle;
    evt.getByLabel(fRegCNNDirModuleLabel, fRegCNNDirResultLabel, RegCnnDirResultListHandle);

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

                lepEng_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().E();
                lepPDG_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().PdgCode();

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
        if(fabs(nuvtxx_truth[0]) < 310. && fabs(nuvtxy_truth[0]) < 550. && nuvtxz_truth[0] > 50. && nuvtxz_truth[0] < 1250.) FidCut = 1;
        if (fabs(nupdg_truth[0])==14) {
            mutrk_contain = engrecoHandle->longestTrackContained; // 1 = contained, 0 = exiting, -1 = not set
        }
    }

    // Get RecoE from DUNE
    if (!engrecoHandle.failedToGet())
    {
        ErecoNu           = engrecoHandle->fNuLorentzVector.E();
        RecoLepEnNu       = engrecoHandle->fLepLorentzVector.E();
        RecoHadEnNu       = engrecoHandle->fHadLorentzVector.E();
        RecoMethodNu      = engrecoHandle->recoMethodUsed;  // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
        RecoTrackMomMethod = engrecoHandle->trackMomMethod; // 1 = range, 0 = MCS, -1 = not set
        std::cout<< "EnergyReco: "<<ErecoNu << std::endl;
    }

    // Get RegCNN Results
    if (!cnnresultListHandle.failedToGet())
    {
        if (!cnnresultListHandle->empty())
        {
            const std::vector<float>& v = (*cnnresultListHandle)[0].fOutput;
            regcnn_energy = v[0];
            for (unsigned int ii = 0; ii < 3; ii++){
                regcnn_vertex[ii] = v[ii];
            }
        }
    }

    if (!RegCnnDirResultListHandle.failedToGet()) {
        if (!RegCnnDirResultListHandle->empty()) {
            const std::vector<float>& v = (*RegCnnDirResultListHandle)[0].fOutput;
            for (unsigned int i= 0; i< v.size(); ++i) {
                regcnn_dir[i] = v[i];
            }
        }
    }

    // fill entry
    fTree->Fill();
}

void myana::RegCNNAna::reset()
{
    // set default to nonsense value
    ievt               = -9999;
    mutrk_contain      = -999;
    regcnn_energy      = -99999;
    ErecoNu            = -99999;
    RecoLepEnNu        = -99999;
    RecoHadEnNu        = -99999;
    RecoMethodNu       = -99999;
    RecoTrackMomMethod = -99999;
    
    fUncorrectedHadEn    = -99999;
    RecoMuTrackLength    = -99999;
    fUncorrectedMuMomMCS = -99999;

    for (int ii = 0; ii < 3; ii++){ 
        regcnn_vertex[ii] = -99999;
        regcnn_dir[ii] = -99999;
    }

    nu_truth_N = 0;
    n_track_pad = 0;
    for (int ii = 0; ii < kMax; ++ii)
    {
        nupdg_truth[ii]  = -99999;
        numode_truth[ii] = -99999;
        nuccnc_truth[ii] = -99999;
        nueng_truth[ii]  = -99999;
        lepPDG_truth[ii] = -99999;
        lepEng_truth[ii] = -99999;

        nuvtxx_truth[ii] = -99999;
        nuvtxy_truth[ii] = -99999;
        nuvtxz_truth[ii] = -99999;

        track_id[ii] = -99999;
	    all_track_px[ii] = -99999;
	    all_track_py[ii] = -99999;
	    all_track_pz[ii] = -99999;
	    all_track_true_px[ii] = -99999;
	    all_track_true_py[ii] = -99999;
	    all_track_true_pz[ii] = -99999;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Track> myana::RegCNNAna::GetLongestTrack(const art::Event &event)
{
    art::Ptr<recob::Track> pTrack(art::Ptr<recob::Track>(art::ProductID("nullTrack")));
    const std::vector<art::Ptr<recob::Track> > tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
    if (0 == tracks.size())
        return pTrack;

    double longestLength(std::numeric_limits<double>::lowest());
    for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
    {
        const double length(tracks[iTrack]->Length());
        if (length-longestLength > std::numeric_limits<double>::epsilon())
        {
            longestLength = length;
            pTrack = tracks[iTrack];
        }
    }

    /***********************
    // get all PFParticles first, and then get tracks by looping all the PFParticles
    art::Ptr<recob::Track> pTrack(art::Ptr<recob::Track>(art::ProductID("nullTrack")));
    // Get all PFParticles
    const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(event, fParticleLabel);
    if (0 == particles.size())
        return pTrack;

    double longestLength(std::numeric_limits<double>::lowest());
    for (const art::Ptr<recob::PFParticle> &particle : particles) {
        // Get the track if this particle is track-like
        if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(particle, event, fParticleLabel, fTrackLabel)) {
            const art::Ptr<recob::Track> trk = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle, event, fParticleLabel, fTrackLabel);
            const double length(trk->Length());
            if (length-longestLength > std::numeric_limits<double>::epsilon()) {
                longestLength = length;
                pTrack = trk;
            }
        }
    } // end of PFParticles
    ************************/

    return pTrack;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Shower> myana::RegCNNAna::GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                           detinfo::DetectorPropertiesData const& detProp,
                                                           const art::Event &event)
{
    art::Ptr<recob::Shower> pShower(art::Ptr<recob::Shower>(art::ProductID("nullShower")));
    const std::vector<art::Ptr<recob::Shower> > showers(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));
    if (0 == showers.size())
        return pShower;

    double maxCharge(std::numeric_limits<double>::lowest());
    for (unsigned int iShower = 0; iShower < showers.size(); ++iShower)
    {
        const std::vector<art::Ptr<recob::Hit> > showerHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(showers[iShower],
            event,fShowerToHitLabel),2));
        const double showerCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, showerHits));
        if (showerCharge-maxCharge > std::numeric_limits<double>::epsilon())
        {
            maxCharge = showerCharge;
            pShower = showers[iShower];
        }
    }

    /*************************
    // get all PFParticles first, and then get tracks by looping all the PFParticles
    art::Ptr<recob::Shower> pShower(art::Ptr<recob::Shower>(art::ProductID("nullShower")));
    // Get all PFParticles
    const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(event, fParticleLabel);
    if (0 == particles.size())
        return pShower;

    double maxCharge(std::numeric_limits<double>::lowest());
    for (const art::Ptr<recob::PFParticle> &particle : particles) {
        // Get the shower if this particle is shower-like
        if (dune_ana::DUNEAnaPFParticleUtils::IsShower(particle, event, fParticleLabel, fShowerLabel)) {
            const art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(particle, event, fParticleLabel, fShowerLabel);
            const std::vector<art::Ptr<recob::Hit> > showerHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(shower, event, fShowerToHitLabel), 2));
            const double showerCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(showerHits));
            if (showerCharge-maxCharge > std::numeric_limits<double>::epsilon()) {
                maxCharge = showerCharge;
                pShower = shower;
            }
        }
    } // end of PFParticles
    ************************/

    return pShower;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

double myana::RegCNNAna::CalculateEnergyFromCharge(const double charge)
{
    return fCalorimetryAlg.ElectronsFromADCArea(charge,2)*1./fRecombFactor/util::kGeVToElectrons;
}

std::vector<std::pair<const simb::MCParticle*, double> >  
myana::RegCNNAna::get_sortedMCParticle(std::unordered_map<const simb::MCParticle*, double> mcEMap) {
    std::vector<std::pair<const simb::MCParticle*, double> > outVec;
    double total_E = 0;
    for (std::pair<const simb::MCParticle*, double> const& p: mcEMap) {
        outVec.push_back(p);
        total_E += p.second;
    }
    std::sort(outVec.begin(), outVec.end(), [](std::pair<const simb::MCParticle*, double> a, std::pair<const simb::MCParticle*, double> b){return a.second > b.second;});
    if (abs(total_E) < 1e-6) {total_E = 1;} // protect agains zero division
    for (std::pair<const simb::MCParticle*, double> & p: outVec) {
        p.second /= total_E;
    }

    return outVec;
}

DEFINE_ART_MODULE(myana::RegCNNAna)
