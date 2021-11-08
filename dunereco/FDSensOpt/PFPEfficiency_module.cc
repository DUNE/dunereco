////////////////////////////////////////////////////////////////////////
// Class:       PFPEfficiency
// Plugin Type: analyzer (art v3_03_01)
// File:        PFPEfficiency_module.cc
//
// Generated at Tue Jan 14 20:15:19 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "TH1D.h"
#include "TEfficiency.h"

#include <map>
#include <string>

namespace dune {
  class PFPEfficiency;
}


class dune::PFPEfficiency : public art::EDAnalyzer {
public:
  explicit PFPEfficiency(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PFPEfficiency(PFPEfficiency const&) = delete;
  PFPEfficiency(PFPEfficiency&&) = delete;
  PFPEfficiency& operator=(PFPEfficiency const&) = delete;
  PFPEfficiency& operator=(PFPEfficiency&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

private:
  
  std::string fMCTruthModuleLabel;
  std::string fPFPModuleLabel;
  std::string fHitModuleLabel;

  float fFidVolCutX;
  float fFidVolCutY;
  float fFidVolCutZ;

  float fFidVolXmin;
  float fFidVolXmax;
  float fFidVolYmin;
  float fFidVolYmax;
  float fFidVolZmin;
  float fFidVolZmax;
  
  int    MC_isCC;
  int    MC_incoming_PDG;
  double MC_incoming_P[4];
  double MC_vertex[4];
  double MC_lepton_startMomentum[4];

  bool insideFV(double vertex[4]);

  //efficiency
  TH1D *h_den[6];
  TH1D *h_num[6];
  TEfficiency *eff[6];

  //efficiency for completeness and purity > 0.5
  TH1D *h_num_good[6];
  TEfficiency *eff_good[6];

  std::string names[6];
  int pdg[6];

};


dune::PFPEfficiency::PFPEfficiency(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCTruthModuleLabel(p.get<std::string>("MCTruthModuleLabel")),
  fPFPModuleLabel(p.get<std::string>("PFPModuleLabel")),
  fHitModuleLabel(p.get<std::string>("HitModuleLabel")),
  fFidVolCutX(p.get<float>("FidVolCutX")),
  fFidVolCutY(p.get<float>("FidVolCutY")),
  fFidVolCutZ(p.get<float>("FidVolCutZ"))
{
  names[0] = "piplus";
  names[1] = "piminus";
  names[2] = "proton";
  names[3] = "muon";
  names[4] = "el";
  names[5] = "gamma";
  
  pdg[0] = 211;
  pdg[1] = -211;
  pdg[2] = 2212;
  pdg[3] = 13;
  pdg[4] = 11;
  pdg[5] = 22;
}

void dune::PFPEfficiency::analyze(art::Event const& event)
{

  //!save neutrino's interaction info
  auto MCtruthHandle = event.getHandle<std::vector<simb::MCTruth> >(fMCTruthModuleLabel);
  std::vector<art::Ptr<simb::MCTruth>> MCtruthlist;
  art::fill_ptr_vector(MCtruthlist, MCtruthHandle);
  art::Ptr<simb::MCTruth> MCtruth; 

  MCtruth = MCtruthlist[0];
  if( MCtruth->NeutrinoSet() ){
    simb::MCNeutrino nu = MCtruth->GetNeutrino();
    if( nu.CCNC() == 0 ) MC_isCC = 1;
    else if ( nu.CCNC() == 1 ) MC_isCC = 0;
    simb::MCParticle neutrino = nu.Nu();
    MC_incoming_PDG = nu.Nu().PdgCode();
    const TLorentzVector& nu_momentum = nu.Nu().Momentum(0);
    nu_momentum.GetXYZT(MC_incoming_P);
    const TLorentzVector& vertex =neutrino.Position(0);
    vertex.GetXYZT(MC_vertex);
  }

  if (!MC_isCC) return;
  if (std::abs(MC_incoming_PDG)!=16) return;
  
  bool isFiducial =insideFV( MC_vertex );
  if( !isFiducial ) return;

  // * hits
  std::vector<art::Ptr<recob::Hit> > hitlist;
  auto hitListHandle = event.getHandle< std::vector<recob::Hit> >(fHitModuleLabel);
  if (hitListHandle)
    art::fill_ptr_vector(hitlist, hitListHandle);

  std::vector < art::Ptr < recob::PFParticle > > pfpList;
  auto pfpListHandle = event.getHandle < std::vector < recob::PFParticle > >(fPFPModuleLabel);
  if (pfpListHandle) {
    art::fill_ptr_vector(pfpList, pfpListHandle);
  }

  // Get all clusters
  std::vector < art::Ptr < recob::Cluster > > cluList;
  auto cluListHandle = event.getHandle < std::vector < recob::Cluster > >(fPFPModuleLabel);
  if (cluListHandle) {
    art::fill_ptr_vector(cluList, cluListHandle);
  }

  // Get cluster-PFParticle association
  art::FindManyP<recob::Cluster> fmcpfp(pfpListHandle, event, fPFPModuleLabel);

  // Get hit-cluster association
  art::FindManyP<recob::Hit> fmhc(cluListHandle, event, fPFPModuleLabel);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  std::map<int, int> hitmap;

  for (auto const& hit : hitlist){
    if (hit->WireID().Plane!=2) continue;
    std::map<int,double> trkide;
    std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
    for(size_t e = 0; e < TrackIDs.size(); ++e){
      trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
    } 
    double maxe = -1;
    double tote = 0;
    int TrackID = 0;
    for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
      tote += ii->second;
      if ((ii->second)>maxe){
        maxe = ii->second;
        TrackID = ii->first;
      }
    }
    const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
    if (particle){
      ++hitmap[TrackID];
    }
  }

  std::map<int, float> completeness;
  std::map<int, float> purity;

  for (size_t i = 0; i<pfpList.size(); ++i){
    std::map<int, int> hitmap0;
    if (fmcpfp.isValid()){
      // Get clusters associated with pfparticle
      auto const& clusters = fmcpfp.at(i);
      for (auto const & cluster : clusters){
        if (fmhc.isValid()){
          // Get hits associated with cluster
          auto const& hits = fmhc.at(cluster.key());
          //auto const& hits = hitsFromSlice.at(sliceid);
          //std::cout<<hits.size()<<std::endl;
          for (auto const& hit : hits){
            if (hit->WireID().Plane!=2) continue;
            std::map<int,double> trkide;
            std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
            for(size_t e = 0; e < TrackIDs.size(); ++e){
              trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
            } 
            double maxe = -1;
            double tote = 0;
            int TrackID = 0;
            for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
              tote += ii->second;
              if ((ii->second)>maxe){
                maxe = ii->second;
                TrackID = ii->first;
              }
            }
            const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
            if (particle){
              ++hitmap0[TrackID];
            }
          }
        }
      }
    }
    int maxhits = -1;
    int tothits = 0;
    int TrackID = 0;
    for (std::map<int,int>::iterator ii = hitmap0.begin(); ii!=hitmap0.end(); ++ii){
      tothits += ii->second;
      if ((ii->second)>maxhits){
        maxhits = ii->second;
        TrackID = ii->first;
      }
    }
    if (TrackID){
      double new_completeness = 1.0*maxhits/hitmap[TrackID];
      double new_purity = 1.0*maxhits/tothits;
      //const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
      //if (particle->PdgCode() == 13) std::cout<<new_completeness<<" "<<new_purity<<std::endl;
      if (new_completeness*new_purity > completeness[TrackID]*purity[TrackID]){ 
        completeness[TrackID] = new_completeness;
        purity[TrackID] = new_purity;
      }
    }
  }

  const sim::ParticleList& plist = pi_serv->ParticleList();
  for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    simb::MCParticle *particle = ipar->second;
    if (particle->Mother() == 0){ //primary particle
      for (int i = 0; i<6; ++i){
        if (particle->PdgCode() == pdg[i]){
          h_den[i]->Fill(particle->P());
          if (completeness[particle->TrackId()]>0){
            h_num[i]->Fill(particle->P());
            //if (i==3) std::cout<<completeness[particle->TrackId()]<<" "<<purity[particle->TrackId()]<<std::endl;
            if (completeness[particle->TrackId()]>0.5&&
                purity[particle->TrackId()]>0.5){
              h_num_good[i]->Fill(particle->P());
            }
          }
        }
      }
    }
  }
}

void dune::PFPEfficiency::beginJob(){
  // Get geometry.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  std::cout<<"Fiducial volume:"<<"\n"
	   <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
	   <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
	   <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";

  art::ServiceHandle<art::TFileService const> tfs;
  
  for (int i = 0; i<6; ++i){
    h_den[i] = tfs->make<TH1D>(Form("h_den_%s",names[i].c_str()), Form("%s;Momentum (GeV/c)", names[i].c_str()), 40, 0, 20);
    h_num[i] = tfs->make<TH1D>(Form("h_num_%s",names[i].c_str()), Form("%s;Momentum (GeV/c)", names[i].c_str()), 40, 0, 20);
    h_num_good[i] = tfs->make<TH1D>(Form("h_num_good_%s",names[i].c_str()), Form("%s;Momentum (GeV/c)", names[i].c_str()), 40, 0, 20);
    h_den[i]->Sumw2();
    h_num[i]->Sumw2();
    h_num_good[i]->Sumw2();
  }
}

void dune::PFPEfficiency::endJob(){

  art::ServiceHandle<art::TFileService const> tfs;

  for (int i = 0; i<6; ++i){
    if (TEfficiency::CheckConsistency(*h_num[i], *h_den[i])){
      eff[i] = tfs->make<TEfficiency>(*h_num[i], *h_den[i]);
      eff[i]->SetTitle(names[i].c_str());
      eff[i]->Write(Form("eff_%s",names[i].c_str()));
    }
    if (TEfficiency::CheckConsistency(*h_num_good[i], *h_den[i])){
      eff_good[i] = tfs->make<TEfficiency>(*h_num_good[i], *h_den[i]);
      eff_good[i]->SetTitle((names[i]+", completeness>0.5, purity>0.5").c_str());
      eff_good[i]->Write(Form("eff_good_%s",names[i].c_str()));
    }
  }
}
bool dune::PFPEfficiency::insideFV( double vertex[4]){

     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     if (x>fFidVolXmin && x<fFidVolXmax&&
	 y>fFidVolYmin && y<fFidVolYmax&&
	 z>fFidVolZmin && z<fFidVolZmax)
       return true;
     else
       return false;
}

DEFINE_ART_MODULE(dune::PFPEfficiency)
