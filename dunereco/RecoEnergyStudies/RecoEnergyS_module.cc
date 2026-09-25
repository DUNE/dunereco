////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       RecoEnergyS
// Module Type: analyzer
// File:        RecoEnergyS_module.cc
// Author:      Ilsoo Seong
//
////////////////////////////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/Exceptions.h"// geo::InvalidWireError
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/Filters/ChannelFilter.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/sim.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/CVN/func/InteractionType.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/CVN/art/PixelMapProducer.h"
//#include "dune/RegCNN/func/RegCNNResult.h"

// ROOT & C++ includes
#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector3.h"

const int kMaxHits = 50000;
const int kMax     = 3000;
//const int nADCvec  = 4492;
const int kMaxClst = 100;
const int daughterHitCut = 100;

namespace lar_pandora{class LArPandoraHelper;}

namespace recoestudy {
  class RecoEnergyS;
}

class recoestudy::RecoEnergyS : public art::EDAnalyzer {
public:

  RecoEnergyS(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& evt);
  void reconfigure(fhicl::ParameterSet const& p);
  void reset();
  int FindTrackID(art::Ptr<recob::Hit> const& hit, detinfo::DetectorClocksData const& clockData);
  int FindShowerVtx(void);
  double GetPixelMapWireIdx(const geo::WireID& wireID);
  int ShiftGlobalWireIdx(double *global_idx, double *global_tick, int *tpc_idx, int *plane_idx, double *offsets);
  void GetDUNEGlobalWireTDC(const geo::WireID& wireID, double localTDC,
                           unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC, double driftV) const;
  bool ConvertXYZtoWireTDC(const double *vtx_loc, double *TPC, double *LocalWire, double *LocalTick,
			   double *GlobalPlane, double *GlobalWire, double *GlobalTick, detinfo::DetectorPropertiesData const& detProp);

  const simb::MCParticle* get_bestMCParticle(const art::Ptr<recob::Hit>& hit, detinfo::DetectorClocksData const& clockData);
  bool passesNHitsCut(const simb::MCParticle *part);
  void addDaughters(const simb::MCParticle *part, const std::vector<int> finalStatePdgs, std::vector< const simb::MCParticle * > &allParts);


private:
  calo::CalorimetryAlg fCaloAlg;

  // Variables for the tree
  TTree* fTree;
  int ievt;
  int run;
  int subrun;
  int    ntpcs_evt;
  double det_t0;
  double driftvelocity;
  double SampleRate;
  int nwires_p0;
  int nwires_p1;
  int nwires_p2;

  int primpdg;
  double trueEnergy;
  double trueVertex_X;
  double trueVertex_Y;
  double trueVertex_Z;
  double trueEnd_X;
  double trueEnd_Y;
  double trueEnd_Z;
  double truePx;
  double truePy;
  double truePz;
  int    vtx_in_tpc;
  double trueVtxWire[3];
  double trueVtxGlobalWire[3];
  double trueVtxGlobalTick[3];
  double trueVtxGlobalPlane[3];
  double trueVtxTPC[3];
  double trueVtxTime[3];
  double trueVtxChan[3];
  double trueEndWire[3];
  double trueEndTime[3];
  double trueEndGlobalWire[3];
  double trueEndGlobalTime[3];
  double trueMidWire[3];
  double trueMidTime[3];

  double CVNVtxWire[3];
  double CVNVtxTPC[3];
  double CVNVtxTime[3];
  double CVNVtxGlobalWire[3];
  double CVNVtxGlobalTick[3];
  double CVNVtxGlobalPlane[3];

  double CVNVtx2Wire[3];
  double CVNVtx2TPC[3];
  double CVNVtx2Time[3];
  double CVNVtx2GlobalWire[3];
  double CVNVtx2GlobalTick[3];
  double CVNVtx2GlobalPlane[3];


  double EndProngGlobalWire[3];
  double EndProngGlobalTick[3];
  double EndProngGlobalPlane[3];

  int    nutrue_fid;
  int    nu_truth_N;
  double nueng_truth[kMax];
  int    nupdg_truth[kMax];
  double nuvtxx_truth[kMax];
  double nuvtxy_truth[kMax];
  double nuvtxz_truth[kMax];
  int    numode_truth[kMax];
  int    nuccnc_truth[kMax];

  int    leppdg_truth[kMax];
  double lepp0_truth[kMax];
  double lepp1_truth[kMax];
  double lepp2_truth[kMax];
  double lepeng_truth[kMax];

  int    mutrk_contain;

  int    evthaspripf;
  double recoVertex_X;
  double recoVertex_Y;
  double recoVertex_Z;
  double recoVtxWire[3];
  double recoVtxTime[3];

  double depositU;
  double depositV;
  double depositZ;
  double ChargeU;
  double ChargeV;
  double ChargeZ;
  double correctedChargeU;
  double correctedChargeV;
  double correctedChargeZ;
  double EnergyU;
  double EnergyV;
  double EnergyZ;
  double correctedEnergyU;
  double correctedEnergyV;
  double correctedEnergyZ;

  double ErecoNu;
  double RecoLepEnNu;
  double RecoHadEnNu;
  int    RecoMethodNu;

  float  CVN_vertex[3];
  float  CVN_vertex1ststep[3];
  float  CVN_vertex2ndstep[3];

  float  CVN_vertex2[3];

  double vertexDetectorDist;
  int    nhits;
  int    ntpcs;

  int    hit_cryst      [kMaxHits];
  int    hit_tpc        [kMaxHits];
  int    hit_plane      [kMaxHits];
  int    hit_wire       [kMaxHits];
  int    hit_channel    [kMaxHits];
  int    hit_multi      [kMaxHits];
  float  hit_peakT      [kMaxHits];
  float  hit_amp        [kMaxHits];
  float  hit_minusT     [kMaxHits];
  float  hit_plusT      [kMaxHits];
  float  hit_rmsT       [kMaxHits];
  float  hit_startT     [kMaxHits];
  float  hit_endT       [kMaxHits];
  float  hit_charge     [kMaxHits];
  float  hit_corrcharge [kMaxHits];
  int    hit_truetrackid[kMaxHits];
  double hit_energy     [kMaxHits];
  float  hit_good       [kMaxHits];
  int    hit_ndf        [kMaxHits];
  int    hit_prong_tag  [kMaxHits];
  int    hit_prong_pdg  [kMaxHits];
  int    hit_prong_pdg_mom [kMaxHits];

  int    hit_ntrackID   [kMaxHits];
  double hit_trueE      [kMaxHits];

  int    hit_nID_Pri    [kMaxHits];
  int    hit_PDG        [kMaxHits];
  int    hit_Mo_PDG     [kMaxHits];
  int    hit_mu_flag    [kMaxHits];
  int    hit_ele_flag    [kMaxHits];

  double hit_x          [kMaxHits];
  double hit_y          [kMaxHits];
  double hit_z          [kMaxHits];

  double hit_sp_x          [kMaxHits];
  double hit_sp_y          [kMaxHits];
  double hit_sp_z          [kMaxHits];


  int    hit_shwid      [kMaxHits];
  int    hit_ledshwid   [kMaxHits];
  double hit_shwdirX    [kMaxHits];
  double hit_shwdirY    [kMaxHits];
  double hit_shwdirZ    [kMaxHits];
  double hit_shwcosth   [kMaxHits];
  double hit_shwphi     [kMaxHits];
  int    hit_shwstartT  [kMaxHits];
  int    hit_shwendT    [kMaxHits];
  double hit_shwevteng  [kMaxHits];
  double hit_shwevtdedx [kMaxHits];

  int n_tracks_pand;
  double  all_track_length[kMax];
  double all_track_mom[kMax];
  double all_track_calmom[kMax];
  double all_track_startx[kMax];
  double all_track_starty[kMax];
  double all_track_startz[kMax];
  double all_track_px[kMax];
  double all_track_py[kMax];
  double all_track_pz[kMax];

  int all_track_true_pdg[kMax];
  int all_track_true_pdg_mom[kMax];
  double all_track_true_Efrac[kMax];
  double all_track_dist_reco_vtx[kMax];
  double all_track_dist_true_vtx[kMax];
  float  all_track_min_angle[kMax];
  double all_track_true_eng[kMax];
  double all_track_true_px[kMax];
  double all_track_true_py[kMax];
  double all_track_true_pz[kMax];
  double all_track_true_mom_px[kMax];
  double all_track_true_mom_py[kMax];
  double all_track_true_mom_pz[kMax];

  double all_track_true_endx[kMax];
  double all_track_true_endy[kMax];
  double all_track_true_endz[kMax];

  double all_track_endx[kMax];
  double all_track_endy[kMax];
  double all_track_endz[kMax];

  int all_track2_true_pdg[kMax];
  int all_track2_true_pdg_mom[kMax];
  double all_track2_true_Efrac[kMax];
  double all_track2_true_eng[kMax];
  double all_track2_true_px[kMax];
  double all_track2_true_py[kMax];
  double all_track2_true_pz[kMax];
  double all_track2_true_mom_px[kMax];
  double all_track2_true_mom_py[kMax];
  double all_track2_true_mom_pz[kMax];


  double shw_ChargeU;
  double shw_ChargeV;
  double shw_ChargeZ;
  double shw_correctedChargeU;
  double shw_correctedChargeV;
  double shw_correctedChargeZ;
  double shw_EnergyU;
  double shw_EnergyV;
  double shw_EnergyZ;
  double shw_correctedEnergyU;
  double shw_correctedEnergyV;
  double shw_correctedEnergyZ;

  int n_showers;
  float shw_tot_chg;
  double shw_chg[kMax];
  int   shw_nhits[kMax];

  int    led_shw_idx;
  int    led_shw_id;
  double led_shw_eng;
  double led_shw_dirx;
  double led_shw_diry;
  double led_shw_dirz;
  double led_shw_startx;
  double led_shw_starty;
  double led_shw_startz;
  double led_shw_costh;
  double shwVtxTPC[3];
  double shwVtxWire[3];
  double shwVtxTime[3];
  double shwVtxGlobalWire[3];
  double shwVtxGlobalTick[3];
  double shwVtxGlobalPlane[3];

  double led_shw_trueP0;
  double led_shw_trueP1;
  double led_shw_trueP2;
  double led_shw_trueP3;
  double led_shw_truedang;
  int    led_shw_truePDG;
  int    led_shw_truePDGMom;
  double led_shw_truecosth;

  float  shw_first_wire;
  float  shw_first_time;

  int all_shw_idx[kMax];
  double all_shw_eng[kMax];
  double all_shw_dedx[kMax];
  float all_shw_length[kMax];
  double all_shw_startx[kMax];
  double all_shw_starty[kMax];
  double all_shw_startz[kMax];
  double all_shw_px[kMax];
  double all_shw_py[kMax];
  double all_shw_pz[kMax];

  double all_shw_costh[kMax];
  double all_shw_phi[kMax];
  double all_shw_true_Efrac[kMax];
  int all_shw_true_pdg[kMax];
  int all_shw_true_pdg_mom[kMax];
  double all_shw_dist_reco_vtx[kMax];
  double all_shw_dist_true_vtx[kMax];
  float  all_shw_min_angle[kMax];
  double all_shw_true_eng[kMax];
  double all_shw_true_px[kMax];
  double all_shw_true_py[kMax];
  double all_shw_true_pz[kMax];
  double all_shw_true_mom_px[kMax];
  double all_shw_true_mom_py[kMax];
  double all_shw_true_mom_pz[kMax];

  double all_shw2_true_Efrac[kMax];
  int all_shw2_true_pdg[kMax];
  int all_shw2_true_pdg_mom[kMax];
  double all_shw2_true_eng[kMax];
  double all_shw2_true_px[kMax];
  double all_shw2_true_py[kMax];
  double all_shw2_true_pz[kMax];
  double all_shw2_true_mom_px[kMax];
  double all_shw2_true_mom_py[kMax];
  double all_shw2_true_mom_pz[kMax];


  int n_tracks_pmtrack;

  int n_pand_particles;
  int n_prims;
  int m_pf_vtx;
  double m_pf_vtx_x;
  double m_pf_vtx_y;
  double m_pf_vtx_z;
  double PFVtxTPC[3];
  double PFVtxWire[3];
  double PFVtxTime[3];
  double PFVtxGlobalWire[3];
  double PFVtxGlobalTick[3];
  double PFVtxGlobalPlane[3];

  int clst_lead_idx;
  int n_clusters;
  int clst_plane[kMax];
  float clst_integral[kMax];
  float clst_sumadc[kMax];
  float clst_width[kMax];
  int clst_nhits[kMax];
  float clst_startwire[kMax];
  float clst_starttick[kMax];
  float clst_endwire[kMax];
  float clst_endtick[kMax];

  int hit_shw_count;

  int shw_idx[kMaxHits];
  double shw_startX[kMaxHits];
  double shw_startY[kMaxHits];
  double shw_startZ[kMaxHits];
  double shw_dirX[kMaxHits];
  double shw_dirY[kMaxHits];
  double shw_dirZ[kMaxHits];

  int hit_shw_tpc[kMaxHits];
  int hit_shw_plane[kMaxHits];
  int hit_shw_wire[kMaxHits];
  int hit_shw_channel[kMaxHits];
  int hit_shw_multi[kMaxHits];
  int hit_shw_truetrackid[kMaxHits];
  int hit_shw_ntrackID[kMaxHits];
  float hit_shw_peakT[kMaxHits];
  int hit_shw_startT[kMaxHits];
  int hit_shw_endT[kMaxHits];
  double hit_shw_charge[kMaxHits];
  double hit_shw_corrcharge[kMaxHits];
  double hit_shw_energy[kMaxHits];
  double hit_shw_trueE[kMaxHits];

  int hit_shw_pdg[kMaxHits];
  int hit_shw_pdg_mom[kMaxHits];
  double hit_shw_global_plane[kMaxHits];
  double hit_shw_global_wire[kMaxHits];
  double hit_shw_global_tick[kMaxHits];

  int hit_Nallraws[kMaxHits];
  int hit_Nadcvec[kMaxHits];
  int hit_raw_wireID[kMaxHits];
  int hit_raw_chID[kMaxHits];
  std::vector<std::vector<int> > hit_rawadc;
  std::vector<std::vector<int> > hit_rawtick;

  // make new wire index for the pixel map

  int n_wires;
  int wire_assn[kMaxHits];
  int wire_min_tick[kMaxHits];
  int wire_max_tick[kMaxHits];
  int wire_min_shwtick[kMaxHits][kMaxClst];
  int wire_max_shwtick[kMaxHits][kMaxClst];
  int wire_min_trktick[kMaxHits][kMaxClst];
  int wire_max_trktick[kMaxHits][kMaxClst];
  int wire_max_shwrms[kMaxHits][kMaxClst];
  int wire_max_trkrms[kMaxHits][kMaxClst];


  int wire_min_ledshwtick[kMaxHits];
  int wire_max_ledshwtick[kMaxHits];
  float wire_max_rmsT[kMaxHits];


  std::vector<std::vector<double> > wire_adc;
  std::vector<std::vector<double> > wire_corradc;
  std::vector<std::vector<int> > wire_tick;
  std::vector<std::vector<double> > wire_x;
  std::vector<std::vector<int> > wire_shwflag;
  std::vector<std::vector<int> > wire_ledshwflag;
  std::vector<std::vector<int> > wire_pdg;
  std::vector<std::vector<int> > wire_pdg_mom;
  std::vector<std::vector<double> > wire_dr;
  std::vector<std::vector<double> > wire_dist;
  std::vector<std::vector<int> > wire_tag_pdg;
  std::vector<std::vector<int> > wire_tag_pdg_mom;
  std::vector<std::vector<int> > wire_tag_is_shw;
  std::vector<std::vector<int> > wire_tag_ith;
  std::vector<std::vector<double> > wire_tag_true_px;
  std::vector<std::vector<double> > wire_tag_true_py;
  std::vector<std::vector<double> > wire_tag_true_pz;
  std::vector<std::vector<double> > wire_tag_reco_px;
  std::vector<std::vector<double> > wire_tag_reco_py;
  std::vector<std::vector<double> > wire_tag_reco_pz;
  std::vector<std::vector<double> > wire_tag_reco_eng;

  std::vector<std::vector<double> > ch_chg;
  std::vector<std::vector<int> > ch_tick;

  double hit_global_wire[kMaxHits];
  double hit_global_tick[kMaxHits];
  double hit_global_plane[kMaxHits];
  double hit_offset[3];

  int wire_roifirst[kMaxHits];
  int wire_roiend[kMaxHits];
//  float wirecharge;
//  int wiretick;
  double sum_ch_chg[kMaxHits];

  int leadingShowerID;

  std::vector<int> trueDaughterPDGs;
  std::vector<float> trueDaughterE;
  std::vector<float> trueDaughterDepoE;
  std::vector<float> trueDaughterPx;
  std::vector<float> trueDaughterPy;
  std::vector<float> trueDaughterPz;


//  int hit_IDEs_size[kMaxHits];

  int TrueParticle_size;
  std::vector<int> TruePart_PDG;
  std::vector<double> TruePart_E;

  int simCh_size;
  int tdc_size[kMaxHits];
  std::vector<double> ChgU;
  std::vector<double> ChgV;
  std::vector<double> ChgZ;

  TVector3 *v3_on_planes[3];

  float max_shwtick[kMaxClst];
  float min_shwtick[kMaxClst];
  float max_trktick[kMaxClst];
  float min_trktick[kMaxClst];
  float max_shwrms[kMaxClst];
  float max_trkrms[kMaxClst];

  std::string fHitsModuleLabel;
  std::string fClusterModuleLabel;
  std::string fShowerModuleLabel;
  std::string fTrackModuleLabel1;
  std::string fTrackModuleLabel2;
  std::string fRawDigitModuleLabel;
  std::string fVertexModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fPandoraNuVertexModuleLabel;
  std::string fMCGenModuleLabel;
  std::string fEnergyRecoNueLabel;
  std::string fEnergyRecoNumuLabel;
  std::string fEnergyRecoNuelseLabel;
  /*std::string fRegCNNModuleLabel;
  std::string fRegCNNResultLabel;

  std::string fRegCNN1ststepModuleLabel;
  std::string fRegCNN1ststepResultLabel;
  std::string fRegCNN2ndstepModuleLabel;
  std::string fRegCNN2ndstepResultLabel;

  std::string fRegCNNModule2Label;
  std::string fRegCNNResult2Label;*/
  bool fIsVD;

  int fGlobalWireMethod;

  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  art::ServiceHandle<geo::Geometry> geom;
  geo::Geometry* m_geom = &*geom;
  geo::WireReadoutGeom const* wireReadout = &art::ServiceHandle<geo::WireReadout>()->Get();
  detinfo::DetectorProperties const* detprop = nullptr;

  double fCVNResultIsAntineutrino;
  double fCVNResultNue, fCVNResultNumu, fCVNResultNutau, fCVNResultNC; // flavour
  double fCVNResult0Protons, fCVNResult1Protons, fCVNResult2Protons, fCVNResultNProtons; // #protons
  double fCVNResult0Pions, fCVNResult1Pions, fCVNResult2Pions, fCVNResultNPions; // #pions
  double fCVNResult0Pizeros, fCVNResult1Pizeros, fCVNResult2Pizeros, fCVNResultNPizeros; // #pizeros
  double fCVNResult0Neutrons, fCVNResult1Neutrons, fCVNResult2Neutrons, fCVNResultNNeutrons; // #neutrons

};

recoestudy::RecoEnergyS::RecoEnergyS(fhicl::ParameterSet const& pset) : EDAnalyzer(pset),
	fCaloAlg (pset.get<fhicl::ParameterSet>("CalorimetryAlg")){
std::cout<<"recoenergy RecoEnergyS::RecoEnergyS(fhicl::ParameterSet const& pset) "<<std::endl;
  this->reconfigure(pset);

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("WC","WC");
  fTree->Branch("ievent",            &ievt, "ievent/I");
  fTree->Branch("run",               &run, "run/I");
  fTree->Branch("subrun",            &subrun, "subrun/I");
  fTree->Branch("det_t0",            &det_t0, "det_t0/D");
  fTree->Branch("driftvelocity",     &driftvelocity, "driftvelocity/D");
  fTree->Branch("SampleRate",        &SampleRate, "SampleRate/D");
  fTree->Branch("NTPCs_evt",         &ntpcs_evt, "NTPCs_evt/I");
  fTree->Branch("NWires_P0",         &nwires_p0, "NWires_P0/I");
  fTree->Branch("NWires_P1",         &nwires_p1, "NWires_P1/I");
  fTree->Branch("NWires_P2",         &nwires_p2, "NWires_P2/I");

  fTree->Branch("NuTruthFid",        &nutrue_fid, "NuTruthFid/I");
  fTree->Branch("NuTruthN",          &nu_truth_N, "NuTruthN/I");
  fTree->Branch("NuEngTruth",        nueng_truth, "NuEngTruth[NuTruthN]/D");
  fTree->Branch("NuPDGTruth",        nupdg_truth, "NuPDGTruth[NuTruthN]/I");
  fTree->Branch("NuVtxXTruth",       nuvtxx_truth, "NuVtxXTruth[NuTruthN]/D");
  fTree->Branch("NuVtxYTruth",       nuvtxy_truth, "NuVtxYTruth[NuTruthN]/D");
  fTree->Branch("NuVtxZTruth",       nuvtxz_truth, "NuVtxZTruth[NuTruthN]/D");
  fTree->Branch("NuModeTruth",       numode_truth, "NuModeTruth[NuTruthN]/I");
  fTree->Branch("NuCCNCTruth",       nuccnc_truth, "NuCCNCTruth[NuTruthN]/I");

  fTree->Branch("LepPDGTruth",       leppdg_truth, "LepPDGTruth[NuTruthN]/I");
  fTree->Branch("LepP0Truth",        lepp0_truth,  "LepP0Truth[NuTruthN]/D");
  fTree->Branch("LepP1Truth",        lepp1_truth,  "LepP1Truth[NuTruthN]/D");
  fTree->Branch("LepP2Truth",        lepp2_truth,  "LepP2Truth[NuTruthN]/D");
  fTree->Branch("LepEngTruth",       lepeng_truth,  "LepEngTruth[NuTruthN]/D");
  fTree->Branch("MuTrackCont",       &mutrk_contain, "MuTrackCont/I");

  fTree->Branch("PrimPDG",           &primpdg,   "PrimPDG/I");
  fTree->Branch("TrueEnergy",        &trueEnergy,   "TrueEnergy/D");
  fTree->Branch("TrueVertexX",       &trueVertex_X, "TrueVertexX/D");
  fTree->Branch("TrueVertexY",       &trueVertex_Y, "TrueVertexY/D");
  fTree->Branch("TrueVertexZ",       &trueVertex_Z, "TrueVertexZ/D");
  fTree->Branch("TrueEndX",          &trueEnd_X, "TrueEndX/D");
  fTree->Branch("TrueEndY",          &trueEnd_Y, "TrueEndY/D");
  fTree->Branch("TrueEndZ",          &trueEnd_Z, "TrueEndZ/D");
  fTree->Branch("TruePx",            &truePx, "TruePx/D");
  fTree->Branch("TruePy",            &truePy, "TruePy/D");
  fTree->Branch("TruePz",            &truePz, "TruePz/D");
  fTree->Branch("TrueVtxin",         &vtx_in_tpc, "TrueVtxin/I");
  fTree->Branch("TrueVtxWire",       trueVtxWire, "TrueVtxWire[3]/D");
  fTree->Branch("TrueVtxGlobalWire", trueVtxGlobalWire, "TrueVtxGlobalWire[3]/D");
  fTree->Branch("TrueVtxGlobalTick", trueVtxGlobalTick, "TrueVtxGlobalTick[3]/D");
  fTree->Branch("TrueVtxGlobalPlane", trueVtxGlobalPlane, "TrueVtxGlobalPlane[3]/D");
  fTree->Branch("TrueVtxTPC",        trueVtxTPC,   "TrueVtxTPC[3]/D");
  fTree->Branch("TrueVtxTime",       trueVtxTime,  "TrueVtxTime[3]/D");
  fTree->Branch("TrueVtxChan",       trueVtxChan,  "TrueVtxChan[3]/D");
  fTree->Branch("TrueEndWire",       trueEndWire,  "TrueEndWire[3]/D");
  fTree->Branch("TrueEndTime",       trueEndTime,  "TrueEndTime[3]/D");
  fTree->Branch("trueEndGlobalWire",       trueEndGlobalWire,  "trueEndGlobalWire[3]/D");
  fTree->Branch("trueEndGlobalTime",       trueEndGlobalTime,  "trueEndGlobalTime[3]/D");
  fTree->Branch("TrueMidWire",       trueMidWire,  "TrueMidWire[3]/D");
  fTree->Branch("TrueMidTime",       trueMidTime,  "TrueMidTime[3]/D");

  fTree->Branch("CVNVertex1ststep",  CVN_vertex1ststep,  "CVNVertex1ststep[3]/F");
  fTree->Branch("CVNVertex2ndstep",  CVN_vertex2ndstep,  "CVNVertex2ndstep[3]/F");

  fTree->Branch("CVNVertex2",        CVN_vertex2,  "CVNVertex2[3]/F");
  fTree->Branch("CVNVtx2Wire",       CVNVtx2Wire, "CVNVtx2Wire[3]/D");
  fTree->Branch("CVNVtx2TPC",        CVNVtx2TPC,   "CVNVtx2TPC[3]/D");
  fTree->Branch("CVNVtx2Time",       CVNVtx2Time,  "CVNVtx2Time[3]/D");
  fTree->Branch("CVNVtx2GlobalWire", CVNVtx2GlobalWire, "CVNVtx2GlobalWire[3]/D");
  fTree->Branch("CVNVtx2GlobalTick", CVNVtx2GlobalTick, "CVNVtx2GlobalTick[3]/D");
  fTree->Branch("CVNVtx2GlobalPlane",CVNVtx2GlobalPlane, "CVNVtx2GlobalPlane[3]/D");


  fTree->Branch("CVNVertex",         CVN_vertex,  "CVNVertex[3]/F");
  fTree->Branch("CVNVtxWire",        CVNVtxWire, "CVNVtxWire[3]/D");
  fTree->Branch("CVNVtxTPC",         CVNVtxTPC,   "CVNVtxTPC[3]/D");
  fTree->Branch("CVNVtxTime",        CVNVtxTime,  "CVNVtxTime[3]/D");
  fTree->Branch("CVNVtxGlobalWire",  CVNVtxGlobalWire, "CVNVtxGlobalWire[3]/D");
  fTree->Branch("CVNVtxGlobalTick",  CVNVtxGlobalTick, "CVNVtxGlobalTick[3]/D");
  fTree->Branch("CVNVtxGlobalPlane", CVNVtxGlobalPlane, "CVNVtxGlobalPlane[3]/D");

  fTree->Branch("EndProngGlobalWire", EndProngGlobalWire, "EndProngGlobalWire[3]/D");
  fTree->Branch("EndProngGlobalTick", EndProngGlobalTick, "EndProngGlobalTick[3]/D");
  fTree->Branch("EndProngGlobalPlane", EndProngGlobalPlane, "EndProngGlobalPlane[3]/D");


  fTree->Branch("EvthasPriPF",       &evthaspripf, "EvthasPriPF/I");
  fTree->Branch("RecoVertexX",       &recoVertex_X, "RecoVertexX/D");
  fTree->Branch("RecoVertexY",       &recoVertex_Y, "RecoVertexY/D");
  fTree->Branch("RecoVertexZ",       &recoVertex_Z, "RecoVertexZ/D");
  fTree->Branch("RecoVtxWire",       recoVtxWire,   "RecoVtxWire[3]/D");
  fTree->Branch("RecoVtxTime",       recoVtxTime,   "RecoVtxTime[3]/D");


  fTree->Branch("DepositU",          &depositU, "DepositU/D");
  fTree->Branch("DepositV",          &depositV, "DepositV/D");
  fTree->Branch("DepositZ",          &depositZ, "DepositZ/D");
  fTree->Branch("ChargeU",           &ChargeU, "ChargeU/D");
  fTree->Branch("ChargeV",           &ChargeV, "ChargeV/D");
  fTree->Branch("ChargeZ",           &ChargeZ, "ChargeZ/D");
  fTree->Branch("CorrectedChargeU",  &correctedChargeU, "CorrectedChargeU/D");
  fTree->Branch("CorrectedChargeV",  &correctedChargeV, "CorrectedChargeV/D");
  fTree->Branch("CorrectedChargeZ",  &correctedChargeZ, "CorrectedChargeZ/D");
  fTree->Branch("VertexDetectorDist",&vertexDetectorDist, "VertexDetectorDist/D");
  fTree->Branch("EnergyU",           &EnergyU, "EnergyU/D");
  fTree->Branch("EnergyV",           &EnergyV, "EnergyV/D");
  fTree->Branch("EnergyZ",           &EnergyZ, "EnergyZ/D");
  fTree->Branch("CorrectedEnergyU",  &correctedEnergyU, "CorrectedEnergyU/D");
  fTree->Branch("CorrectedEnergyV",  &correctedEnergyV, "CorrectedEnergyV/D");
  fTree->Branch("CorrectedEnergyZ",  &correctedEnergyZ, "CorrectedEnergyZ/D");

  fTree->Branch("ErecoNu",           &ErecoNu,      "ErecoNu/D");
  fTree->Branch("RecoLepEnNu",       &RecoLepEnNu,  "RecoLepEnNu/D");
  fTree->Branch("RecoHadEnNu",       &RecoHadEnNu,  "RecoHadEnNu/D");
  fTree->Branch("RecoMethodNu",      &RecoMethodNu, "RecoMethodNu/I");

  fTree->Branch("NHits",             &nhits, "NHits/I");
  fTree->Branch("NTPCs",             &ntpcs, "NTPCs/I");

  fTree->Branch("Hit_Cryst",         hit_cryst,       "Hit_Cryst[NHits]/I");
  fTree->Branch("Hit_TPC",           hit_tpc,         "Hit_TPC[NHits]/I");
  fTree->Branch("Hit_Plane",         hit_plane,       "Hit_Plane[NHits]/I");
  fTree->Branch("Hit_Wire",          hit_wire,        "Hit_Wire[NHits]/I");
  fTree->Branch("Hit_Channel",       hit_channel,     "Hit_Channel[NHits]/I");
  fTree->Branch("Hit_Multi",         hit_multi,       "Hit_Multi[NHits]/I");
  fTree->Branch("Hit_PeakT",         hit_peakT,       "Hit_PeakT[NHits]/F");
  fTree->Branch("Hit_Amp",           hit_amp,         "Hit_Amp[NHits]/F");
  fTree->Branch("Hit_MinusT",        hit_minusT,      "Hit_MinusT[NHits]/F");
  fTree->Branch("Hit_PlusT",         hit_plusT,       "Hit_PlusT[NHits]/F");
  fTree->Branch("Hit_RMST",          hit_rmsT,        "Hit_RMST[NHits]/F");
  fTree->Branch("Hit_StartT",        hit_startT,      "Hit_StartT[NHits]/F");
  fTree->Branch("Hit_EndT",          hit_endT,        "Hit_EndT[NHits]/F");
  fTree->Branch("Hit_Charge",        hit_charge,      "Hit_Charge[NHits]/F");
  fTree->Branch("Hit_CorrCharge",    hit_corrcharge,  "Hit_CorrCharge[NHits]/F");
  fTree->Branch("Hit_Energy",        hit_energy,      "Hit_Energy[NHits]/D");
  fTree->Branch("Hit_TrueTrackID",   hit_truetrackid, "Hit_TrueTrackID[NHits]/I");
  fTree->Branch("Hit_X",             hit_x,           "Hit_X[NHits]/D");
  fTree->Branch("Hit_Y",             hit_y,           "Hit_Y[NHits]/D");
  fTree->Branch("Hit_Z",             hit_z,           "Hit_Z[NHits]/D");
  fTree->Branch("Hit_SP_X",             hit_sp_x,           "Hit_SP_X[NHits]/D");
  fTree->Branch("Hit_SP_Y",             hit_sp_y,           "Hit_SP_Y[NHits]/D");
  fTree->Branch("Hit_SP_Z",             hit_sp_z,           "Hit_SP_Z[NHits]/D");
  fTree->Branch("Hit_Global_Wire",   hit_global_wire, "Hit_Global_Wire[NHits]/D");
  fTree->Branch("Hit_Global_Tick",   hit_global_tick, "Hit_Global_Tick[NHits]/D");
  fTree->Branch("Hit_Global_Plane",  hit_global_plane, "Hit_Global_Plane[NHits]/D");
  fTree->Branch("Hit_Good",          hit_good,        "Hit_Good[NHits]/F");
  fTree->Branch("Hit_Ndf",           hit_ndf,         "Hit_Ndf[NHits]/I");
  fTree->Branch("Hit_Prong_Tag",     hit_prong_tag,   "Hit_Prong_Tag[NHits]/I");
  fTree->Branch("Hit_Prong_PDG",     hit_prong_pdg,   "Hit_Prong_PDG[NHits]/I");
  fTree->Branch("Hit_Prong_PDG_Mom", hit_prong_pdg_mom, "Hit_Prong_PDG_Mom[NHits]/I");

  fTree->Branch("Hit_ShwID",         hit_shwid,       "Hit_ShwID[NHits]/I");
  fTree->Branch("Hit_LedShwID",      hit_ledshwid,    "Hit_LedShwID[NHits]/I");
  fTree->Branch("Hit_ShwDirX",       hit_shwdirX,     "Hit_ShwDirX[NHits]/D");
  fTree->Branch("Hit_ShwDirY",       hit_shwdirY,     "Hit_ShwDirY[NHits]/D");
  fTree->Branch("Hit_ShwDirZ",       hit_shwdirZ,     "Hit_ShwDirZ[NHits]/D");
  fTree->Branch("Hit_ShwCoth",       hit_shwcosth,    "Hit_ShwCoth[NHits]/D");
  fTree->Branch("Hit_ShwPhi",        hit_shwphi,      "Hit_ShwPhi[NHits]/D");
  fTree->Branch("Hit_ShwStartT",     hit_shwstartT,   "Hit_ShwStartT[NHits]/I");
  fTree->Branch("Hit_ShwEndT",       hit_shwendT,     "Hit_ShwEndT[NHits]/I");
  fTree->Branch("Hit_ShwEvtEng",     hit_shwevteng,   "Hit_ShwEvtEng[NHits]/D");
  fTree->Branch("Hit_ShwEvtdEdx",    hit_shwevtdedx,  "Hit_ShwEvtdEdx[NHits]/D");


  fTree->Branch("Hit_ntrackID",      hit_ntrackID,   "Hit_ntrackID[NHits]/I");
  fTree->Branch("Hit_trueE",         hit_trueE,      "Hit_trueE[NHits]/D");

  fTree->Branch("Hit_nID_Pri",       hit_nID_Pri,    "Hit_nID_Pri[NHits]/I");
  fTree->Branch("Hit_PDG",           hit_PDG,        "Hit_PDG[NHits]/I");
  fTree->Branch("Hit_Mo_PDG",        hit_Mo_PDG,     "Hit_Mo_PDG[NHits]/I");
  fTree->Branch("Hit_mu_flag",       hit_mu_flag,    "Hit_mu_flag[NHits]/I");
  fTree->Branch("Hit_ele_flag",      hit_ele_flag,    "Hit_ele_flag[NHits]/I");

  fTree->Branch("Hit_offset",        hit_offset,    "Hit_offset[3]/D");

  fTree->Branch("n_tracks_pad",       &n_tracks_pand, "n_tracks_pad/I");
  fTree->Branch("all_track_length",   all_track_length,  "all_track_length[n_tracks_pad]/D");
  fTree->Branch("all_track_mom",      all_track_mom,  "all_track_mom[n_tracks_pad]/D");
  fTree->Branch("all_track_calmom",   all_track_calmom,  "all_track_calmom[n_tracks_pad]/D");
  fTree->Branch("all_track_startx",   all_track_startx,  "all_track_startx[n_tracks_pad]/D");
  fTree->Branch("all_track_starty",   all_track_starty,  "all_track_starty[n_tracks_pad]/D");
  fTree->Branch("all_track_startz",   all_track_startz,  "all_track_startz[n_tracks_pad]/D");
  fTree->Branch("all_track_px",       all_track_px,      "all_track_px[n_tracks_pad]/D");
  fTree->Branch("all_track_py",       all_track_py,      "all_track_py[n_tracks_pad]/D");
  fTree->Branch("all_track_pz",       all_track_pz,      "all_track_pz[n_tracks_pad]/D");
  fTree->Branch("all_track_true_Efrac",     all_track_true_Efrac,     "all_track_true_Efrac[n_tracks_pad]/D");
  fTree->Branch("all_track_true_pdg",       all_track_true_pdg,     "all_track_true_pdg[n_tracks_pad]/I");
  fTree->Branch("all_track_true_pdg_mom",   all_track_true_pdg_mom,     "all_track_true_pdg_mom[n_tracks_pad]/I");
  fTree->Branch("all_track_dist_reco_vtx",  all_track_dist_reco_vtx,     "all_track_dist_reco_vtx[n_tracks_pad]/D");
  fTree->Branch("all_track_dist_true_vtx",  all_track_dist_true_vtx,     "all_track_dist_true_vtx[n_tracks_pad]/D");
  fTree->Branch("all_track_min_angle",      all_track_min_angle,     "all_track_min_angle[n_tracks_pad]/F");
  fTree->Branch("all_track_true_eng",       all_track_true_eng,     "all_track_true_eng[n_tracks_pad]/D");
  fTree->Branch("all_track_true_px",        all_track_true_px,      "all_track_true_px[n_tracks_pad]/D");
  fTree->Branch("all_track_true_py",        all_track_true_py,      "all_track_true_py[n_tracks_pad]/D");
  fTree->Branch("all_track_true_pz",        all_track_true_pz,      "all_track_true_pz[n_tracks_pad]/D");
  fTree->Branch("all_track_true_mom_px",        all_track_true_mom_px,      "all_track_true_mom_px[n_tracks_pad]/D");
  fTree->Branch("all_track_true_mom_py",        all_track_true_mom_py,      "all_track_true_mom_py[n_tracks_pad]/D");
  fTree->Branch("all_track_true_mom_pz",        all_track_true_mom_pz,      "all_track_true_mom_pz[n_tracks_pad]/D");

  fTree->Branch("all_track_true_endx",        all_track_true_endx,      "all_track_true_endx[n_tracks_pad]/D");
  fTree->Branch("all_track_true_endy",        all_track_true_endy,      "all_track_true_endy[n_tracks_pad]/D");
  fTree->Branch("all_track_true_endz",        all_track_true_endz,      "all_track_true_endz[n_tracks_pad]/D");

  fTree->Branch("all_track_endx",        all_track_endx,      "all_track_endx[n_tracks_pad]/D");
  fTree->Branch("all_track_endy",        all_track_endy,      "all_track_endy[n_tracks_pad]/D");
  fTree->Branch("all_track_endz",        all_track_endz,      "all_track_endz[n_tracks_pad]/D");

  fTree->Branch("all_track2_true_Efrac",     all_track2_true_Efrac,   "all_track2_true_Efrac[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_pdg",       all_track2_true_pdg,     "all_track2_true_pdg[n_tracks_pad]/I");
  fTree->Branch("all_track2_true_pdg_mom",   all_track2_true_pdg_mom, "all_track2_true_pdg_mom[n_tracks_pad]/I");
  fTree->Branch("all_track2_true_eng",       all_track2_true_eng,   "all_track2_true_eng[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_px",        all_track2_true_px,      "all_track2_true_px[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_py",        all_track2_true_py,      "all_track2_true_py[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_pz",        all_track2_true_pz,      "all_track2_true_pz[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_mom_px",    all_track2_true_mom_px,  "all_track2_true_mom_px[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_mom_py",    all_track2_true_mom_py,  "all_track2_true_mom_py[n_tracks_pad]/D");
  fTree->Branch("all_track2_true_mom_pz",    all_track2_true_mom_pz,  "all_track2_true_mom_pz[n_tracks_pad]/D");


  fTree->Branch("n_showers",         &n_showers,     "n_showers/I");
  fTree->Branch("shw_tot_chg",       &shw_tot_chg,   "shw_tot_chg/F");
  fTree->Branch("shw_chg",           shw_chg,        "shw_chg[n_showers]/D");
  fTree->Branch("shw_nhits",         shw_nhits,      "shw_nhits[n_showers]/I");
  fTree->Branch("all_shw_idx",       all_shw_idx,    "all_shw_idx[n_showers]/I");
  fTree->Branch("all_shw_eng",       all_shw_eng,    "all_shw_eng[n_showers]/D");
  fTree->Branch("all_shw_dedx",      all_shw_dedx,   "all_shw_dedx[n_showers]/D");
  fTree->Branch("all_shw_length",    all_shw_length,  "all_shw_length[n_showers]/F");
  fTree->Branch("all_shw_startx",    all_shw_startx,  "all_shw_startx[n_showers]/D");
  fTree->Branch("all_shw_starty",    all_shw_starty,  "all_shw_starty[n_showers]/D");
  fTree->Branch("all_shw_startz",    all_shw_startz,  "all_shw_startz[n_showers]/D");
  fTree->Branch("all_shw_px",        all_shw_px,      "all_shw_px[n_showers]/D");
  fTree->Branch("all_shw_py",        all_shw_py,      "all_shw_py[n_showers]/D");
  fTree->Branch("all_shw_pz",        all_shw_pz,      "all_shw_pz[n_showers]/D");
  fTree->Branch("all_shw_true_eng",  all_shw_true_eng,     "all_shw_true_eng[n_showers]/D");
  fTree->Branch("all_shw_true_px",   all_shw_true_px,      "all_shw_true_px[n_showers]/D");
  fTree->Branch("all_shw_true_py",   all_shw_true_py,      "all_shw_true_py[n_showers]/D");
  fTree->Branch("all_shw_true_pz",   all_shw_true_pz,      "all_shw_true_pz[n_showers]/D");
  fTree->Branch("all_shw_true_mom_px",   all_shw_true_mom_px,      "all_shw_true_mom_px[n_showers]/D");
  fTree->Branch("all_shw_true_mom_py",   all_shw_true_mom_py,      "all_shw_true_mom_py[n_showers]/D");
  fTree->Branch("all_shw_true_mom_pz",   all_shw_true_mom_pz,      "all_shw_true_mom_pz[n_showers]/D");

  fTree->Branch("all_shw_costh",     all_shw_costh,   "all_shw_costh[n_showers]/D");
  fTree->Branch("all_shw_phi",       all_shw_phi,     "all_shw_phi[n_showers]/D");
  fTree->Branch("all_shw_true_Efrac",     all_shw_true_Efrac,     "all_shw_true_Efrac[n_showers]/D");
  fTree->Branch("all_shw_true_pdg",       all_shw_true_pdg,     "all_shw_true_pdg[n_showers]/I");
  fTree->Branch("all_shw_true_pdg_mom",   all_shw_true_pdg_mom,     "all_shw_true_pdg_mom[n_showers]/I");
  fTree->Branch("all_shw_dist_reco_vtx",  all_shw_dist_reco_vtx,     "all_shw_dist_reco_vtx[n_showers]/D");
  fTree->Branch("all_shw_dist_true_vtx",  all_shw_dist_true_vtx,     "all_shw_dist_true_vtx[n_showers]/D");
  fTree->Branch("all_shw_min_angle",       all_shw_min_angle,     "all_shw_min_angle[n_showers]/F");

  fTree->Branch("all_shw2_true_Efrac",     all_shw2_true_Efrac,     "all_shw2_true_Efrac[n_showers]/D");
  fTree->Branch("all_shw2_true_pdg",       all_shw2_true_pdg,       "all_shw2_true_pdg[n_showers]/I");
  fTree->Branch("all_shw2_true_pdg_mom",   all_shw2_true_pdg_mom,   "all_shw2_true_pdg_mom[n_showers]/I");
  fTree->Branch("all_shw2_true_eng",       all_shw2_true_eng,        "all_shw2_true_eng[n_showers]/D");
  fTree->Branch("all_shw2_true_px",        all_shw2_true_px,        "all_shw2_true_px[n_showers]/D");
  fTree->Branch("all_shw2_true_py",        all_shw2_true_py,        "all_shw2_true_py[n_showers]/D");
  fTree->Branch("all_shw2_true_pz",        all_shw2_true_pz,        "all_shw2_true_pz[n_showers]/D");
  fTree->Branch("all_shw2_true_mom_px",    all_shw2_true_mom_px,    "all_shw2_true_mom_px[n_showers]/D");
  fTree->Branch("all_shw2_true_mom_py",    all_shw2_true_mom_py,    "all_shw2_true_mom_py[n_showers]/D");
  fTree->Branch("all_shw2_true_mom_pz",    all_shw2_true_mom_pz,    "all_shw2_true_mom_pz[n_showers]/D");


  fTree->Branch("led_shw_idx",       &led_shw_idx, "led_shw_idx/I");
  fTree->Branch("led_shw_eng",       &led_shw_eng, "led_shw_eng/D");
  fTree->Branch("led_shw_dirx",      &led_shw_dirx, "led_shw_dirx/D");
  fTree->Branch("led_shw_diry",      &led_shw_diry, "led_shw_diry/D");
  fTree->Branch("led_shw_dirz",      &led_shw_dirz, "led_shw_dirz/D");
  fTree->Branch("led_shw_startx",    &led_shw_startx, "led_shw_startx/D");
  fTree->Branch("led_shw_starty",    &led_shw_starty, "led_shw_starty/D");
  fTree->Branch("led_shw_startz",    &led_shw_startz, "led_shw_startz/D");
  fTree->Branch("led_shw_costh",     &led_shw_costh, "led_shw_costh/D");
  fTree->Branch("ShwVtxTPC",         shwVtxTPC, "ShwVtxTPC[3]/D");
  fTree->Branch("ShwVtxWire",        shwVtxWire, "ShwVtxWire[3]/D");
  fTree->Branch("ShwVtxTime",        shwVtxTime,  "ShwVtxTime[3]/D");
  fTree->Branch("ShwVtxGlobalWire",  shwVtxGlobalWire, "ShwVtxGlobalWire[3]/D");
  fTree->Branch("ShwVtxGlobalTick",  shwVtxGlobalTick, "ShwVtxGlobalTick[3]/D");
  fTree->Branch("ShwVtxGlobalPlane", shwVtxGlobalPlane, "ShwVtxGlobalPlane[3]/D");

  fTree->Branch("led_shw_trueP0",       &led_shw_trueP0, "led_shw_trueP0/D");
  fTree->Branch("led_shw_trueP1",       &led_shw_trueP1, "led_shw_trueP1/D");
  fTree->Branch("led_shw_trueP2",       &led_shw_trueP2, "led_shw_trueP2/D");
  fTree->Branch("led_shw_trueP3",       &led_shw_trueP3, "led_shw_trueP3/D");
  fTree->Branch("led_shw_truedang",     &led_shw_truedang, "led_shw_truedang/D");
  fTree->Branch("led_shw_truePDG",      &led_shw_truePDG, "led_shw_truePDG/I");
  fTree->Branch("led_shw_truePDGMom",   &led_shw_truePDGMom, "led_shw_truePDGMom/I");
  fTree->Branch("led_shw_truecosth",    &led_shw_truecosth, "led_shw_truecosth/D");

  fTree->Branch("ShwFirstWire",       &shw_first_wire, "ShwFirstWire/F");
  fTree->Branch("ShwFirstTime",       &shw_first_time, "ShwFirstTime/F");

  fTree->Branch("n_clusters",       &n_clusters, "n_clusters/I");
  fTree->Branch("Clst_Lead_Idx",    &clst_lead_idx,  "Clst_Lead_idx/I");
  fTree->Branch("Clst_Plane",       clst_plane,      "Clst_Plane[n_clusters]/I");
  fTree->Branch("Clst_Integral",    clst_integral,   "Clst_Integral[n_clusters]/F");
  fTree->Branch("Clst_SumADC",      clst_sumadc,     "Clst_SumADC[n_clusters]/F");
  fTree->Branch("Clst_Width",       clst_width,      "Clst_Width[n_clusters]/F");
  fTree->Branch("Clst_Nhits",       clst_nhits,      "Clst_Nhits[n_clusters]/I");
  fTree->Branch("Clst_StartWire",    clst_startwire, "Clst_StartWire[n_clusters]/F");
  fTree->Branch("Clst_StartTick",    clst_starttick, "Clst_StartTick[n_clusters]/F");
  fTree->Branch("Clst_EndWire",      clst_endwire,   "Clst_EndWire[n_clusters]/F");
  fTree->Branch("Clst_EndTick",      clst_endtick,   "Clst_EndTick[n_clusters]/F");

  fTree->Branch("Shw_ChargeU",           &shw_ChargeU, "Shw_ChargeU/D");
  fTree->Branch("Shw_ChargeV",           &shw_ChargeV, "Shw_ChargeV/D");
  fTree->Branch("Shw_ChargeZ",           &shw_ChargeZ, "Shw_ChargeZ/D");
  fTree->Branch("Shw_CorrectedChargeU",  &shw_correctedChargeU, "Shw_CorrectedChargeU/D");
  fTree->Branch("Shw_CorrectedChargeV",  &shw_correctedChargeV, "Shw_CorrectedChargeV/D");
  fTree->Branch("Shw_CorrectedChargeZ",  &shw_correctedChargeZ, "Shw_CorrectedChargeZ/D");
  fTree->Branch("Shw_EnergyU",           &shw_EnergyU, "Shw_EnergyU/D");
  fTree->Branch("Shw_EnergyV",           &shw_EnergyV, "Shw_EnergyV/D");
  fTree->Branch("Shw_EnergyZ",           &shw_EnergyZ, "Shw_EnergyZ/D");
  fTree->Branch("Shw_CorrectedEnergyU",  &shw_correctedEnergyU, "Shw_CorrectedEnergyU/D");
  fTree->Branch("Shw_CorrectedEnergyV",  &shw_correctedEnergyV, "Shw_CorrectedEnergyV/D");
  fTree->Branch("Shw_CorrectedEnergyZ",  &shw_correctedEnergyZ, "Shw_CorrectedEnergyZ/D");


  fTree->Branch("NHits_Shw",             &hit_shw_count, "NHits_Shw/I");

  fTree->Branch("Shw_Idx",               shw_idx,             "Shw_Idx[NHits_Shw]/I" );
  fTree->Branch("Shw_StartX",            shw_startX,          "Shw_StartX[NHits_Shw]/D" );
  fTree->Branch("Shw_StartY",            shw_startY,          "Shw_StartY[NHits_Shw]/D" );
  fTree->Branch("Shw_StartZ",            shw_startZ,          "Shw_StartZ[NHits_Shw]/D" );
  fTree->Branch("Shw_DirX",              shw_dirX,            "Shw_DirX[NHits_Shw]/D" );
  fTree->Branch("Shw_DirY",              shw_dirY,            "Shw_DirY[NHits_Shw]/D" );
  fTree->Branch("Shw_DirZ",              shw_dirZ,            "Shw_DirZ[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_TPC",           hit_shw_tpc,         "Hit_Shw_TPC[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_Plane",         hit_shw_plane,       "Hit_Shw_Plane[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_Wire",          hit_shw_wire,        "Hit_Shw_Wire[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_Channel",       hit_shw_channel,     "Hit_Shw_Channel[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_Multi",         hit_shw_multi,       "Hit_Shw_Multi[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_TrueTrackID",   hit_shw_truetrackid, "Hit_Shw_TrueTrackID[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_ntrackID",      hit_shw_ntrackID,    "Hit_Shw_ntrackID[NHits_Shw]/I" );

  fTree->Branch("Hit_Shw_PeakT",         hit_shw_peakT,       "Hit_Shw_PeakT[NHits_Shw]/F" );
  fTree->Branch("Hit_Shw_StartT",        hit_shw_startT,      "Hit_Shw_StartT[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_EndT",          hit_shw_endT,        "Hit_Shw_EndT[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_Charge",        hit_shw_charge,      "Hit_Shw_Charge[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_CorrCharge",    hit_shw_corrcharge,  "Hit_Shw_CorrCharge[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_Energy",        hit_shw_energy,      "Hit_Shw_Energy[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_trueE",         hit_shw_trueE,       "Hit_Shw_trueE[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_PDG",           hit_shw_pdg,         "Hit_Shw_PDG[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_PDG_Mom",       hit_shw_pdg_mom,     "Hit_Shw_PDG_Mom[NHits_Shw]/I" );
  fTree->Branch("Hit_Shw_Global_Plane",  hit_shw_global_plane,      "Hit_Shw_Global_Plane[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_Global_Wire",   hit_shw_global_wire,      "Hit_Shw_Global_Wire[NHits_Shw]/D" );
  fTree->Branch("Hit_Shw_Global_Tick",   hit_shw_global_tick,      "Hit_Shw_Global_Tick[NHits_Shw]/D" );

  fTree->Branch("n_wires",           &n_wires, "n_wires/I");

  fTree->Branch("Hit_Nallraws",      hit_Nallraws,   "Hit_Nallraws[n_wires]/I");
  fTree->Branch("Hit_Nadcvec",       hit_Nadcvec,    "Hit_Nadcvec[n_wires]/I");
  fTree->Branch("Hit_raw_wireID",    hit_raw_wireID, "Hit_raw_wireID[n_wires]/I");
  fTree->Branch("Hit_raw_chID",      hit_raw_chID,   "Hit_raw_chID[n_wires]/I");
  fTree->Branch("Hit_rawadc",        &hit_rawadc);
  fTree->Branch("Hit_rawtick",       &hit_rawtick);


  fTree->Branch("Wire_assn",         wire_assn,      "Wire_assn[n_wires]/I");
  fTree->Branch("Wire_min_tick",     wire_min_tick,  "Wire_min_tick[n_wires]/I");
  fTree->Branch("Wire_max_tick",     wire_max_tick,  "Wire_max_tick[n_wires]/I");
  fTree->Branch("Wire_max_rmsT",     wire_max_rmsT,  "Wire_max_rmsT[n_wires]/F");
  fTree->Branch("Wire_min_shwtick",    wire_min_shwtick,  "Wire_min_shwtick[n_wires][100]/I");
  fTree->Branch("Wire_max_shwtick",    wire_max_shwtick,  "Wire_max_shwtick[n_wires][100]/I");
  fTree->Branch("Wire_min_trktick",    wire_min_trktick,  "Wire_min_trktick[n_wires][100]/I");
  fTree->Branch("Wire_max_trktick",    wire_max_trktick,  "Wire_max_trktick[n_wires][100]/I");
  fTree->Branch("Wire_max_shwrms",     wire_max_shwrms,  "Wire_max_shwrms[n_wires][100]/I");
  fTree->Branch("Wire_max_trkrms",     wire_max_trkrms,  "Wire_max_trkrms[n_wires][100]/I");

  fTree->Branch("Wire_min_ledshwtick", wire_min_ledshwtick,  "Wire_min_ledshwtick[n_wires]/I");
  fTree->Branch("Wire_max_ledshwtick", wire_max_ledshwtick,  "Wire_max_ledshwtick[n_wires]/I");


  fTree->Branch("Wire_adc",          &wire_adc);
  fTree->Branch("Wire_corradc",      &wire_corradc);
  fTree->Branch("Wire_tick",         &wire_tick);
  fTree->Branch("Wire_shwflag",      &wire_shwflag);
  fTree->Branch("Wire_ledshwflag",   &wire_ledshwflag);
  fTree->Branch("Wire_x",            &wire_x);
  fTree->Branch("Wire_pdg",          &wire_pdg);
  fTree->Branch("Wire_pdg_mom",      &wire_pdg_mom);
  fTree->Branch("Wire_dr",           &wire_dr);
  fTree->Branch("Wire_dist",         &wire_dist);
  fTree->Branch("Wire_tag_pdg",      &wire_tag_pdg);
  fTree->Branch("Wire_tag_pdg_mom",  &wire_tag_pdg_mom);
  fTree->Branch("Wire_tag_is_shw",   &wire_tag_is_shw);
  fTree->Branch("Wire_tag_ith",      &wire_tag_ith);
  fTree->Branch("Wire_tag_true_px",  &wire_tag_true_px);
  fTree->Branch("Wire_tag_true_py",  &wire_tag_true_py);
  fTree->Branch("Wire_tag_true_pz",  &wire_tag_true_pz);
  fTree->Branch("Wire_tag_reco_px",  &wire_tag_reco_px);
  fTree->Branch("Wire_tag_reco_py",  &wire_tag_reco_py);
  fTree->Branch("Wire_tag_reco_pz",  &wire_tag_reco_pz);
  fTree->Branch("Wire_tag_reco_eng",  &wire_tag_reco_eng);


  fTree->Branch("Ch_charge",         &ch_chg);
  fTree->Branch("Ch_tick",           &ch_tick);
  fTree->Branch("Wire_ROIfirst",     wire_roifirst,  "Wire_ROIfirst[n_wires]/I");
  fTree->Branch("Wire_ROIend",       wire_roiend,    "Wire_ROIend[n_wires]/I");
  fTree->Branch("Sum_chchg",         sum_ch_chg,     "Sum_chchg[n_wires]/D");

  fTree->Branch("TrueParticle_size", &TrueParticle_size, "TrueParticle_size/I");
  fTree->Branch("TruePart_PDG",      &TruePart_PDG);
  fTree->Branch("TruePart_E",        &TruePart_E);

  fTree->Branch("trueDaughterPDGs",    &trueDaughterPDGs);
  fTree->Branch("trueDaughterE",       &trueDaughterE);
  fTree->Branch("trueDaughterDepoE",   &trueDaughterDepoE);
  fTree->Branch("trueDaughterPx",      &trueDaughterPx);
  fTree->Branch("trueDaughterPy",      &trueDaughterPy);
  fTree->Branch("trueDaughterPz",      &trueDaughterPz);

  fTree->Branch("ChgU",              &ChgU);
  fTree->Branch("ChgV",              &ChgV);
  fTree->Branch("ChgZ",              &ChgZ);

  fTree->Branch("n_tracks_pm",        &n_tracks_pmtrack, "n_tracks_pm/I");

  fTree->Branch("leadingShowerID",    &leadingShowerID,  "leadingShowerID/I");


  fTree->Branch("n_pand_par",        &n_pand_particles, "n_pand_par/I");
  fTree->Branch("n_prims",           &n_prims, "n_prims/I");
  fTree->Branch("m_pf_vtx",          &m_pf_vtx, "m_pf_vtx/D");
  fTree->Branch("m_pf_vtx_x",        &m_pf_vtx_x, "m_pf_vtx_x/D");
  fTree->Branch("m_pf_vtx_y",        &m_pf_vtx_y, "m_pf_vtx_y/D");
  fTree->Branch("m_pf_vtx_z",        &m_pf_vtx_z, "m_pf_vtx_z/D");

  fTree->Branch("PFVtxTPC",          PFVtxTPC, "PFVtxTPC[3]/D");
  fTree->Branch("PFVtxWire",         PFVtxWire, "PFVtxWire[3]/D");
  fTree->Branch("PFVtxTime",         PFVtxTime, "PFVtxTime[3]/D");
  fTree->Branch("PFVtxGlobalWire",   PFVtxGlobalWire, "PFVtxGlobalWire[3]/D");
  fTree->Branch("PFVtxGlobalTick",   PFVtxGlobalTick, "PFVtxGlobalTick[3]/D");
  fTree->Branch("PFVtxGlobalPlane",   PFVtxGlobalPlane, "PFVtxGlobalPlane[3]/D");

  fTree->Branch("cvnisantineutrino", &fCVNResultIsAntineutrino, "cvnisantineutrino/D");
  fTree->Branch("cvnnue",            &fCVNResultNue,            "cvnnue/D");
  fTree->Branch("cvnnumu",           &fCVNResultNumu,           "cvnnumu/D");
  fTree->Branch("cvnnutau",          &fCVNResultNutau,          "cvnnutau/D");
  fTree->Branch("cvnnc",             &fCVNResultNC,             "cvnnc/D");
  fTree->Branch("cvn0protons",       &fCVNResult0Protons,       "cvn0protons/D");
  fTree->Branch("cvn1protons",       &fCVNResult1Protons,       "cvn1protons/D");
  fTree->Branch("cvn2protons",       &fCVNResult2Protons,       "cvn2protons/D");
  fTree->Branch("cvnNprotons",       &fCVNResultNProtons,       "cvnNprotons/D");
  fTree->Branch("cvn0pions",         &fCVNResult0Pions,         "cvn0pions/D");
  fTree->Branch("cvn1pions",         &fCVNResult1Pions,         "cvn1pions/D");
  fTree->Branch("cvn2pions",         &fCVNResult2Pions,         "cvn2pions/D");
  fTree->Branch("cvnNpions",         &fCVNResultNPions,         "cvnNpions/D");
  fTree->Branch("cvn0pizeros",       &fCVNResult0Pizeros,       "cvn0pizeros/D");
  fTree->Branch("cvn1pizeros",       &fCVNResult1Pizeros,       "cvn1pizeros/D");
  fTree->Branch("cvn2pizeros",       &fCVNResult2Pizeros,       "cvn2pizeros/D");
  fTree->Branch("cvnNpizeros",       &fCVNResultNPizeros,       "cvnNpizeros/D");
  fTree->Branch("cvn0neutrons",      &fCVNResult0Neutrons,      "cvn0neutrons/D");
  fTree->Branch("cvn1neutrons",      &fCVNResult1Neutrons,      "cvn1neutrons/D");
  fTree->Branch("cvn2neutrons",      &fCVNResult2Neutrons,      "cvn2neutrons/D");
  fTree->Branch("cvnNneutrons",      &fCVNResultNNeutrons,      "cvnNneutrons/D");


}

void recoestudy::RecoEnergyS::reconfigure(fhicl::ParameterSet const& pset) {
  fHitsModuleLabel     = pset.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel  = pset.get<std::string>("ClusterModuleLabel");
  fShowerModuleLabel   = pset.get<std::string>("ShowerModuleLabel");
  fTrackModuleLabel1   = pset.get<std::string>("TrackModuleLabel1");
  fTrackModuleLabel2   = pset.get<std::string>("TrackModuleLabel2");
  fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
  fVertexModuleLabel   = pset.get<std::string>("VertexModuleLabel");
  fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
  fPandoraNuVertexModuleLabel = pset.get<std::string>("PandoraNuVertexModuleLabel");
  fMCGenModuleLabel = pset.get<std::string>("MCGenModuleLabel");
  fGlobalWireMethod = pset.get<int>("GlobalWireMethod");
  fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
  fEnergyRecoNumuLabel = pset.get<std::string>("EnergyRecoNumuLabel");
  fEnergyRecoNuelseLabel = pset.get<std::string>("EnergyRecoNuelseLabel");
  /*fRegCNNModuleLabel = pset.get<std::string>("RegCNNModuleLabel");
  fRegCNNResultLabel = pset.get<std::string>("RegCNNResultLabel");
  fRegCNN1ststepModuleLabel = "regcnneval";
  fRegCNN1ststepResultLabel = "regcnnresult";
  fRegCNN2ndstepModuleLabel = "regcnnevalk";
  fRegCNN2ndstepResultLabel = "regcnnresultk";
  fRegCNNModule2Label = pset.get<std::string>("RegCNNModule2Label");
  fRegCNNResult2Label = pset.get<std::string>("RegCNNResult2Label");*/
  fIsVD = pset.get<bool>("IsVD");

  ievt = -999;
  subrun = -999;
  run = -999;

  v3_on_planes[0] = new TVector3();
  v3_on_planes[1] = new TVector3();
  v3_on_planes[2] = new TVector3();

}

std::vector< int> getUniques(std::vector< int> coll)
{
  std::vector< int> uniques;
  for ( int tpc : coll)
  {
    if (std::find(uniques.begin(), uniques.end(), tpc) == uniques.end())
        uniques.push_back(tpc);
  }
 return uniques;
}

std::vector<std::pair<const simb::MCParticle*,double>> get_sortedMCParticle(std::unordered_map<const simb::MCParticle*,double> mcEMap){
     std::vector<std::pair<const simb::MCParticle*, double>> outVec;
     double total_E = 0;
     for (const std::pair<const simb::MCParticle*,double> p : mcEMap){
  	   outVec.push_back(p);
	   total_E += p.second;
     }
     std::sort(outVec.begin(), outVec.end(),
  	   [](std::pair<const simb::MCParticle*,double>a, std::pair<const simb::MCParticle*,double>b)
  	   {return a.second > b.second;});
     if (abs(total_E) < 1e-6) { total_E = 1; } // Protect against zero division
     for (std::pair<const simb::MCParticle*,double> & p : outVec) {
	p.second /= total_E;
     }

     //std::cout << "size of MCParticle associated with reco shower " << outVec.size() << std::endl;
     return outVec;

}
const simb::MCParticle* recoestudy::RecoEnergyS::get_bestMCParticle(const art::Ptr<recob::Hit>& hit, detinfo::DetectorClocksData const& clockData){
    std::unordered_map<const simb::MCParticle*, double> mcEMap;
    std::vector<sim::TrackIDE> TrackIDs = backtracker->HitToTrackIDEs(clockData, hit);
    for(size_t e = 0; e < TrackIDs.size(); ++e){
        mcEMap[particleinventory->TrackIdToParticle_P(TrackIDs[e].trackID)] += TrackIDs[e].energy;
    }
    auto outVec = get_sortedMCParticle(mcEMap);
    if (outVec.size()>0) return outVec[0].first;
    else return NULL;

}

bool recoestudy::RecoEnergyS::ConvertXYZtoWireTDC(const double *vtx_loc,
  	 double *TPC, double *LocalWire, double *LocalTick,
	 double *GlobalPlane, double *GlobalWire, double *GlobalTick, detinfo::DetectorPropertiesData const& detProp){

  geo::Point_t vtx(vtx_loc[0], vtx_loc[1], vtx_loc[2]);
  if (!(geom->FindTPCAtPosition(vtx).isValid)) return false;
  else {
    int rawcrys = 0;
    for (int iplane = 0; iplane<3; iplane++){
      int rawtpc = (int) (geom->FindTPCAtPosition(vtx)).TPC;
      geo::PlaneID planeID{(unsigned int) rawcrys, (unsigned int) rawtpc, (unsigned int) iplane};
	    geo::PlaneGeo const& planegeo_temp = wireReadout->Plane(planeID);
	    geo::WireID w1;
	    try {
  	          w1 = wireReadout->Plane(planeID).NearestWireID(vtx);
	    }
	    catch (geo::InvalidWireError const& e) {
	          if (!e.hasSuggestedWire()) throw;
	          w1 = planegeo_temp.ClosestWireID(e.suggestedWireID());
	    }
	    //double time1 = detprop->ConvertXToTicks(regvtx_loc[0]+trueParticle->T()*detprop->DriftVelocity()*1e-3, iplane, rawtpc, rawcrys);
	    double time1 = detProp.ConvertXToTicks(vtx_loc[0], iplane, rawtpc, rawcrys);
	    TPC[iplane]  = (double)rawtpc;
            LocalWire[iplane] = (double)w1.Wire;
	    LocalTick[iplane] = time1;
	    if (fGlobalWireMethod == 1){
		    GlobalPlane[iplane] = iplane;
		    GlobalWire[iplane]  = GetPixelMapWireIdx(w1);
    	    	    GlobalTick[iplane]  = (rawtpc%2==0) ? -time1 : time1;
	    	    if(rawtpc%2==1) GlobalWire[iplane] += hit_offset[iplane];
	    } else if (fGlobalWireMethod == 2){
		    unsigned int globalWire  = w1.Wire;
		    unsigned int globalPlane = w1.Plane;
		    double globalTDC = time1;
		    GetDUNEGlobalWireTDC(w1, time1, globalWire, globalPlane, globalTDC, detProp.DriftVelocity());
		    GlobalPlane[iplane]= (double) globalPlane;
		    GlobalWire[iplane] = (double) globalWire;
		    GlobalTick[iplane] = globalTDC;
	    } else {
		   std::cout << "Wrong Global Wire Method" << std::endl;
	    }
    }
    return true;
  } // end of inTPC

} // end of ConvertXYZtoWireTDC

bool recoestudy::RecoEnergyS::passesNHitsCut(const simb::MCParticle *part) {
    int nHits = backtracker->TrackIdToSimIDEs_Ps(part->TrackId()).size();
    if (part->PdgCode() == 2112) {
        for(int d = 0; d < part->NumberDaughters(); ++d){
          nHits += backtracker->TrackIdToSimIDEs_Ps(part->Daughter(d)).size();
        }
    }
    std::cout << "Particle with PDG " << part->PdgCode() << " has " << nHits << " hits" << std::endl;
    return (nHits >= daughterHitCut);
}

void recoestudy::RecoEnergyS::addDaughters(const simb::MCParticle *part, const std::vector<int> finalStatePdgs, std::vector< const simb::MCParticle * > &allParts) {
    int parentTrackID = part->TrackId();
    std::cout << "Looking for daughters..." << std::endl;
    for (const simb::MCParticle * potDaughterPart: allParts) {
        if (parentTrackID != potDaughterPart->Mother()) continue;
        int currentPdg = potDaughterPart->PdgCode();
        if (currentPdg > 1000000) continue; // get rid of nuclear products
        if (abs(currentPdg) == 12 || abs(currentPdg) == 14 || abs(currentPdg) == 16) continue;
        std::cout << "Found daughter " << currentPdg << std::endl;
        if (std::find(finalStatePdgs.begin(), finalStatePdgs.end(), currentPdg) != finalStatePdgs.end()) { // It's in the finalStatePdgs list
            if (!passesNHitsCut(potDaughterPart)) continue;
            trueDaughterPDGs.push_back(currentPdg);
            trueDaughterE.push_back(potDaughterPart->E());
            trueDaughterPx.push_back(potDaughterPart->Px());
            trueDaughterPy.push_back(potDaughterPart->Py());
            trueDaughterPz.push_back(potDaughterPart->Pz());
            std::vector< const sim::IDE * > ides = backtracker->TrackIdToSimIDEs_Ps( potDaughterPart->TrackId(), geo::kZ );
            float Edep = 0;
            for (auto ide: ides) {
              Edep += ide->energy / 1000;
            }
            trueDaughterDepoE.push_back( Edep );
        }
        else {
            addDaughters(potDaughterPart, finalStatePdgs, allParts);
        }
    }
    std::cout << "Done looking for daughters" << std::endl;
}


void recoestudy::RecoEnergyS::analyze(art::Event const& evt) {
  //std::cout << "Start Ana" << std::endl;

  /// Analyse function to save information for calibrating shower energies
  /// This is written assuming single particle per event
std::cout<<"recoenergy RecoEnergyS::analyze(art::Event const& evt):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-):-)"<<std::endl;
std::cout<<"fGlobalWireMethod = "<<fGlobalWireMethod<<std::endl;
  this->reset();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detPropEvt = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  //auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  det_t0        = clockData.TriggerOffsetTPC();
  driftvelocity = detPropEvt.DriftVelocity(detPropEvt.Efield(),detPropEvt.Temperature())/1000.;
  SampleRate    = sampling_rate(clockData);
  bool isMC = !evt.isRealData();

  ievt = evt.id().event();
  run  = evt.id().run();
  subrun  = evt.id().subRun();
  //                     planeID, tpc#, CryostatID
  nwires_p0 = wireReadout->Nwires(geo::PlaneID(0,0,0));
  nwires_p1 = wireReadout->Nwires(geo::PlaneID(0,0,1));
  nwires_p2 = wireReadout->Nwires(geo::PlaneID(0,0,2));

  art::InputTag itag1("cvneva", "cvnresult");
  art::Handle<std::vector<cvn::Result>> cvnin;
  evt.getByLabel(itag1, cvnin);
  if( !cvnin.failedToGet() ) {
    //using i = cvn::Interaction;
    //if(cvnin->empty() || (*cvnin)[0].fOutput.size() <= i::kNutauOther){
    //if(cvnin->empty()){
      fCVNResultIsAntineutrino = fCVNResultNue = fCVNResultNumu = fCVNResultNutau = fCVNResultNC = \
      fCVNResult0Protons = fCVNResult1Protons = fCVNResult2Protons = fCVNResultNProtons = \
      fCVNResult0Pions = fCVNResult1Pions = fCVNResult2Pions = fCVNResultNPions = \
      fCVNResult0Pizeros = fCVNResult1Pizeros = fCVNResult2Pizeros = fCVNResultNPizeros = \
      fCVNResult0Neutrons = fCVNResult1Neutrons = fCVNResult2Neutrons = fCVNResultNNeutrons = -3;
    //}
    if(!cvnin->empty()){
      //const std::vector<float>& v = (*cvnin)[0].fOutput;
      //fCVNResultNue = v[i::kNueQE] + v[i::kNueRes] + v[i::kNueDIS] + v[i::kNueOther];
      //fCVNResultNumu = v[i::kNumuQE] + v[i::kNumuRes] + v[i::kNumuDIS] + v[i::kNumuOther];
      //fCVNResultNutau = v[i::kNutauQE] + v[i::kNutauRes] + v[i::kNutauDIS] + v[i::kNutauOther]
  
      //fCVNResultIsAntineutrino = (*cvnin)[0].GetIsAntineutrinoProbability();
      fCVNResultNC = (*cvnin)[0].fOutput[0][0];
      fCVNResultNue = (*cvnin)[0].fOutput[0][1];
      fCVNResultNumu = (*cvnin)[0].fOutput[0][2];

/*
      fCVNResultNue = (*cvnin)[0].GetNueProbability();
      fCVNResultNumu = (*cvnin)[0].GetNumuProbability();
      fCVNResultNutau = (*cvnin)[0].GetNutauProbability();
      fCVNResultNC = (*cvnin)[0].GetNCProbability();
  
      fCVNResult0Protons = (*cvnin)[0].Get0protonsProbability();
      fCVNResult1Protons = (*cvnin)[0].Get1protonsProbability();
      fCVNResult2Protons = (*cvnin)[0].Get2protonsProbability();
      fCVNResultNProtons = (*cvnin)[0].GetNprotonsProbability();
  
      fCVNResult0Pions = (*cvnin)[0].Get0pionsProbability();
      fCVNResult1Pions = (*cvnin)[0].Get1pionsProbability();
      fCVNResult2Pions = (*cvnin)[0].Get2pionsProbability();
      fCVNResultNPions = (*cvnin)[0].GetNpionsProbability();
  
      fCVNResult0Pizeros = (*cvnin)[0].Get0pizerosProbability();
      fCVNResult1Pizeros = (*cvnin)[0].Get1pizerosProbability();
      fCVNResult2Pizeros = (*cvnin)[0].Get2pizerosProbability();
      fCVNResultNPizeros = (*cvnin)[0].GetNpizerosProbability();
  
      fCVNResult0Neutrons = (*cvnin)[0].Get0neutronsProbability();
      fCVNResult1Neutrons = (*cvnin)[0].Get1neutronsProbability();
      fCVNResult2Neutrons = (*cvnin)[0].Get2neutronsProbability();
      fCVNResultNNeutrons = (*cvnin)[0].GetNneutronsProbability();
*/
    }
  }

  // Get the hits out of the event record
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Get the clusters out of the event record
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > clusters;
  if (evt.getByLabel(fClusterModuleLabel,clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  art::FindManyP<recob::Cluster> fmc(hitHandle, evt, fClusterModuleLabel);

  // Get the emshower out of the event record
  art::Handle<std::vector<recob::Shower> > showerHandle;
  std::vector<art::Ptr<recob::Shower> > showers;
  if (evt.getByLabel(fShowerModuleLabel,showerHandle))
    art::fill_ptr_vector(showers, showerHandle);

  leadingShowerID = -999;
  long unsigned int nLeadingShowerHits = 0;
  for (art::Ptr<recob::Shower> show: showers) {
      if (show->dEdx().size() > nLeadingShowerHits) {
          nLeadingShowerHits = show->dEdx().size();
          leadingShowerID = show->ID();
      }
  }

  art::FindManyP<recob::Cluster> fmshc(showerHandle, evt, fClusterModuleLabel);

  // Get the pmtrack out of the event record
  art::Handle<std::vector<recob::Track> > trackHandle1;
  std::vector<art::Ptr<recob::Track> > tracks_pmtrack;
  if (evt.getByLabel(fTrackModuleLabel1,trackHandle1))
    art::fill_ptr_vector(tracks_pmtrack, trackHandle1);

  // Get the pandoratrack out of the event record
  art::Handle<std::vector<recob::Track> > trackHandle2;
  std::vector<art::Ptr<recob::Track> > tracks_pand;
  if (evt.getByLabel(fTrackModuleLabel2,trackHandle2))
    art::fill_ptr_vector(tracks_pand, trackHandle2);

  // Get the vertex out of the event record
  art::Handle<std::vector<recob::Vertex> > vertexHandle;
  std::vector<art::Ptr<recob::Vertex> > vertex_list;
  if (evt.getByLabel(fVertexModuleLabel,vertexHandle))
    art::fill_ptr_vector(vertex_list, vertexHandle);

  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (isMC){
    if (evt.getByLabel(fMCGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
  }

  // Get RegCNN Results
  /*
  art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle;
  evt.getByLabel(fRegCNNModuleLabel, fRegCNNResultLabel, cnnresultListHandle);

  art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle1ststep;
  evt.getByLabel(fRegCNN1ststepModuleLabel, fRegCNN1ststepResultLabel, cnnresultListHandle1ststep);

  art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle2ndstep;
  evt.getByLabel(fRegCNN2ndstepModuleLabel, fRegCNN2ndstepResultLabel, cnnresultListHandle2ndstep);

  art::Handle<std::vector<cnn::RegCNNResult>> cnnresultListHandle2;
  evt.getByLabel(fRegCNNModule2Label, fRegCNNResult2Label, cnnresultListHandle2);
  */


  art::FindManyP<recob::Wire> fmwire(hitHandle, evt, fHitsModuleLabel);
  art::FindManyP<raw::RawDigit> fmraw(hitHandle, evt, fHitsModuleLabel);
  art::FindMany<recob::PFParticle> fmPFParticle(vertexHandle, evt, fPFParticleModuleLabel);

  art::FindManyP<recob::SpacePoint> spFromHits(hitHandle,evt,fPFParticleModuleLabel);


  art::FindManyP<recob::Hit> fmsh(showerHandle, evt, fShowerModuleLabel);
  art::FindManyP<recob::Hit> fmtrk(trackHandle2, evt, fTrackModuleLabel2);

  double fIntShwEnergy = -0.02;
  double fGradShwEnergy = 0.985;


  //std::cout << "Loop simchannels" << std::endl;
  // Find the energy deposited on each plane in the TPC
  int icount = 0;
  //const std::vector<const sim::SimChannel*>& simChannels = backtracker->SimChannels();
  const std::vector<art::Ptr< sim::SimChannel > >& simChannels = backtracker->SimChannels();
  simCh_size = simChannels.size();
  for (auto channelIt = simChannels.begin(); channelIt != simChannels.end(); ++channelIt) {
    if (icount >= kMaxHits) continue;
    int plane = wireReadout->View((*channelIt)->Channel());
    tdc_size[icount] = ((*channelIt)->TDCIDEMap()).size();
    for (auto const& tdcIt : (*channelIt)->TDCIDEMap()) {
      for (auto const& ideIt : tdcIt.second) {
        switch (plane) {
          case geo::kU:
            depositU += ideIt.energy;
            ChgU.push_back(ideIt.energy);
            break;
          case geo::kV:
            depositV += ideIt.energy;
            ChgV.push_back(ideIt.energy);
            break;
          case geo::kZ:
            depositZ += ideIt.energy;
            ChgZ.push_back(ideIt.energy);
            break;
        }
      }
    }
    icount ++;
  }

 //-----------------------------------------------------------------------------------------
  //std::cout << "True Info" << std::endl;
  // use simChannel to extract energy and charge per TDC tick
  // true info


  nu_truth_N = 0;
  if (mclist.size()>0){
        int neutrino_i = 0;
        for(unsigned int iList = 0; (iList < mclist.size()) && (neutrino_i < kMax) ; ++iList){
          if (mclist[iList]->NeutrinoSet()){
            //nueng_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().E();
            nueng_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().Momentum().E();
            nuvtxx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vx();
            nuvtxy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vy();
            nuvtxz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vz();
	    nupdg_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().PdgCode();
            nuccnc_truth[neutrino_i] = mclist[iList]->GetNeutrino().CCNC();
            numode_truth[neutrino_i] = mclist[iList]->GetNeutrino().Mode();

            leppdg_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().PdgCode();
            lepp0_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Lepton().Momentum().X();
            lepp1_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Lepton().Momentum().Y();
            lepp2_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Lepton().Momentum().Z();
	    lepeng_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Momentum().T();

	    neutrino_i++;
          }
	}
	nu_truth_N = neutrino_i;
  }
  // skip NC events
  //if (nuccnc_truth[0] == 1) return;

   // Get DUNE energy Reco
  art::Handle<dune::EnergyRecoOutput> engrecoHandle;
  if (fabs(nupdg_truth[0])==12){
    evt.getByLabel(fEnergyRecoNueLabel,engrecoHandle);
  } else if (fabs(nupdg_truth[0])==14){
    evt.getByLabel(fEnergyRecoNumuLabel,engrecoHandle);
  } else {
    evt.getByLabel(fEnergyRecoNuelseLabel,engrecoHandle);
  }
  trkf::TrackMomentumCalculator trkm{};

  //true CC event with true vertex in fiducial volume

  nutrue_fid = 0;
  if (nu_truth_N>0){
    if(fabs(nuvtxx_truth[0]) < 310. && fabs(nuvtxy_truth[0]) < 550. && nuvtxz_truth[0] > 50. && nuvtxz_truth[0] < 1250.) nutrue_fid = 1;
    if (fabs(nupdg_truth[0])==14){
        // 1 logest reco track contained
        // 0 logest reco track exiting
        mutrk_contain = engrecoHandle->longestTrackContained;
    }

  }

  //const sim::ParticleList& trueParticles = backtracker->ParticleList();
  const sim::ParticleList& trueParticles = particleinventory->ParticleList();
  const simb::MCParticle* trueParticle = trueParticles.Primary(0);
//  const simb::MCParticle* trueParticle = particleIt->second;

  primpdg    = trueParticle->PdgCode();
  trueEnergy = trueParticle->Momentum().E();
  trueVertex_X = trueParticle->Vx();
  trueVertex_Y = trueParticle->Vy();
  trueVertex_Z = trueParticle->Vz();
  trueEnd_X = trueParticle->EndX();
  trueEnd_Y = trueParticle->EndY();
  trueEnd_Z = trueParticle->EndZ();
  std::cout << "TrackID 0? " << primpdg << " " << trueEnergy << std::endl;
  std::cout << "   " << nu_truth_N << " " << nupdg_truth[0] << " " << nueng_truth[0] << std::endl;

  truePx       = trueParticle->Px();
  truePy       = trueParticle->Py();
  truePz       = trueParticle->Pz();

  const std::vector<int> finalStatePdgs = {11, -11, 13, -13, 22, 211, -211, 321, -321, 2212, -2212, 2112, -2112};
  if (mclist.size() > 0) {
      std::vector< const simb::MCParticle * > parts = particleinventory->MCTruthToParticles_Ps(mclist[0]);
      std::cout << "Looking at primary daughters:" << std::endl;
      for (long unsigned int i = 0; i < parts.size(); i++) {
          const simb::MCParticle *part = parts.at(i);
          int currentPdg = part->PdgCode();
          if (part->Mother() != 0) continue;
          if (part->StatusCode() != 1) continue;
          if (currentPdg > 1000000) continue; // get rid of nuclear products
          if (abs(currentPdg) == 12 || abs(currentPdg) == 14 || abs(currentPdg) == 16) continue;
          std::cout << currentPdg << std::endl;
          if (std::find(finalStatePdgs.begin(), finalStatePdgs.end(), currentPdg) != finalStatePdgs.end()) { // It's in the finalStatePdgs list
              if (!passesNHitsCut(part)) continue;
              trueDaughterPDGs.push_back(currentPdg);
              trueDaughterE.push_back(part->E());
              trueDaughterPx.push_back(part->Px());
              trueDaughterPy.push_back(part->Py());
              trueDaughterPz.push_back(part->Pz());
              std::vector< const sim::IDE * > ides = backtracker->TrackIdToSimIDEs_Ps( part->TrackId(), geo::kZ );
              float Edep = 0;
              for (auto ide: ides) {
                Edep += ide->energy / 1000;
              }
              trueDaughterDepoE.push_back( Edep );
          }
          else {
              addDaughters(part, finalStatePdgs, parts);
          }
      }
  }



  // Find the distance between the particle vertex and the edge of the detector
  TVector3 vertex = TVector3(trueParticle->Vx(),trueParticle->Vy(),trueParticle->Vz());
  TVector3 end = TVector3(trueParticle->EndX(),trueParticle->EndY(),trueParticle->EndZ());
  TVector3 direction = TVector3(trueParticle->Px(),trueParticle->Py(),trueParticle->Pz()).Unit();
  bool inTPC = false;
  TVector3 pos;

  double const vertex_loc[3] = {trueVertex_X, trueVertex_Y, trueVertex_Z};
  double const end_loc[3]    = {trueEnd_X, trueEnd_Y, trueEnd_Z};
  vtx_in_tpc = 0;
  double temp_tpc[3];
  double temp_plane[3];

  inTPC = ConvertXYZtoWireTDC(vertex_loc, trueVtxTPC, trueVtxWire, trueVtxTime,
	     trueVtxGlobalPlane, trueVtxGlobalWire, trueVtxGlobalTick, detPropEvt);
  ConvertXYZtoWireTDC(end_loc, temp_tpc, trueEndWire, trueEndTime,
         temp_plane, trueEndGlobalWire, trueEndGlobalTime, detPropEvt);
  if (inTPC) vtx_in_tpc = 1;

  /*
  // Get RegCNN Results
  if (!cnnresultListHandle.failedToGet())
  {
      if (!cnnresultListHandle->empty())
      {
          const std::vector<float>& v = (*cnnresultListHandle)[0].fOutput;
          for (unsigned int ii = 0; ii < 3; ii++){
              CVN_vertex[ii] = v[ii];
          }
      }
  }

  if (!cnnresultListHandle1ststep.failedToGet())
  {
      if (!cnnresultListHandle1ststep->empty())
      {
          const std::vector<float>& v = (*cnnresultListHandle1ststep)[0].fOutput;
          for (unsigned int ii = 0; ii < 3; ii++){
              CVN_vertex1ststep[ii] = v[ii];
          }
      }
  }
  if (!cnnresultListHandle2ndstep.failedToGet())
  {
      if (!cnnresultListHandle2ndstep->empty())
      {
          const std::vector<float>& v = (*cnnresultListHandle2ndstep)[0].fOutput;
          for (unsigned int ii = 0; ii < 3; ii++){
              CVN_vertex2ndstep[ii] = v[ii];
          }
      }
  }

  if (!cnnresultListHandle2.failedToGet())
  {
      if (!cnnresultListHandle2->empty())
      {
          const std::vector<float>& v = (*cnnresultListHandle2)[0].fOutput;
          for (unsigned int ii = 0; ii < 3; ii++){
              CVN_vertex2[ii] = v[ii];
          }
      }
  }
  */


  //-----------------------------------------------------------------------------------------
  /*
  // find the wire and tick of reconstructed vertices
  double const regvtx_loc[3] = {(double)CVN_vertex[0], (double)CVN_vertex[1], (double)CVN_vertex[2]};
  inTPC = ConvertXYZtoWireTDC(regvtx_loc, CVNVtxTPC, CVNVtxWire, CVNVtxTime,
		  CVNVtxGlobalPlane, CVNVtxGlobalWire, CVNVtxGlobalTick, detPropEvt);
  if (inTPC){
  }

  double const regvtx_loc2[3] = {(double)CVN_vertex2[0], (double)CVN_vertex2[1], (double)CVN_vertex2[2]};
  inTPC = ConvertXYZtoWireTDC(regvtx_loc2, CVNVtx2TPC, CVNVtx2Wire, CVNVtx2Time,
		  CVNVtx2GlobalPlane, CVNVtx2GlobalWire, CVNVtx2GlobalTick, detPropEvt);
  */

  //-----------------------------------------------------------------------------------------
  // for track info
  n_tracks_pmtrack = 0;
  n_tracks_pand = 0;
  if (trackHandle1.isValid()){
    n_tracks_pmtrack = tracks_pmtrack.size();
  }
  if (trackHandle2.isValid()){
    n_tracks_pand = tracks_pand.size();
  }

  //-----------------------------------------------------------------------------------------
  // Pandora Nu Vertex
  //std::cout << "Pandora Vertex Info" << std::endl;
  lar_pandora::PFParticleVector particleVector;
  lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
  lar_pandora::VertexVector vertexVector;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);

  n_prims = 0;
  n_pand_particles = particleVector.size();
  m_pf_vtx = 0, m_pf_vtx_x = -10000, m_pf_vtx_y = -10000, m_pf_vtx_z = -10000;

  double xyz_temp[3] = {0.0, 0.0, 0.0} ;
  for (unsigned int ipfp = 0; ipfp < particleVector.size(); ipfp++){
    const art::Ptr<recob::PFParticle> particle = particleVector.at(ipfp);
    if (!particle->IsPrimary()) continue;
    n_prims ++;

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

	   m_pf_vtx = 1;
           m_pf_vtx_x = xyz_temp[0];
           m_pf_vtx_y = xyz_temp[1];
           m_pf_vtx_z = xyz_temp[2];
       }
   } // end of if
  } // end of loop particleVector
  inTPC = ConvertXYZtoWireTDC(xyz_temp, PFVtxTPC, PFVtxWire, PFVtxTime,
	     PFVtxGlobalPlane, PFVtxGlobalWire, PFVtxGlobalTick, detPropEvt);
  //-----------------------------------------------------------------------------------------
  // loop over for reconstructed track

  if (trackHandle2.isValid()){

      for (int itrack = 0; itrack < n_tracks_pand; itrack++){
	art::Ptr<recob::Track> ptrack = tracks_pand.at(itrack);
	all_track_length[itrack] = ptrack->Length();
	all_track_startx[itrack] = ptrack->Start().X();
	all_track_starty[itrack] = ptrack->Start().Y();
	all_track_startz[itrack] = ptrack->Start().Z();
    all_track_endx[itrack] = ptrack->End().X();
    all_track_endy[itrack] = ptrack->End().Y();
    all_track_endz[itrack] = ptrack->End().Z();
	all_track_px[itrack] = ptrack->VertexDirection().X();
	all_track_py[itrack] = ptrack->VertexDirection().Y();
	all_track_pz[itrack] = ptrack->VertexDirection().Z();
	all_track_mom[itrack] = ptrack->StartMomentum();
	all_track_calmom[itrack] = (double)trkm.GetTrackMomentum(ptrack->Length(),13);
	TVector3 vv1(ptrack->VertexDirection().X(), ptrack->VertexDirection().Y(), ptrack->VertexDirection().Z());

	std::cout << "Track..." << std::endl;
	vv1.Print();
	std::cout << trkm.GetTrackMomentum(ptrack->Length(),13) << std::endl;

	double dist = pow(all_track_startx[itrack]-m_pf_vtx_x, 2)+pow(all_track_starty[itrack]-m_pf_vtx_y,2)+pow(all_track_startz[itrack]-m_pf_vtx_z,2);
	dist = sqrt(dist);
	double dist_true = pow(all_track_startx[itrack]-trueVertex_X, 2)+pow(all_track_starty[itrack]-trueVertex_Y,2)+pow(all_track_startz[itrack]-trueVertex_Z,2);
	dist_true = sqrt(dist_true);

        if (fmtrk.isValid()){
           std::vector< art::Ptr<recob::Hit> > vhit = fmtrk.at(itrack);
	   std::unordered_map<const simb::MCParticle*, double> mcEMap;

           for (size_t h = 0; h < vhit.size(); ++h){
    	       std::vector<sim::TrackIDE> TrackIDs = backtracker->HitToTrackIDEs(clockData, vhit[h]);
               for(size_t e = 0; e < TrackIDs.size(); ++e){
		 mcEMap[particleinventory->TrackIdToParticle_P(TrackIDs[e].trackID)] += TrackIDs[e].energy;
               }
	   }


	   std::vector<std::pair<const simb::MCParticle*,double>> trktrue = get_sortedMCParticle(mcEMap);
           all_track_dist_true_vtx[itrack] = dist_true;
           all_track_dist_reco_vtx[itrack] = dist;

	   if (trktrue.size()>0){
	        int pdg = trktrue[0].first->PdgCode();
		int trkid = trktrue[0].first->TrackId();
		int mom_trkid = -1;
                int mom_pdg = trktrue[0].first->Mother()==0 ? -1 : trueParticles[trktrue[0].first->Mother()]->PdgCode();
		for (auto itPart: trueParticles){
			const simb::MCParticle* pPart = itPart.second;
			if (pPart->Process() != "primary") continue;
			int iid = pPart->TrackId();
			if (trkid==iid) continue;
			for (int imm = 0; imm < pPart->NumberDaughters(); imm++){
				int jid = pPart->Daughter(imm);
				if (jid == trkid){
					std::cout << "Track ---> " <<  imm << " " << trkid << " " << iid << " : " << pdg << " " << pPart->PdgCode() << " " << mom_pdg << " " << pPart->PdgCode() << std::endl;
					mom_trkid = iid;
					mom_pdg = pPart->PdgCode();
					break;
				}
			}
		}
	        int pdg2 = -1;
		int trkid2 = -1;
		int mom_trkid2 = -1;
		int mom_pdg2 = -1;
		if (trktrue.size()>1){
			pdg2 = trktrue[1].first->PdgCode();
			trkid2 = trktrue[1].first->TrackId();
                	mom_pdg2 = trktrue[1].first->Mother()==0 ? -1 : trueParticles[trktrue[1].first->Mother()]->PdgCode();
			for (auto itPart: trueParticles){
				const simb::MCParticle* pPart = itPart.second;
				if (pPart->Process() != "primary") continue;
				int iid = pPart->TrackId();
				if (trkid2==iid) continue;
				for (int imm = 0; imm < pPart->NumberDaughters(); imm++){
					int jid = pPart->Daughter(imm);
					if (jid == trkid2){
						mom_trkid2 = iid;
						mom_pdg2 = pPart->PdgCode();
						break;
					}
				}
			}
		} // end of trktrue.size().1


	        if (true){
	            TVector3 v3_true(trktrue[0].first->Momentum().Vect());
	            float dang_true = (float)(vv1.Angle(v3_true));
	            all_track_true_Efrac[itrack] = trktrue[0].second;
	            all_track_true_eng[itrack] = trktrue[0].first->Momentum().E();
	            all_track_true_pdg[itrack] = pdg;
	            all_track_true_pdg_mom[itrack] = mom_pdg;
	            all_track_min_angle[itrack] = dang_true;
	            all_track_true_px[itrack] = v3_true.X();
	            all_track_true_py[itrack] = v3_true.Y();
	            all_track_true_pz[itrack] = v3_true.Z();
                all_track_true_endx[itrack] = trktrue[0].first->EndX();
                all_track_true_endy[itrack] = trktrue[0].first->EndY();
                all_track_true_endz[itrack] = trktrue[0].first->EndZ();
		    if (mom_trkid>0){
	  	          all_track_true_mom_px[itrack] = trueParticles[mom_trkid]->Momentum().Vect().X();
	  	          all_track_true_mom_py[itrack] = trueParticles[mom_trkid]->Momentum().Vect().Y();
	  	          all_track_true_mom_pz[itrack] = trueParticles[mom_trkid]->Momentum().Vect().Z();
		    }
		    if (trktrue.size()>1){
	                TVector3 v3_true(trktrue[1].first->Momentum().Vect());
	                all_track2_true_Efrac[itrack] = trktrue[1].second;
	                all_track2_true_eng[itrack] = trktrue[1].first->Momentum().E();
	                all_track2_true_pdg[itrack] = pdg2;
	                all_track2_true_pdg_mom[itrack] = mom_pdg2;
	                all_track2_true_px[itrack] = v3_true.X();
	                all_track2_true_py[itrack] = v3_true.Y();
	                all_track2_true_pz[itrack] = v3_true.Z();
		        if (mom_trkid2>0){
	  	              all_track2_true_mom_px[itrack] = trueParticles[mom_trkid2]->Momentum().Vect().X();
	  	              all_track2_true_mom_py[itrack] = trueParticles[mom_trkid2]->Momentum().Vect().Y();
	  	              all_track2_true_mom_pz[itrack] = trueParticles[mom_trkid2]->Momentum().Vect().Z();
		        }
		    } // end of trktrue.size()>1

		}
	    }
	} // end of fmtrk
      } // end of itrack
  }
  //-----------------------------------------------------------------------------------------
  // loop over for reconstructed shower
  std::map<int, double> trkide;
  float showerCharge = 0;
  int icount_shower = 0;
  hit_shw_count = 0;
  shw_correctedChargeU = 0;
  shw_correctedChargeV = 0;
  shw_correctedChargeZ = 0;
  shw_ChargeU = 0;
  shw_ChargeV = 0;
  shw_ChargeZ = 0;
  n_showers = 0;
  ntpcs = 0;
  int tmp_tpc = -10;

  //std::cout << "EM Shower" << std::endl;
  double max_shw_energy = -1000;
  int electron_candidate = -1;
  int e_candidate_id = -1;
  std::cout << electron_candidate << std::endl;
  if (showerHandle.isValid()){
     n_showers = showers.size();
     TVector3 v3_lep(lepp0_truth[0], lepp1_truth[0], lepp2_truth[0]);

     for (int ishower = 0; ishower < n_showers; ++ishower){
	if (ishower >= kMax) continue;
         art::Ptr<recob::Shower> emshower = showers.at(ishower);

 	 all_shw_idx[ishower] = emshower->ID();
	 all_shw_eng[ishower] = emshower->Energy()[2];
	 all_shw_dedx[ishower] = emshower->dEdx()[2];
	 all_shw_length[ishower] = emshower->Length();
         all_shw_startx[ishower] = emshower->ShowerStart().X();
         all_shw_starty[ishower] = emshower->ShowerStart().Y();
         all_shw_startz[ishower] = emshower->ShowerStart().Z();
         all_shw_px[ishower] = emshower->Direction().X();
         all_shw_py[ishower] = emshower->Direction().Y();
         all_shw_pz[ishower] = emshower->Direction().Z();

	 TVector3 vv1(emshower->Direction().X(), emshower->Direction().Y(),emshower->Direction().Z());
         all_shw_costh[ishower]  = cos(vv1.Theta());
         all_shw_phi[ishower]    = vv1.Phi();

	 double dist = pow(all_shw_startx[ishower]-m_pf_vtx_x, 2)+pow(all_shw_starty[ishower]-m_pf_vtx_y,2)+pow(all_shw_startz[ishower]-m_pf_vtx_z,2);
	 dist = sqrt(dist);
	 double dist_true = pow(all_shw_startx[ishower]-trueVertex_X, 2)+pow(all_shw_starty[ishower]-trueVertex_Y,2)+pow(all_shw_startz[ishower]-trueVertex_Z,2);
	 dist_true = sqrt(dist_true);


         if (fmsh.isValid()){
           std::vector< art::Ptr<recob::Hit> > vhit = fmsh.at(ishower);
	   std::unordered_map<const simb::MCParticle*, double> mcEMap;

	   std::cout << " N Hits in Shower: " << vhit.size() << std::endl;
           for (size_t h = 0; h < vhit.size(); ++h){
    	       std::vector<sim::TrackIDE> TrackIDs = backtracker->HitToTrackIDEs(clockData, vhit[h]);
               for(size_t e = 0; e < TrackIDs.size(); ++e){
                 trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
		 mcEMap[particleinventory->TrackIdToParticle_P(TrackIDs[e].trackID)] += TrackIDs[e].energy;
               }
	   }
           // get best MCParticle candidate based on deposited energy
	   std::vector<std::pair<const simb::MCParticle*,double>> shwtrue = get_sortedMCParticle(mcEMap);
           all_shw_dist_true_vtx[ishower] = dist_true;
           all_shw_dist_reco_vtx[ishower] = dist;
	   if (shwtrue.size()>0){
	        int pdg = shwtrue[0].first->PdgCode();
		int trkid = shwtrue[0].first->TrackId();
		int mom_trkid = -1;
	        int mom_pdg = shwtrue[0].first->Mother()==0 ? -1 : trueParticles[shwtrue[0].first->Mother()]->PdgCode();
		for (auto itPart: trueParticles){
			const simb::MCParticle* pPart = itPart.second;
			if (pPart->Process() != "primary") continue;
			int iid = pPart->TrackId();
			if (trkid==iid) continue;
			for (int imm = 0; imm < pPart->NumberDaughters(); imm++){
				int jid = pPart->Daughter(imm);
				if (jid == trkid){
					std::cout << "Shower ---> " <<  imm << " " << trkid << " " << iid << " : " << pdg << " " << pPart->PdgCode() << " " << mom_pdg << " " << pPart->PdgCode() << std::endl;
					mom_trkid = iid;
					mom_pdg = pPart->PdgCode();
					break;
				}
			}
		}
	        int pdg2 = -1;
		int trkid2 = -1;
		int mom_trkid2 = -1;
		int mom_pdg2 = -1;
		if (shwtrue.size()>1){
			pdg2 = shwtrue[1].first->PdgCode();
			trkid2 = shwtrue[1].first->TrackId();
                	mom_pdg2 = shwtrue[1].first->Mother()==0 ? -1 : trueParticles[shwtrue[1].first->Mother()]->PdgCode();
			for (auto itPart: trueParticles){
				const simb::MCParticle* pPart = itPart.second;
				if (pPart->Process() != "primary") continue;
				int iid = pPart->TrackId();
				if (trkid2==iid) continue;
				for (int imm = 0; imm < pPart->NumberDaughters(); imm++){
					int jid = pPart->Daughter(imm);
					if (jid == trkid2){
						mom_trkid2 = iid;
						mom_pdg2 = pPart->PdgCode();
						break;
					}
				}
			}
		} // end of shwtrue.size()>1


	        if (true){
	            TVector3 v3_true(shwtrue[0].first->Momentum().Vect());
	            float dang_true = (float)(vv1.Angle(v3_true));
	            if (abs(pdg) == 11){
	                electron_candidate = ishower;
	                e_candidate_id = emshower->ID();
	            }
	            all_shw_true_Efrac[ishower] = shwtrue[0].second;
	            all_shw_true_eng[ishower] = shwtrue[0].first->Momentum().E();
	            all_shw_true_pdg[ishower] = pdg;
	            all_shw_true_pdg_mom[ishower] = mom_pdg;
	            all_shw_min_angle[ishower] = dang_true;
	            all_shw_true_px[ishower] = v3_true.X();
	            all_shw_true_py[ishower] = v3_true.Y();
	            all_shw_true_pz[ishower] = v3_true.Z();
		    if (mom_trkid>0){
	  	          all_shw_true_mom_px[ishower] = trueParticles[mom_trkid]->Momentum().Vect().X();
	  	          all_shw_true_mom_py[ishower] = trueParticles[mom_trkid]->Momentum().Vect().Y();
	  	          all_shw_true_mom_pz[ishower] = trueParticles[mom_trkid]->Momentum().Vect().Z();
		    }
	            std::cout << ishower << " PDG " << pdg << " " << mom_pdg << std::endl;
	            std::cout << "Angle: " << dang_true*TMath::RadToDeg() << " " << all_shw_eng[ishower] << std::endl;
		    v3_true.Print();
		    vv1.Print();
		    if (shwtrue.size()>1){
	                TVector3 v3_true(shwtrue[1].first->Momentum().Vect());
	                all_shw2_true_Efrac[ishower] = shwtrue[1].second;
	                all_shw2_true_eng[ishower] = shwtrue[1].first->Momentum().E();
	                all_shw2_true_pdg[ishower] = pdg2;
	                all_shw2_true_pdg_mom[ishower] = mom_pdg2;
	                all_shw2_true_px[ishower] = v3_true.X();
	                all_shw2_true_py[ishower] = v3_true.Y();
	                all_shw2_true_pz[ishower] = v3_true.Z();
		        if (mom_trkid2>0){
	  	              all_shw2_true_mom_px[ishower] = trueParticles[mom_trkid2]->Momentum().Vect().X();
	  	              all_shw2_true_mom_py[ishower] = trueParticles[mom_trkid2]->Momentum().Vect().Y();
	  	              all_shw2_true_mom_pz[ishower] = trueParticles[mom_trkid2]->Momentum().Vect().Z();
		        }
		    } // end of shwtrue.size()>1

		}
	    }


           for (size_t h = 0; h < vhit.size(); ++h){
             if (h<kMaxHits){
               double correctedHitCharge = vhit[h]->Integral() * fCaloAlg.LifetimeCorrection(clockData, detPropEvt, vhit[h]->PeakTime(), det_t0);
               if (vhit[h]->WireID().Plane == 2){
                 showerCharge += (float) correctedHitCharge;
		 shw_chg[ishower] += correctedHitCharge;
		 shw_nhits[ishower] += 1;
               } // if plane==2
               //-------------------------------------------------------------------

               switch (vhit[h]->WireID().Plane) {
               case 0:
                 shw_ChargeU += vhit[h]->Integral();
                 shw_correctedChargeU += correctedHitCharge;
                 break;
               case 1:
                 shw_ChargeV += vhit[h]->Integral();
                 shw_correctedChargeV += correctedHitCharge;
                 break;
               case 2:
                 shw_ChargeZ += vhit[h]->Integral();
                 shw_correctedChargeZ += correctedHitCharge;
                 break;
               }

               // Fill hit level info
	       if ((int)vhit[h]->WireID().TPC != tmp_tpc){
		   tmp_tpc = (int)vhit[h]->WireID().TPC;
		   ntpcs ++;
  	       }
               hit_shw_pdg[hit_shw_count]     = all_shw_true_pdg[ishower];
               hit_shw_pdg_mom[hit_shw_count] = all_shw_true_pdg_mom[ishower];
               hit_shw_tpc[hit_shw_count]     = vhit[h]->WireID().TPC;
               hit_shw_plane[hit_shw_count]   = vhit[h]->WireID().Plane;
               hit_shw_wire[hit_shw_count]    = vhit[h]->WireID().Wire ;
               hit_shw_peakT[hit_shw_count]   = vhit[h]->PeakTime();
               hit_shw_startT[hit_shw_count]  = vhit[h]->StartTick() ;
               hit_shw_endT[hit_shw_count]    = vhit[h]->EndTick() ;
               hit_shw_charge[hit_shw_count]  = vhit[h]->Integral() ;
               hit_shw_channel[hit_shw_count] = vhit[h]->Channel() ;
               hit_shw_multi[hit_shw_count]   = vhit[h]->Multiplicity() ;
               hit_shw_corrcharge[hit_shw_count] = correctedHitCharge ;
               hit_shw_energy[hit_shw_count]  = ((correctedHitCharge * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;

               // Fill true info
               //hit_shw_ntrackID[hit_shw_count] =  TrackIDs.size() ;
               //hit_shw_trueE[hit_shw_count]    = sum_energy ;
               hit_shw_ntrackID[hit_shw_count] =  -1 ;
               hit_shw_trueE[hit_shw_count]    = -1 ;
	       // fill shower level
               shw_idx[hit_shw_count]  = ishower;
               shw_dirX[hit_shw_count] = emshower->Direction().X();
               shw_dirY[hit_shw_count] = emshower->Direction().Y();
               shw_dirZ[hit_shw_count] = emshower->Direction().Z();
               shw_startX[hit_shw_count] = emshower->ShowerStart().X();
               shw_startY[hit_shw_count] = emshower->ShowerStart().Y();
               shw_startZ[hit_shw_count] = emshower->ShowerStart().Z();

	       unsigned int globalWire = vhit[h]->WireID().Wire;
	       unsigned int globalPlane = vhit[h]->WireID().Plane;
	       double globalTDC = vhit[h]->PeakTime();
	       GetDUNEGlobalWireTDC(vhit[h]->WireID(), (double)vhit[h]->PeakTime(), globalWire, globalPlane, globalTDC, detPropEvt.DriftVelocity());
               hit_shw_global_plane[hit_shw_count]   = (double)globalPlane;
               hit_shw_global_wire[hit_shw_count]    = (double)globalWire;
               hit_shw_global_tick[hit_shw_count]    = globalTDC;


	       hit_shw_count ++;
               //-------------------------------------------------------------------
             } // end of kMaxHits
           } // end of shower hit loop


 	   double shw_eng_zplane = emshower->Energy()[2];
 	   if (shw_eng_zplane > max_shw_energy){
	          max_shw_energy = shw_eng_zplane;
		  led_shw_idx  = ishower;
		  led_shw_id   = emshower->ID();
		  double tmp_shw_eng = ((fCaloAlg.ElectronsFromADCArea(max_shw_energy, 2) * (1.0 / 0.63) * util::kGeVToElectrons) - fIntShwEnergy) / fGradShwEnergy ;
		  tmp_shw_eng = std::sqrt(tmp_shw_eng*tmp_shw_eng+0.0005109989461*0.0005109989461);

	          led_shw_eng  = max_shw_energy; // fixme
                  led_shw_dirx = emshower->Direction().X();
                  led_shw_diry = emshower->Direction().Y();
                  led_shw_dirz = emshower->Direction().Z();
                  led_shw_startx = emshower->ShowerStart().X();
                  led_shw_starty = emshower->ShowerStart().Y();
                  led_shw_startz = emshower->ShowerStart().Z();
		  led_shw_costh  = cos(vv1.Theta());


   	          // get unique primary particles
	          //std::vector<int> uni_trackid = getUniques(trackid_col);
	          //std::vector<int> count_trackid;
	          //count_trackid.resize(uni_trackid.size());

		  std::cout << "Reco Shw Dir: " << all_shw_costh[ishower] << " " << all_shw_phi[ishower] << std::endl;
		  std::cout << "Reco vs Pri0  " << vv1.Angle(trueParticle->Momentum().Vect())*TMath::RadToDeg() << std::endl;
	          //std::cout << "size of uni/all trackids : " << uni_trackid.size() << " / " << trackid_col.size() << std::endl;
		  if (shwtrue.size()>0){
    	              led_shw_trueP0 = shwtrue[0].first->Momentum().Px();
       	              led_shw_trueP1 = shwtrue[0].first->Momentum().Py();
       	              led_shw_trueP2 = shwtrue[0].first->Momentum().Pz();
       	              led_shw_trueP3 = shwtrue[0].first->Momentum().E();
       	              led_shw_truedang  = vv1.Angle(shwtrue[0].first->Momentum().Vect());
       	              led_shw_truePDG   = shwtrue[0].first->PdgCode();
       	              led_shw_truePDGMom = all_shw_true_pdg_mom[ishower];
	              led_shw_truecosth = cos(shwtrue[0].first->Momentum().Vect().Theta());
		  }

	   } // collect leading emshower

         } // end of fmsh
	 if (showerCharge>0) icount_shower += 1;
         shw_tot_chg = showerCharge;
     } // end of loop for ishower
  } // end of if showerHandle

  shw_correctedEnergyU = ((shw_correctedChargeU * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  shw_correctedEnergyV = ((shw_correctedChargeV * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  shw_correctedEnergyZ = ((shw_correctedChargeZ * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  shw_EnergyU = ((shw_ChargeU * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  shw_EnergyV = ((shw_ChargeV * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  shw_EnergyZ = ((shw_ChargeZ * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  // end of emshower
  //-----------------------------------------------------------------------------------------


  // Lifetime-corrected charge
  correctedChargeU = 0;
  correctedChargeV = 0;
  correctedChargeZ = 0;
  ChargeU = 0;
  ChargeV = 0;
  ChargeZ = 0;

  tmp_tpc = -10;
  ntpcs_evt = 0;

  float max_tick = 0;
  float min_tick = 1e5;
  double max_rms = 0;
  int temp_wireID = -99999;

  int max_wireID = 0;
  int min_wireID = 1e5;
  float max_ledshwtick = -1;
  float min_ledshwtick = 1e5;

  // Look at the hits
  n_wires = 0;
  std::cout << e_candidate_id << min_wireID << max_wireID << std::endl;

  art::FindManyP<recob::Shower> fmshwhit(hitHandle, evt, fShowerModuleLabel);
  art::FindManyP<recob::Track> fmtrkhit(hitHandle, evt, fTrackModuleLabel2);
  //shower::EMShowerAlg emshoweralg;
  std::cout << "Loop Hits: " << hits.size() << std::endl;
  trkide.clear();
  nhits = 0;
  // loop hit level first
  for (unsigned int hitIt = 0; hitIt < hits.size(); ++hitIt) {
    if (hitIt >= kMaxHits) break;
    nhits += 1;
    // Get the hit
    art::Ptr<recob::Hit> hit = hits.at(hitIt);
    std::vector<art::Ptr<recob::SpacePoint>> sp = spFromHits.at(hitIt);

    if(!sp.empty()) {
        hit_sp_x[hitIt] = sp[0]->XYZ()[0];
        hit_sp_y[hitIt] = sp[0]->XYZ()[1];
        hit_sp_z[hitIt] = sp[0]->XYZ()[2];
    }

    double correctedHitCharge = ( hit->Integral() * TMath::Exp( (SampleRate * hit->PeakTime()) / (detPropEvt.ElectronLifetime()*1e3) ) );
    switch (hit->WireID().Plane) {
    case 0:
      ChargeU += hit->Integral();
      correctedChargeU += correctedHitCharge;
      break;
    case 1:
      ChargeV += hit->Integral();
      correctedChargeV += correctedHitCharge;
      break;
    case 2:
      ChargeZ += hit->Integral();
      correctedChargeZ += correctedHitCharge;
      break;
    }

    // Fill hit level info
    hit_cryst   [hitIt] = hit->WireID().Cryostat;
    hit_tpc     [hitIt] = hit->WireID().TPC;
    hit_plane   [hitIt] = hit->WireID().Plane;
    hit_wire    [hitIt] = hit->WireID().Wire;
    hit_peakT   [hitIt] = hit->PeakTime();
    hit_rmsT    [hitIt] = hit->RMS();
    hit_startT  [hitIt] = hit->StartTick();
    hit_endT    [hitIt] = hit->EndTick();
    hit_charge  [hitIt] = hit->Integral();
    hit_channel [hitIt] = hit->Channel();
    hit_multi   [hitIt] = hit->Multiplicity();
    hit_amp     [hitIt] = hit->PeakAmplitude();
    hit_minusT  [hitIt] = hit->PeakTimeMinusRMS();
    hit_plusT   [hitIt] = hit->PeakTimePlusRMS();
    hit_good    [hitIt] = hit->GoodnessOfFit();
    hit_ndf     [hitIt] = hit->DegreesOfFreedom();

    hit_corrcharge  [hitIt] = correctedHitCharge;
    hit_energy  [hitIt] = ((correctedHitCharge * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy;




  // Workspace geometry has two drift regions
  //   //                  |-----|-----| /  /
       //      y ^         |  3  |  2  |/  /
       //        | -| z    |-----|-----|  /
       //        | /       |  1  |  0  | /
       //  x <---|/        |-----|-----|/
      //

    if (fGlobalWireMethod == 1){
	    hit_global_wire[hitIt] = GetPixelMapWireIdx(hit->WireID());
    	    hit_global_tick[hitIt] = (hit_tpc[hitIt]%2==0) ? -(double)hit_peakT[hitIt] : (double)hit_peakT[hitIt];
	    hit_global_plane[hitIt] = (double) hit->WireID().Plane;
    } else if (fGlobalWireMethod == 2){
	    unsigned int globalWire = hit->WireID().Wire;
	    unsigned int globalPlane = hit->WireID().Plane;
	    double globalTDC = hit->PeakTime();
	    GetDUNEGlobalWireTDC(hit->WireID(), (double)hit->PeakTime(), globalWire, globalPlane, globalTDC, detPropEvt.DriftVelocity());
	    hit_global_wire[hitIt]  = (double) globalWire;
	    hit_global_tick[hitIt]  = globalTDC;
	    //hit_global_tick[hitIt]  = (hit_tpc[hitIt]%2==0) ? -globalTDC : -globalTDC;
	    hit_global_plane[hitIt] = (double) globalPlane;
    } else {
	   std::cout << "Wrong Global Wire Method" << std::endl;
    }
    // 09/11 here add typecast (double)
    //std::cout<<"fGlobalWireMethod, hit_peakT[hitIt], hit_global_wire[hitIt]  "<<fGlobalWireMethod<<" "<<hit_peakT[hitIt]<<" "<<hit_global_wire[hitIt]<<std::endl;
    // Fill true info
    hit_ntrackID[hitIt] =  -1;
    hit_trueE   [hitIt] =  -1;

    //auto xyz = geom -> Cryostat(geo::CryostatID(hit->WireID().Cryostat)).TPC(hit->WireID().TPC).Plane(hit->WireID().Plane).Wire(hit->WireID().Wire).GetCenter();
    //auto xyz = m_geom->WireIDToWireGeo(hit->WireID()).GetCenter(); 
    auto xyz = wireReadout->Wire(hit->WireID()).GetCenter();

    hit_y[hitIt] = xyz.Y(); // cm
    hit_z[hitIt] = xyz.Z(); // cm
    hit_x[hitIt] = detPropEvt.ConvertTicksToX(hit_peakT[hitIt], hit_plane[hitIt], hit_tpc[hitIt], hit->WireID().Cryostat);

   if ((int)hit->WireID().TPC != tmp_tpc){
	tmp_tpc = (int)hit->WireID().TPC;
        ntpcs_evt ++;
    }
    if (temp_wireID < -9999) temp_wireID = hit_wire[hitIt];
    if (temp_wireID != hit_wire[hitIt]){ // make index for wires
	    wire_assn[n_wires] = (int)(hitIt-1);
	    wire_min_tick[n_wires] = (int)min_tick;
	    wire_max_tick[n_wires] = (int)max_tick;
	    wire_max_rmsT[n_wires] = (float)max_rms;
	    if (max_ledshwtick >=0){
	      wire_min_ledshwtick[n_wires] = (int)round(min_ledshwtick);
	      wire_max_ledshwtick[n_wires] = (int)round(max_ledshwtick);
	    }
	    for (int iclst = 0; iclst < kMaxClst; iclst++){
		    wire_min_shwtick[n_wires][iclst] = (int)round(min_shwtick[iclst]);
		    wire_max_shwtick[n_wires][iclst] = (int)round(max_shwtick[iclst]);
		    wire_min_trktick[n_wires][iclst] = (int)round(min_trktick[iclst]);
		    wire_max_trktick[n_wires][iclst] = (int)round(max_trktick[iclst]);
		    wire_max_shwrms[n_wires][iclst] = (int)round(max_shwrms[iclst]);
		    wire_max_trkrms[n_wires][iclst] = (int)round(max_trkrms[iclst]);
		    max_shwtick[iclst] = -999;
		    min_shwtick[iclst] = 1e+5;
		    max_trktick[iclst] = -999;
		    min_trktick[iclst] = 1e+5;
		    max_shwrms[iclst] = -999;
		    max_trkrms[iclst] = -999;
	    }

	    n_wires ++;

	    temp_wireID = hit_wire[hitIt];
	    max_wireID = 0; min_wireID = 1e5; max_tick = 0; min_tick = 1e5;
	    max_ledshwtick = -1; min_ledshwtick = 1e5;
	    max_rms = 0;
    }

    if (hit_peakT[hitIt] > max_tick) {max_tick = hit_peakT[hitIt];}
    if (hit_peakT[hitIt] < min_tick) {min_tick = hit_peakT[hitIt];}

    if (hit_wire[hitIt] > max_wireID) max_wireID = hit_wire[hitIt];
    if (hit_wire[hitIt] < min_wireID) min_wireID = hit_wire[hitIt];
    if (hit_rmsT[hitIt] > max_rms) max_rms = (double)hit_rmsT[hitIt];


    // Find the true track this hit is associated with
    //hit_truetrackid[hitIt] = this->FindTrackID(hit);
    hit_truetrackid[hitIt] = -1;

    const simb::MCParticle *truemc = get_bestMCParticle(hit, clockData);
    int pdg = -10;
    int mom_pdg = -10;
    if (truemc){
	pdg = (int)abs(truemc->PdgCode());
	mom_pdg = truemc->Mother() == 0 ? -1 : trueParticles[truemc->Mother()]->PdgCode();
    }

   hit_prong_pdg[hitIt] = pdg;
   hit_prong_pdg_mom[hitIt] = mom_pdg;
   if (showerHandle.isValid()){
       if (fmshwhit.isValid() && fmshwhit.at(hitIt).size()!=0){
           hit_prong_tag[hitIt] = fmshwhit.at(hitIt)[0]->ID();
	   if (fmshwhit.at(hitIt)[0]->ID() < kMaxClst){
	       if (hit_peakT[hitIt] > max_shwtick[fmshwhit.at(hitIt)[0]->ID()])
		       max_shwtick[fmshwhit.at(hitIt)[0]->ID()] = hit_peakT[hitIt];
	       if (hit_peakT[hitIt] < min_shwtick[fmshwhit.at(hitIt)[0]->ID()])
		       min_shwtick[fmshwhit.at(hitIt)[0]->ID()] = hit_peakT[hitIt];
	       if (hit_rmsT[hitIt] > max_shwrms[fmshwhit.at(hitIt)[0]->ID()])
		       max_shwrms[fmshwhit.at(hitIt)[0]->ID()] = hit_rmsT[hitIt];
	   }

       }
   }
   if (trackHandle2.isValid()){
       if (fmtrkhit.isValid() && fmtrkhit.at(hitIt).size()!=0){
           hit_prong_tag[hitIt] = kMaxClst+fmtrkhit.at(hitIt)[0]->ID();
	   if (fmtrkhit.at(hitIt)[0]->ID() < kMaxClst){
	       if (hit_peakT[hitIt] > max_trktick[fmtrkhit.at(hitIt)[0]->ID()])
		       max_trktick[fmtrkhit.at(hitIt)[0]->ID()] = hit_peakT[hitIt];
	       if (hit_peakT[hitIt] < min_trktick[fmtrkhit.at(hitIt)[0]->ID()])
		       min_trktick[fmtrkhit.at(hitIt)[0]->ID()] = hit_peakT[hitIt];
	       if (hit_rmsT[hitIt] > max_trkrms[fmtrkhit.at(hitIt)[0]->ID()])
		       max_trkrms[fmtrkhit.at(hitIt)[0]->ID()] = hit_rmsT[hitIt];
	   }

       }
   }

  } // end of hitIt
  std::cout << "End of Hit Loop " << std::endl;
  //Shift Global Wire Idx
  if (fGlobalWireMethod == 1){
    ShiftGlobalWireIdx(hit_global_wire, hit_global_tick, hit_tpc, hit_plane, hit_offset);
  }

  correctedEnergyU = ((correctedChargeU * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  correctedEnergyV = ((correctedChargeV * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  correctedEnergyZ = ((correctedChargeZ * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  EnergyU = ((ChargeU * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  EnergyV = ((ChargeV * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;
  EnergyZ = ((ChargeZ * (1.0 / 0.63) * (23.6e-9 / 4.966e-3)) - fIntShwEnergy) / fGradShwEnergy ;


  // Put energies in GeV units
  depositU /= 1000;
  depositV /= 1000;
  depositZ /= 1000;

  // Get RecoE from DUNE
  if (!engrecoHandle.failedToGet())
  {
      ErecoNu          = engrecoHandle->fNuLorentzVector.E();
      RecoLepEnNu      = engrecoHandle->fLepLorentzVector.E();
      RecoHadEnNu      = engrecoHandle->fHadLorentzVector.E();
      RecoMethodNu     = engrecoHandle->recoMethodUsed;
      std::cout << "Reco nuE: " << ErecoNu << std::endl;
      std::cout << "Reco lepE: " << RecoLepEnNu << std::endl;
      std::cout << "Reco hadE: " << RecoHadEnNu << std::endl;
      std::cout << "Reco Track Cont: " << engrecoHandle->longestTrackContained << std::endl;
  }


  //-----------------------------------------------------------------------------------------
  // find vertex of leading shower
   if (showerHandle.isValid()){
     double shw_xyz[3] = {led_shw_startx, led_shw_starty, led_shw_startz};
     inTPC = ConvertXYZtoWireTDC(shw_xyz, shwVtxTPC, shwVtxWire, shwVtxTime,
 		     shwVtxGlobalPlane, shwVtxGlobalWire, shwVtxGlobalTick, detPropEvt);
     //if (inTPC) FindShowerVtx();
  } // end of showerHandle

  // find positions on wire and tdc

  //std::cout << "Loop Wires" << std::endl;
  //-----------------------------------------------------------------------------------------
  // loop based on hit wires
  for (int iwire = 0; iwire < n_wires; iwire++){
    if (iwire >= kMaxHits) continue;
    int hitIt = wire_assn[iwire];

    // extract more lower level
    int t0_hit = 0;
    int t1_hit = 0;

    // extract rawDigit per hit
    if (fmraw.isValid()){
        std::vector< art::Ptr<raw::RawDigit> > allraws = fmraw.at(hitIt);
	hit_Nallraws[iwire] = allraws.size();
        raw::RawDigit::ADCvector_t adcvec(allraws[0]->Samples());
        raw::Uncompress(allraws[0]->ADCs(), adcvec, allraws[0]->Compression());
	unsigned int chNum = allraws[0]->Channel();
        std::vector< geo::WireID > wireID =
                                    wireReadout -> ChannelToWire(chNum);

	hit_raw_chID[iwire]   = chNum;
	hit_raw_wireID[iwire] = wireID[0].Wire;
        float hit_first_time  = hit_peakT[hitIt] - 3*(hits.at(hitIt)->RMS());
        float hit_end_time    = hit_peakT[hitIt] + 3*(hits.at(hitIt)->RMS());
        t0_hit = (hit_first_time < 0 ) ? 0    : (int)hit_first_time;
        t1_hit = (hit_end_time > 4491) ? 4491 : (int)hit_end_time;
	hit_Nadcvec[iwire] = t1_hit-t0_hit+1;

	std::vector<int> tmpadc;
	std::vector<int> tmptick;
        for (int j = t0_hit; j<=t1_hit; ++j){
	  if (adcvec[j] == 0) continue;
	    int adc_ped = int(adcvec[j])-int(allraws[0]->GetPedestal());
	    tmpadc.push_back(adc_ped);
	    tmptick.push_back(j);
        }
	hit_rawadc.push_back(tmpadc);
	hit_rawtick.push_back(tmptick);
	tmpadc.clear();
        tmptick.clear();

    } // end of fmraw


    // extract reco wire per hit
    if (fmwire.isValid()){
        std::vector< art::Ptr<recob::Wire> > wireptr = fmwire.at(hitIt);

	for (size_t iwireptr = 0; iwireptr < wireptr.size(); iwireptr++){
		std::vector<geo::WireID> wireids = wireReadout->ChannelToWire(wireptr[iwireptr]->Channel());
		bool goodWID = false;
		for (auto const & wid:wireids){
			if (wid.Plane == hits.at(hitIt)->WireID().Plane &&
			    wid.Wire  == hits.at(hitIt)->WireID().Wire &&
			    wid.TPC   == hits.at(hitIt)->WireID().TPC &&
			    wid.Cryostat == hits.at(hitIt)->WireID().Cryostat) goodWID = true;
		}
		if (!goodWID) continue;
	        const recob::Wire::RegionsOfInterest_t& signalROI = wireptr[iwireptr]->SignalROI();

	        raw::TDCtick_t roiFirstBinTick = signalROI.get_ranges()[iwireptr].begin_index();
	        raw::TDCtick_t roiEndBinTick   = signalROI.get_ranges()[iwireptr].end_index();
	        wire_roifirst[iwire] = roiFirstBinTick;
	        wire_roiend[iwire] = roiEndBinTick;

	        std::vector<double> tmp_wireadc;
	        std::vector<double> tmp_wirecorradc;
	        std::vector<int> tmp_wiretick;
	        std::vector<double> tmp_wire_x;
	        std::vector<int> tmp_shwflag;
	        std::vector<int> tmp_ledshwflag;
	        std::vector<double> tmp_dr;
	        std::vector<double> tmp_dist;
	        std::vector<int> tmp_cand_pdg;
	        std::vector<int> tmp_cand_pdg_mom;
	        std::vector<int> tmp_tag_pdg;
	        std::vector<int> tmp_tag_pdg_mom;
	        std::vector<int> tmp_tag_is_shw;
	        std::vector<int> tmp_tag_ith;
	        std::vector<double> tmp_tag_true_px;
	        std::vector<double> tmp_tag_true_py;
	        std::vector<double> tmp_tag_true_pz;
	        std::vector<double> tmp_tag_reco_px;
	        std::vector<double> tmp_tag_reco_py;
	        std::vector<double> tmp_tag_reco_pz;
	        std::vector<double> tmp_tag_reco_eng;

	        const std::vector<float>& signal = wireptr[0]->Signal();

	        float hit_first_time  = wire_min_tick[iwire]- 3*(wire_max_rmsT[iwire]);
		float hit_end_time    = wire_max_tick[iwire]+ 3*(wire_max_rmsT[iwire]);
	        t0_hit = (hit_first_time < 0 ) ? 0    : (int)hit_first_time;
	        t1_hit = (hit_end_time > 4491) ? 4491 : (int)hit_end_time;

	        //for (int j = roiFirstBinTick; j <= roiEndBinTick; ++j){
	        for (int j = t0_hit; j <= t1_hit; ++j){
	        //for (unsigned int j = 0; j < signal.size(); ++j){
		   unsigned int globalWire  = hits.at(hitIt)->WireID().Wire;
		   unsigned int globalPlane = hits.at(hitIt)->WireID().Plane;
		   //unsigned int localWire = hits.at(hitIt)->WireID().Wire;
		   //unsigned int localPlane = hits.at(hitIt)->WireID().Plane;
		   double globalTDC = (double)j;
    		   double correctedadc = ( signal[j] * TMath::Exp( (SampleRate * j) / (detPropEvt.ElectronLifetime()*1e3) ) );

	    	   GetDUNEGlobalWireTDC(hits.at(hitIt)->WireID(), (double)j, globalWire, globalPlane, globalTDC, detPropEvt.DriftVelocity());
		   if (fGlobalWireMethod == 1){
	  	     tmp_wiretick.push_back(j);
		   } else if (fGlobalWireMethod == 2){
	  	     tmp_wiretick.push_back((int)globalTDC);
		   } else { std::cout << "no... way.." << std::endl;}

		   //------------------------------------------------------------------
		   float min_d_tick = 1e+6;
		   int ptag = -1;
		   for (int ishower = 0; ishower < n_showers; ++ishower){
			if (ishower < kMaxClst){
		            if (wire_min_shwtick[iwire][ishower]-wire_max_shwrms[iwire][ishower]*5 <= j && j <= wire_max_shwtick[iwire][ishower]+wire_max_shwrms[iwire][ishower]*5){
				    float d_tick = (float)std::max(abs(j-wire_min_shwtick[iwire][ishower]),abs(j-wire_max_shwtick[iwire][ishower]));
				    //float d_tick = (float)(j-(wire_max_shwtick[iwire][ishower]+wire_min_shwtick[iwire][ishower])/2);
				    if (d_tick < min_d_tick){
					    min_d_tick = d_tick;
				            ptag = ishower;
				    }
		            }
			}
		   } // ishowers
      		   for (int itrack = 0; itrack < n_tracks_pand; itrack++){
			if (itrack < kMaxClst){
		            if (wire_min_trktick[iwire][itrack]-wire_max_trkrms[iwire][itrack]*5 <= j && j <= wire_max_trktick[iwire][itrack]+wire_max_trkrms[iwire][itrack]*5){
				    float d_tick = (float)std::max(abs(j-wire_min_trktick[iwire][itrack]),abs(j-wire_max_trktick[iwire][itrack]));
				    //float d_tick = (float)(j-(wire_max_trktick[iwire][itrack]+wire_min_trktick[iwire][itrack])/2);
				    if (d_tick < min_d_tick){
					    min_d_tick = d_tick;
				            ptag = itrack+kMaxClst;
				    }
		            }
			}

		   } // itrack
		   if (ptag>-1){
			   if (ptag<kMaxClst){
		               tmp_shwflag.push_back(1);
		               tmp_tag_ith.push_back(ptag);
		               tmp_tag_is_shw.push_back(1);
		               tmp_tag_pdg.push_back(all_shw_true_pdg[ptag]);
		               tmp_tag_pdg_mom.push_back(all_shw_true_pdg_mom[ptag]);
		               tmp_tag_true_px.push_back(all_shw_true_px[ptag]);
		               tmp_tag_true_py.push_back(all_shw_true_py[ptag]);
		               tmp_tag_true_pz.push_back(all_shw_true_pz[ptag]);
		               tmp_tag_reco_px.push_back(all_shw_px[ptag]);
		               tmp_tag_reco_py.push_back(all_shw_py[ptag]);
		               tmp_tag_reco_pz.push_back(all_shw_pz[ptag]);
		               tmp_tag_reco_eng.push_back(all_shw_eng[ptag]);
			   } else{
		               tmp_shwflag.push_back(0);
		               tmp_tag_ith.push_back(ptag-kMaxClst);
		               tmp_tag_is_shw.push_back(0);
		               tmp_tag_pdg.push_back(all_track_true_pdg[ptag-kMaxClst]);
		               tmp_tag_pdg_mom.push_back(all_track_true_pdg_mom[ptag-kMaxClst]);
		               tmp_tag_true_px.push_back(all_track_true_px[ptag-kMaxClst]);
		               tmp_tag_true_py.push_back(all_track_true_py[ptag-kMaxClst]);
		               tmp_tag_true_pz.push_back(all_track_true_pz[ptag-kMaxClst]);
		               tmp_tag_reco_px.push_back(all_track_px[ptag-kMaxClst]);
		               tmp_tag_reco_py.push_back(all_track_py[ptag-kMaxClst]);
		               tmp_tag_reco_pz.push_back(all_track_pz[ptag-kMaxClst]);
		               tmp_tag_reco_eng.push_back(all_track_calmom[ptag]);
			   }
		   }else{
	               tmp_shwflag.push_back(-10);
	               tmp_tag_ith.push_back(-10);
		       tmp_tag_pdg.push_back(-10);
		       tmp_tag_pdg_mom.push_back(-10);
		       tmp_tag_is_shw.push_back(-10);
		       tmp_tag_true_px.push_back(-10);
		       tmp_tag_true_py.push_back(-10);
		       tmp_tag_true_pz.push_back(-10);
		       tmp_tag_reco_px.push_back(-10);
		       tmp_tag_reco_py.push_back(-10);
		       tmp_tag_reco_pz.push_back(-10);
	               tmp_tag_reco_eng.push_back(-10);
		   }
		   //------------------------------------------------------------------

	  	   tmp_wireadc.push_back((double)signal[j]);
	  	   tmp_wirecorradc.push_back(correctedadc);
	           tmp_wire_x.push_back( detPropEvt.ConvertTicksToX(j, hit_plane[hitIt], hit_tpc[hitIt], hit_cryst[hitIt]) );
		   if (wire_min_ledshwtick[iwire] <= j && j <= wire_max_ledshwtick[iwire]){
			   tmp_ledshwflag.push_back(1);
		   }
		   else tmp_ledshwflag.push_back(0);

	        } // end of ADCcount
	        wire_adc.push_back(tmp_wireadc);
	        wire_corradc.push_back(tmp_wirecorradc);
	        wire_tick.push_back(tmp_wiretick);
	        wire_x.push_back(tmp_wire_x);
		wire_shwflag.push_back(tmp_shwflag);
		wire_ledshwflag.push_back(tmp_ledshwflag);
		wire_pdg.push_back(tmp_cand_pdg);
		wire_pdg_mom.push_back(tmp_cand_pdg_mom);
		wire_dr.push_back(tmp_dr);
		wire_dist.push_back(tmp_dist);
		wire_tag_pdg.push_back(tmp_tag_pdg);
		wire_tag_pdg_mom.push_back(tmp_tag_pdg_mom);
		wire_tag_is_shw.push_back(tmp_tag_is_shw);
		wire_tag_ith.push_back(tmp_tag_ith);
		wire_tag_true_px.push_back(tmp_tag_true_px);
		wire_tag_true_py.push_back(tmp_tag_true_py);
		wire_tag_true_pz.push_back(tmp_tag_true_pz);
		wire_tag_reco_px.push_back(tmp_tag_reco_px);
		wire_tag_reco_py.push_back(tmp_tag_reco_py);
		wire_tag_reco_pz.push_back(tmp_tag_reco_pz);
		wire_tag_reco_eng.push_back(tmp_tag_reco_eng);


		tmp_wireadc.clear();
		tmp_wirecorradc.clear();
		tmp_wiretick.clear();
		tmp_wire_x.clear();
		tmp_shwflag.clear();
		tmp_ledshwflag.clear();
		tmp_cand_pdg.clear();
		tmp_cand_pdg_mom.clear();
		tmp_dr.clear();
		tmp_dist.clear();
		tmp_tag_pdg.clear();
		tmp_tag_pdg_mom.clear();
		tmp_tag_is_shw.clear();
		tmp_tag_ith.clear();
		tmp_tag_true_px.clear();
		tmp_tag_true_py.clear();
		tmp_tag_true_pz.clear();
		tmp_tag_reco_px.clear();
		tmp_tag_reco_py.clear();
		tmp_tag_reco_pz.clear();
		tmp_tag_reco_eng.clear();

		// simChannel
	        std::vector<int> tmp_chtick;
	        std::vector<double> tmp_chchg;
	        for (int ii  = 0; ii < simCh_size; ++ii){
	            if (simChannels[ii]->Channel() != hits.at(hitIt)->Channel() ) continue;
	            for (int j = t0_hit; j <= t1_hit; ++j){
	    	       tmp_chtick.push_back(j);
	    	       tmp_chchg.push_back(simChannels[ii]->Charge(j));
	               sum_ch_chg[iwire] += simChannels[ii]->Charge(j);
	            }
	        } // simCh_size
	        ch_chg.push_back(tmp_chchg);
        	ch_tick.push_back(tmp_chtick);
	} // end of iwireptr for wireptr
    } // end of fmwire

  } // end of n_wires

  // end of fmwire
  //-----------------------------------------------------------------------------------------

  int rawcrys = 0;
  int rawtpc  = -1;
  if (vtx_in_tpc) {
    std::cout << "True Particle T: " << trueParticle->T() << "  " << detPropEvt.DriftVelocity() << std::endl;
    pos =  vertex + ( (1.0) * direction );
    double currentPos[3]; currentPos[0] = pos.X(); currentPos[1] = pos.Y(); currentPos[2] = pos.Z();
    geo::Point_t pos(currentPos[0], currentPos[1], currentPos[2]);
    if (geom->FindTPCAtPosition(pos).isValid){
	for (int iplane = 0; iplane < 3; iplane++){
          rawtpc = (int) (geom->FindTPCAtPosition(pos)).TPC;
    geo::PlaneID planeID{(unsigned int) rawcrys, (unsigned int) rawtpc, (unsigned int) iplane};
	  geo::PlaneGeo const& planegeo_temp = wireReadout->Plane(planeID);
	  geo::WireID w1;
	  try {
  	        w1 = wireReadout->Plane(planeID).NearestWireID(pos);
	  }
	  catch (geo::InvalidWireError const& e) {
	        if (!e.hasSuggestedWire()) throw;
	        w1 = planegeo_temp.ClosestWireID(e.suggestedWireID());
	  }
	  double time1 = detPropEvt.ConvertXToTicks(currentPos[0]+trueParticle->T()*detPropEvt.DriftVelocity()*1e-3, iplane, rawtpc, rawcrys);
          trueMidWire[iplane] = w1.Wire;
          trueMidTime[iplane] = time1;
	} // end of iplane
     }
  }


  //std::cout << "Cluster Loop" << std::endl;
  n_clusters = clusters.size();
  fTree->Fill();
  return;

}
int recoestudy::RecoEnergyS::FindShowerVtx(void) {
  TVector3 shw_dir (-led_shw_dirx, 0, -led_shw_dirz);
  TVector3 shw_start(led_shw_startx, 0, led_shw_startz);

  float max_vect_prod = -1e+6;
  float Vertex_Wire = 0;
  float Vertex_PeakT = 0;
  TVector3 pos;
  for (int ihit = 0; ihit < nhits; ihit++){
    if (ihit >= kMaxHits) continue;
    if (hit_plane[ihit] != 2) continue;
    TVector3 start_to_hit_pos (hit_x[ihit]-led_shw_startx, 0, hit_z[ihit]-led_shw_startz);
    float vect_prod = shw_dir * start_to_hit_pos;
    if (vect_prod > max_vect_prod){
	max_vect_prod = vect_prod;
		Vertex_PeakT = hit_peakT[ihit];
		Vertex_Wire  = hit_wire[ihit];
    } // end of vect_prod
  } // end of ihit
  shw_first_wire = Vertex_Wire;
  shw_first_time = Vertex_PeakT;

  return 1;
}


int recoestudy::RecoEnergyS::ShiftGlobalWireIdx(double *global_idx, double *global_tick, int *tpc_idx, int *plane_idx, double *offsets){

	std::vector<int> idxlist_of_zeroT;
	float nhits_Z[24] = {0};
	for (int ii = 0; ii < nhits; ii++){
		if (ii >= kMaxHits) continue;
		if (TMath::Abs(global_tick[ii]) < 150) idxlist_of_zeroT.push_back(ii);
		if (plane_idx[ii] == 2) nhits_Z[tpc_idx[ii]] += 1;
	}

	float temp_avg_U[24] = {0}, temp_avg_V[24] = {0};
	float temp_count_U[24] = {0}, temp_count_V[24] = {0};
	std::vector<int> list_of_tpc;
	int count_tpc = 0;
	for (unsigned int ii = 0; ii < idxlist_of_zeroT.size(); ++ii){
		int idx = idxlist_of_zeroT[ii];
		if (plane_idx[idx] == 2) continue;
		int temp_tpc = tpc_idx[idx];
		double temp_tick = (TMath::Abs(global_tick[idx])+1); // add to avoid tick==0
		if (plane_idx[idx] == 0) {
			temp_avg_U[temp_tpc] += global_idx[idx]/temp_tick;
			temp_count_U[temp_tpc] += 1.0/temp_tick;
		}
		if (plane_idx[idx] == 1) {
			temp_avg_V[temp_tpc] += global_idx[idx]/temp_tick;
			temp_count_V[temp_tpc] += 1.0/temp_tick;
		}

		if (list_of_tpc.size()==0) list_of_tpc.push_back(temp_tpc);
		if (list_of_tpc[count_tpc] != temp_tpc) {
			list_of_tpc.push_back(temp_tpc);
			count_tpc++;
		}

	}
	if (list_of_tpc.size() <= 1) return 0;

	int temp_sum_count = 0;
	int nmax_tpcs[2] = {-1, -1};
	for (unsigned int ii = 0; ii < list_of_tpc.size()-1; ii++){
		int temp_tpc1 = list_of_tpc[ii];
		int temp_tpc2 = list_of_tpc[ii+1];
		if (temp_tpc1/4 != temp_tpc2/4) continue;
		if (temp_tpc1%2 == temp_tpc2%2) continue;
		int sum = nhits_Z[temp_tpc1] + nhits_Z[temp_tpc2];
		if (sum > temp_sum_count){
			nmax_tpcs[0] = temp_tpc1;
			nmax_tpcs[1] = temp_tpc2;
			temp_sum_count = sum;
		}

        }
	std::cout << "Here: " << nmax_tpcs[0] << " " << nmax_tpcs[1] << std::endl;
	if (nmax_tpcs[0]<0 || nmax_tpcs[1]<0) return 0;
	double offsetU = 0;
	double offsetV = 0;
	int temp_tpc1 = nmax_tpcs[0];
	int temp_tpc2 = nmax_tpcs[1];


	if (temp_count_U[temp_tpc1] > 0 && temp_count_U[temp_tpc2] > 0) {
		double avg_wire_U1 = temp_avg_U[temp_tpc1]/temp_count_U[temp_tpc1];
		double avg_wire_U2 = temp_avg_U[temp_tpc2]/temp_count_U[temp_tpc2];
		offsetU = avg_wire_U1-avg_wire_U2;
	}
	if (temp_count_V[temp_tpc1] > 0 && temp_count_V[temp_tpc2] > 0) {
		double avg_wire_V1 = temp_avg_V[temp_tpc1]/temp_count_V[temp_tpc1];
		double avg_wire_V2 = temp_avg_V[temp_tpc2]/temp_count_V[temp_tpc2];
		offsetV = avg_wire_V1-avg_wire_V2;
	}
	std::cout << temp_tpc1 << " " << temp_tpc2 << " " << offsetU << " " << offsetV << std::endl;

	std::cout << list_of_tpc.size() << " " << offsetU << " " << offsetV << std::endl;

	for (int ii = 0; ii < nhits; ii++){
		if (plane_idx[ii] == 2) continue;
		if (tpc_idx[ii]%2 == 1) {
			if (plane_idx[ii] == 0) global_idx[ii] += round(offsetU);
			if (plane_idx[ii] == 1) global_idx[ii] += round(offsetV);
		}
	}
	offsets[0] = round(offsetU);
	offsets[1] = round(offsetV);
	offsets[2] = 0;
        std::cout << offsets[0] << " " << offsets[1] << std::endl;
	return 1;
}

// Based on Robert's code in adcutils
void recoestudy::RecoEnergyS::GetDUNEGlobalWireTDC(const geo::WireID& wireID, double localTDC,
                                           unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC, double driftVel) const
{

  if (fIsVD) {
      globalWire = 0;
      globalPlane = 0;
      globalTDC = localTDC;
      cvn::PixelMapProducer producer;
      producer.GetDUNEVertDrift3ViewGlobalWire(wireID.Wire, wireID.Plane,wireID.TPC,globalWire,globalPlane);
  }
  else {
      unsigned int localWire = wireID.Wire;
      unsigned int plane = wireID.Plane;
      unsigned int tpc   = wireID.TPC;
      unsigned int nWiresTPC = 400;
      unsigned int wireGap = 4;
      double driftLen = geom->TPC(geo::TPCID{0,tpc}).DriftDistance();
      double apaLen = geom->TPC(geo::TPCID{0,tpc}).Width() - geom->TPC(geo::TPCID{0,tpc}).ActiveWidth();
      //double driftVel = detprop->DriftVelocity();
    //  unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC
    //  unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC
      double drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC
      double apa_size = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC

       //std::cout << "MEASUREMENTS: driftLen, apaLen, drift_size, apa_size " << driftLen<<" "<<apaLen<<" "<<drift_size << ", " << apa_size << std::endl;

      globalWire = 0;
      globalPlane = 0;

      // Collection plane has more wires
      if(plane == 2){
        nWiresTPC=480;
        wireGap = 5;
        globalPlane = 2;
      }

      bool includeZGap = true;
      if(includeZGap) nWiresTPC += wireGap;

      // Workspace geometry has two drift regions
      //                  |-----|-----| /  /
      //      y ^         |  3  |  2  |/  /
      //        | -| z    |-----|-----|  /
      //        | /       |  1  |  0  | /
      //  x <---|/        |-----|-----|/
      //

      int tpcMod4 = tpc%4;
      int offset = 0;
      // Induction views depend on the drift direction
      if(plane < 2){
        // For TPCs 0 and 3 keep U and V as defined.
        if(tpcMod4 == 0 || tpcMod4 == 3){
          globalPlane = plane;
          // But reverse the TDCs
        }
        // For TPCs 1 and 2, swap U and V.
        else{
          if(plane == 0) globalPlane = 1;
          else globalPlane = 0;
        }
      }
      if(globalPlane != 1){
        globalWire += (tpc/4)*nWiresTPC + (tpcMod4>1)*offset + localWire;
      }
      else{
        globalWire += ((23-tpc)/4)*nWiresTPC + (tpcMod4>1)*offset + localWire;
      }

      if(tpcMod4 == 0 || tpcMod4 == 2){
        globalTDC = drift_size - localTDC;
      }
      else{
        globalTDC = localTDC + drift_size + apa_size;
      }
  }
}


double recoestudy::RecoEnergyS::GetPixelMapWireIdx(const geo::WireID& wireID){

   double globalWire = -9999;
   geo::PlaneID planeID{wireID.Cryostat, 0, wireID.Plane};
   unsigned int nwires = wireReadout->Nwires(planeID);

   // Induction
   if (wireReadout->SignalType(wireID) == geo::kInduction) {
     //auto WireCentre = m_geom->WireIDToWireGeo(wireID).GetCenter();
     auto WireCentre = wireReadout->Wire(wireID).GetCenter();

     int temp_tpc = 0;
     if (wireID.TPC % 2 == 0) {
		temp_tpc = 0;
     }
     else {
		temp_tpc = 1;
     }
     geo::PlaneID const p1(wireID.Cryostat, temp_tpc, wireID.Plane);
     globalWire = wireReadout->Plane(p1).WireCoordinate(WireCentre);
   }

   // Collection
   else {
       int block = wireID.TPC / 4;
       globalWire = (double)( ((int)nwires*block) + wireID.Wire );
   }

    return round(globalWire);

}

int recoestudy::RecoEnergyS::FindTrackID(art::Ptr<recob::Hit> const& hit, detinfo::DetectorClocksData const& clockData) {
  double particleEnergy = 0;
  int likelyTrackID = 0;
  //std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackID(hit);
  std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }
  return likelyTrackID;
}

void recoestudy::RecoEnergyS::reset() {
  trueEnergy = 0;
  trueVertex_X = -99999;
  trueVertex_Y = -99999;
  trueVertex_Z = -99999;
  trueEnd_X = 0;
  trueEnd_Y = 0;
  trueEnd_Z = 0;
  truePx = 0;
  truePy = 0;
  truePz = 0;

  ErecoNu = -1;
  RecoLepEnNu   = -1;
  RecoHadEnNu   = -1;
  RecoMethodNu  = -1;

  for (int ii = 0; ii < 3; ++ii){
  	trueVtxWire[ii] = -999;
  	trueVtxGlobalWire[ii] = -999;
  	trueVtxGlobalTick[ii] = -999;
  	trueVtxGlobalPlane[ii] = 0;
  	trueVtxTPC[ii]  = -999;
  	trueVtxTime[ii] = -999;
  	trueVtxChan[ii] = -999;

  	trueMidWire[ii] = -999;
  	trueMidTime[ii] = -999;

  	trueEndWire[ii] = -999;
  	trueEndTime[ii] = -999;
    trueEndGlobalWire[ii] = -999;
  	trueEndGlobalTime[ii] = -999;

  	recoVtxWire[ii] = -999;
  	recoVtxTime[ii] = -999;
  	shwVtxTPC[ii] = -999;
  	shwVtxWire[ii] = -999;
  	shwVtxTime[ii] = -999;
  	shwVtxGlobalWire[ii] = -999;
  	shwVtxGlobalTick[ii] = -999;

	PFVtxTPC[ii] = -9999;
	PFVtxWire[ii] = -9999;
	PFVtxTime[ii] = -9999;
	PFVtxGlobalWire[ii] = -9999;
	PFVtxGlobalTick[ii] = -9999;
	PFVtxGlobalPlane[ii] = 0;
	hit_offset[ii] = 0;

        CVN_vertex1ststep[ii] = -9999;
        CVN_vertex2ndstep[ii] = -9999;

        CVN_vertex2[ii] = -9999;
  	CVNVtx2Wire[ii] = -9999;
  	CVNVtx2TPC[ii]  = -9999;
  	CVNVtx2Time[ii] = -9999;
  	CVNVtx2GlobalWire[ii] = -9999;
  	CVNVtx2GlobalTick[ii] = -9999;
  	CVNVtx2GlobalPlane[ii] = 0;

        CVN_vertex[ii] = -9999;
  	CVNVtxWire[ii] = -9999;
  	CVNVtxTPC[ii]  = -9999;
  	CVNVtxTime[ii] = -9999;
  	CVNVtxGlobalWire[ii] = -9999;
  	CVNVtxGlobalTick[ii] = -9999;
  	CVNVtxGlobalPlane[ii] = 0;

  	EndProngGlobalWire[ii] = -9999;
  	EndProngGlobalTick[ii] = -9999;
  	EndProngGlobalPlane[ii] = 0;

  }

  evthaspripf = 0;
  recoVertex_X = -9999;
  recoVertex_Y = -9999;
  recoVertex_Z = -9999;

  depositU = 0;
  depositV = 0;
  depositZ = 0;
  vertexDetectorDist = 0;
  nhits = 0;
  shw_tot_chg = 0;
  led_shw_idx = -999;
  led_shw_id = -999;
  led_shw_eng = -999;
  led_shw_dirx = -999;
  led_shw_diry = -999;
  led_shw_dirz = -999;
  led_shw_startx = -999;
  led_shw_starty = -999;
  led_shw_startz = -999;
  led_shw_costh = -999;

  shw_first_wire = -999;
  shw_first_time = -999;
  TrueParticle_size = 0;

  led_shw_trueP0 = 0;
  led_shw_trueP1 = 0;
  led_shw_trueP2 = 0;
  led_shw_trueP3 = 0;
  led_shw_truedang = 0;
  led_shw_truePDG = 0;
  led_shw_truePDGMom = 0;
  led_shw_truecosth = 0;


  for (int hit = 0; hit < kMaxHits; ++hit) {
    hit_cryst[hit] = 0;
    hit_tpc[hit] = 0;
    hit_plane[hit] = 0;
    hit_wire[hit] = 0;
    hit_channel[hit] = 0;
    hit_multi[hit] = 0;
    hit_peakT[hit] = 0;
    hit_amp[hit] = 0;
    hit_minusT[hit] = 0;
    hit_plusT[hit] = 0;
    hit_rmsT[hit] = 0;
    hit_startT[hit] = 0;
    hit_endT[hit] = 0;
    hit_charge[hit] = -999;
    hit_corrcharge[hit] = -999;
    hit_energy[hit] = -999;
    hit_global_wire[hit] = -999;
    hit_global_tick[hit] = -999;
    hit_global_plane[hit] = -999;
    hit_good[hit] = -999;
    hit_ndf[hit] = -999;

    hit_prong_tag[hit] = -1;
    hit_prong_pdg[hit] = -1;
    hit_prong_pdg_mom[hit] = -1;

    hit_shwid[hit] = -999;
    hit_ledshwid[hit] = 0;
    hit_shwdirX[hit] = -999;
    hit_shwdirY[hit] = -999;
    hit_shwdirZ[hit] = -999;
    hit_shwcosth[hit] = -999;
    hit_shwphi[hit] = -999;
    hit_shwstartT[hit] = -999;
    hit_shwendT[hit] = -999;
    hit_shwevteng[hit] = -999;
    hit_shwevtdedx[hit] = -999;

    hit_ntrackID[hit] = 0;
    hit_trueE[hit] = -999;

    hit_nID_Pri[hit] = 0;
    hit_PDG[hit] = 0;
    hit_Mo_PDG[hit] = 0;
    hit_mu_flag[hit] = 0;
    hit_ele_flag[hit] = 0;

    hit_x[hit] = -999;
    hit_y[hit] = -999;
    hit_z[hit] = -999;

    hit_sp_x[hit] = -999;
    hit_sp_y[hit] = -999;
    hit_sp_z[hit] = -999;

    hit_Nallraws[hit] = -999;
    hit_Nadcvec[hit] = -999;
    hit_raw_wireID[hit] = -999;
    hit_raw_chID[hit] = -999;

    wire_roifirst[hit] = -999;
    wire_roiend[hit] = -999;
    sum_ch_chg[hit] = -999;

    shw_idx[hit] = -999;
    hit_shw_tpc[hit] = -999;
    hit_shw_plane[hit] = -999;
    hit_shw_wire[hit] = -999;
    hit_shw_channel[hit] = -999;
    hit_shw_multi[hit] = -999;
    hit_shw_truetrackid[hit] = -999;
    hit_shw_ntrackID[hit] = -999;
    hit_shw_peakT[hit] = -999;
    hit_shw_startT[hit] = -999;
    hit_shw_endT[hit] = -999;
    hit_shw_charge[hit] = -999;
    hit_shw_corrcharge[hit] = -999;
    hit_shw_energy[hit] = -999;
    hit_shw_trueE[hit] = -999;
    hit_shw_pdg[hit] = -999;
    hit_shw_pdg_mom[hit] = -999;
    hit_shw_global_plane[hit] = -999;
    hit_shw_global_wire[hit] = -999;
    hit_shw_global_tick[hit] = -999;
    shw_startX[hit] = -999;
    shw_startY[hit] = -999;
    shw_startZ[hit] = -999;
    shw_dirX[hit] = -999;
    shw_dirY[hit] = -999;
    shw_dirZ[hit] = -999;
    wire_assn[hit] = -999;
    wire_min_tick[hit] = -999;
    wire_max_tick[hit] = -999;
    wire_max_rmsT[hit] = -999;

    wire_min_ledshwtick[hit] = -999;
    wire_max_ledshwtick[hit] = -999;

  }
  for (int ii = 0; ii < kMaxClst; ii++){
    max_shwtick[ii] = -999;
    min_shwtick[ii] = 1e+5;
    max_trktick[ii] = -999;
    min_trktick[ii] = 1e+5;
    max_shwrms[ii] = -999;
    max_trkrms[ii] = -999;
    for (int hit = 0; hit < kMaxHits; hit++){
        wire_max_shwtick[hit][ii] = -999;
        wire_min_shwtick[hit][ii] = 1e+5;
        wire_max_trktick[hit][ii] = -999;
        wire_min_trktick[hit][ii] = 1e+5;
        wire_max_shwrms[hit][ii] = -999;
        wire_max_trkrms[hit][ii] = -999;
    }
  }


  for (int ii = 0; ii < kMax; ++ii) {
    all_track_mom[ii] = -999;
    all_track_calmom[ii] = -999;
    all_track_length[ii] = -999;
    all_track_startx[ii] = -9999;
    all_track_starty[ii] = -9999;
    all_track_startz[ii] = -9999;
    all_track_px[ii] = -9999;
    all_track_py[ii] = -9999;
    all_track_pz[ii] = -9999;
    all_track_true_Efrac[ii] = -999;
    all_track_true_pdg[ii] = -999;
    all_track_true_pdg_mom[ii] = -999;
    all_track_dist_reco_vtx[ii] = -999;
    all_track_dist_true_vtx[ii] = -999;
    all_track_min_angle[ii] = -999;
    all_track_true_eng[ii] = -999;
    all_track_true_px[ii] = -999;
    all_track_true_py[ii] = -999;
    all_track_true_pz[ii] = -999;
    all_track_true_mom_px[ii] = -999;
    all_track_true_mom_py[ii] = -999;
    all_track_true_mom_pz[ii] = -999;

    all_track_endx[ii] = -9999;
    all_track_endy[ii] = -9999;
    all_track_endz[ii] = -9999;

    all_track_true_endx[ii] = -9999;
    all_track_true_endy[ii] = -9999;
    all_track_true_endz[ii] = -9999;

    all_track2_true_Efrac[ii] = -999;
    all_track2_true_pdg[ii] = -999;
    all_track2_true_pdg_mom[ii] = -999;
    all_track2_true_eng[ii] = -999;
    all_track2_true_px[ii] = -999;
    all_track2_true_py[ii] = -999;
    all_track2_true_pz[ii] = -999;
    all_track2_true_mom_px[ii] = -999;
    all_track2_true_mom_py[ii] = -999;
    all_track2_true_mom_pz[ii] = -999;


    shw_chg[ii] = 0;
    shw_nhits[ii] = 0;
    all_shw_idx[ii] = -1;
    all_shw_eng[ii] = -999;
    all_shw_dedx[ii] = -999;
    all_shw_length[ii] = -999;
    all_shw_startx[ii] = -9999;
    all_shw_starty[ii] = -9999;
    all_shw_startz[ii] = -9999;
    all_shw_px[ii] = -9999;
    all_shw_py[ii] = -9999;
    all_shw_pz[ii] = -9999;
    all_shw_true_eng[ii] = -999;
    all_shw_true_px[ii] = -999;
    all_shw_true_py[ii] = -999;
    all_shw_true_pz[ii] = -999;
    all_shw_true_mom_px[ii] = -999;
    all_shw_true_mom_py[ii] = -999;
    all_shw_true_mom_pz[ii] = -999;

    all_shw2_true_Efrac[ii] = -999;
    all_shw2_true_pdg[ii] = -999;
    all_shw2_true_pdg_mom[ii] = -999;
    all_shw2_true_eng[ii] = -999;
    all_shw2_true_px[ii] = -999;
    all_shw2_true_py[ii] = -999;
    all_shw2_true_pz[ii] = -999;
    all_shw2_true_mom_px[ii] = -999;
    all_shw2_true_mom_py[ii] = -999;
    all_shw2_true_mom_pz[ii] = -999;



    all_shw_costh[ii] = -999;
    all_shw_phi[ii] = -999;
    all_shw_true_Efrac[ii] = -999;
    all_shw_true_pdg[ii] = -999;
    all_shw_true_pdg_mom[ii] = -999;
    all_shw_dist_reco_vtx[ii] = -999;
    all_shw_dist_true_vtx[ii] = -999;
    all_shw_min_angle[ii] = -999;
    clst_plane[ii] = -999;
    clst_integral[ii] = -999;
    clst_sumadc[ii]  = -999;
    clst_width[ii]  = -999;
    clst_nhits[ii]  = -999;
    clst_startwire[ii] = -999;
    clst_starttick[ii] = -999;
    clst_endwire[ii]  = -999;
    clst_endtick[ii]  = -999;
    nueng_truth[ii] = -999;
    nupdg_truth[ii] = -999;
    nuvtxx_truth[ii] = -999;
    nuvtxy_truth[ii] = -999;
    nuvtxz_truth[ii] = -999;
    numode_truth[ii] = -999;
    nuccnc_truth[ii] = -999;
    leppdg_truth[ii] = -999;
    lepp0_truth[ii] = -999;
    lepp1_truth[ii] = -999;
    lepp2_truth[ii] = -999;
    lepeng_truth[ii] = -999;

  }
  mutrk_contain = -999;
  ChgU.clear();
  ChgV.clear();
  ChgZ.clear();

  hit_rawadc.clear();
  hit_rawtick.clear();
  wire_adc.clear();
  wire_corradc.clear();
  wire_x.clear();
  wire_tick.clear();
  wire_shwflag.clear();
  wire_ledshwflag.clear();
  wire_pdg.clear();
  wire_pdg_mom.clear();
  wire_dr.clear();
  wire_dist.clear();
  wire_tag_pdg.clear();
  wire_tag_pdg_mom.clear();
  wire_tag_is_shw.clear();
  wire_tag_ith.clear();
  wire_tag_true_px.clear();
  wire_tag_true_py.clear();
  wire_tag_true_pz.clear();
  wire_tag_reco_px.clear();
  wire_tag_reco_py.clear();
  wire_tag_reco_pz.clear();
  wire_tag_reco_eng.clear();



  ch_chg.clear();
  ch_tick.clear();

  TruePart_PDG.clear();
  TruePart_E.clear();
  trueDaughterPDGs.clear();
  trueDaughterE.clear();
  trueDaughterDepoE.clear();
  trueDaughterPx.clear();
  trueDaughterPy.clear();
  trueDaughterPz.clear();
  v3_on_planes[0]->SetXYZ(0,0,0);
  v3_on_planes[1]->SetXYZ(0,0,0);
  v3_on_planes[2]->SetXYZ(0,0,0);

}

DEFINE_ART_MODULE(recoestudy::RecoEnergyS)
