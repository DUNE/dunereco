////////////////////////////////////////////////////////////////////////
// Class:       CCNuSelection
// Plugin Type: analyzer (art v2_07_03)
// File:        CCNuSelection_module.cc
//
// Generated at Tue Jun 27 06:07:56 2017 by Dominic Brailsford using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//STL
#include <iostream>
//ROOT
#include "TTree.h"
#include "TMVA/Reader.h"

//ART
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


//LArSoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
//#include "larreco/RecoAlg/ShowerEnergyAlg.h"

#include "larsim/Utils/TruthMatchUtils.h"

//DUNE
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/TrackPID/algorithms/CTPHelper.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/TrackPID/products/CTPResult.h"

#include "dunereco/CVN/func/Result.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

//Custom
#include "dunereco/FDSelections/pandizzle/PandizzleAlg.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"
#include "FDSelectionUtils.h"
#include "tools/RecoTrackSelector.h"
#include "tools/RecoShowerSelector.h"


constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

constexpr int kDefMaxNTrueVertexParticles = 500;
constexpr int kDefMaxNRecoTracks = 1000;
constexpr int kDefMaxNRecoShowers = 1000;

namespace FDSelection {
  class CCNuSelection;
}


class FDSelection::CCNuSelection : public art::EDAnalyzer {
public:
  explicit CCNuSelection(fhicl::ParameterSet const & p);
  CCNuSelection(CCNuSelection const &) = delete;
  CCNuSelection(CCNuSelection &&) = delete;
  CCNuSelection & operator = (CCNuSelection const &) = delete;
  CCNuSelection & operator = (CCNuSelection &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void beginSubRun(art::SubRun const & sr) override;
  void endSubRun(art::SubRun const & sr) override;
  void endJob() override;

private:
  void Reset();      //Resets all tree vars
  void GetEventInfo(art::Event const & evt); //Grab event-level info that is necessary/handy for selecting events but doesn't really fall under the 'selection' banner
  void GetTruthInfo(art::Event const & evt);  //Grab the truth info from the art record
  void FillVertexInformation(art::Event const & evt); //Grab the vertex parameters 
  void GetRecoTrackInfo(art::Event const & evt); //Fill the reco track arrays
  void RunTrackSelection(art::Event const & evt);  //Run the selection and dump relevant info to the truee
  bool IsTrackContained(art::Ptr<recob::Track> const track, std::vector< art::Ptr<recob::Hit > > const track_hits, art::Event const & evt); // check if the track is contained in the detector
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt);
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToShower(art::Ptr<recob::Shower> const shower, art::Event const & evt);
  art::Ptr<recob::Track> GetTrackMatchedPFParticle(art::Ptr<recob::PFParticle> const pfparticle, art::Event const & evt);

  void FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp);
  void FillChildPFPInformation(art::Ptr<recob::Track> const track, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp);
  void FillChildPFPInformation(art::Ptr<recob::Shower> const shower, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp);


  //Shower stuff
  void GetRecoShowerInfo(art::Event const & evt); //Fill the reco shower arrays
  void RunShowerSelection(art::Event const & evt);  //Run the selection and dump relevant info to the truee
  //double CalculateShowerCharge(art::Ptr<recob::Shower> const track, std::vector< art::Ptr< recob::Hit> > const track_hits); //Calculate hit charge on track as measured by the collection plane
  bool IsPFParticlePrimaryDaughter(art::Ptr<recob::PFParticle> const & pfparticle, std::vector<art::Ptr<recob::PFParticle > > const & pfparticles);

  TVector3 ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector); //Returns a projection o a vector onto a plane

  double DistanceToNuVertex(art::Event const & evt, std::vector< art::Ptr<recob::Hit> > const artHitList, const TVector3 &nuVertex);
  double YZtoU(const geo::WireID &hit_WireID, double y, double x);
  double YZtoV(const geo::WireID &hit_WireID, double y, double x);
  double YZtoW(const geo::WireID &hit_WireID, double y, double x);

  // Declare member data here.


  TTree *fTree; //The selection tree
  //Generic stuff
  int fRun;
  int fSubRun;
  int fEvent;
  int fIsMC;
  //Detector stuff
  double fT0;
  //Neutrino stuff
  int fNuPdg; //Interaction PDG
  int fBeamPdg; //PDG at point of creation
  int fNC;    // 1=is NC, 0=otherwise
  int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  int fTargetZ; //Atomic number of scattering target
  double fQ2; 
  double fENu; 
  double fW; //X-Sec params
  double fX;
  double fY;
  double fNuMomX; //Neutrino momentums
  double fNuMomY;
  double fNuMomZ;
  double fNuMomT;
  double fNuX; //Interaction positions
  double fNuY;
  double fNuZ;
  double fNuT;
  //Outgoing particle count
  int fNPiP; //Number of pi+
  int fNPim; //Number of pi-
  int fNPi0; //Number of pi0
  int fNPhotons; //Number of photons
  int fNProtons; //Number of protons
  int fNNeutrons; //Number of neutrinos
  int fNOther; //Number of other particles
  int fNVertexParticles; //The total number of particles attached to the GHEP vertex
  bool fVertexParticleIsGHEP[kDefMaxNTrueVertexParticles];
  int fVertexParticlePDG[kDefMaxNTrueVertexParticles]; //The PDG of each particle in the GHEP vetex
  int fVertexParticleStatus[kDefMaxNTrueVertexParticles];
  int fVertexParticleNChildren[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomX[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomY[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomZ[kDefMaxNTrueVertexParticles];
  double fVertexParticleMomT[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndX[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndY[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndZ[kDefMaxNTrueVertexParticles];
  double fVertexParticleEndT[kDefMaxNTrueVertexParticles];
  //Outgoing lepton stuff
  int fLepPDG;
  double fMomLepX;
  double fMomLepY;
  double fMomLepZ;
  double fMomLepT;
  double fLepEndX;
  double fLepEndY;
  double fLepEndZ;
  double fLepEndT;
  double fLepNuAngle;
  //Transverse mom stuff
  double  fNuMomTranMag;
  double  fTargNuclMomTranMag;
  double  fInitalMomTranMag;
  double  fLepMomTranMag;
  double  fNuclRemMomTranMag;
  double  fFinalMomTranMagNoLepNoRem;
  double  fFinalMomTranMagNoLepWithRem;
  double  fFinalMomTranMagWithLepNoRem;
  double  fFinalMomTranMagWithLepWithRem;
  //Selection stuff
  //true bits
  int fSelTrackTruePDG;
  bool fSelTrackTruePrimary;
  double fSelTrackTrueMomX;
  double fSelTrackTrueMomY;
  double fSelTrackTrueMomZ;
  double fSelTrackTrueMomT;
  double fSelTrackTrueStartX;
  double fSelTrackTrueStartY;
  double fSelTrackTrueStartZ;
  double fSelTrackTrueStartT;
  double fSelTrackTrueEndX;
  double fSelTrackTrueEndY;
  double fSelTrackTrueEndZ;
  double fSelTrackTrueEndT;
  //reco bits
  int fSelTrackRecoNHits;
  double fSelTrackRecoCompleteness;
  double fSelTrackRecoHitPurity;
  double fSelTrackRecoMomX;
  double fSelTrackRecoMomY;
  double fSelTrackRecoMomZ;
  double fSelTrackRecoMomT;
  double fSelTrackRecoStartX;
  double fSelTrackRecoStartY;
  double fSelTrackRecoStartZ;
  double fSelTrackRecoStartT;
  double fSelTrackRecoEndX;
  double fSelTrackRecoEndY;
  double fSelTrackRecoEndZ;
  double fSelTrackRecoEndT;
  double fSelTrackRecoUpstreamX;
  double fSelTrackRecoUpstreamY;
  double fSelTrackRecoUpstreamZ;
  double fSelTrackRecoUpstreamT;
  double fSelTrackRecoDownstreamX;
  double fSelTrackRecoDownstreamY;
  double fSelTrackRecoDownstreamZ;
  double fSelTrackRecoDownstreamT;
  double fSelTrackRecoEndClosestToVertexX;
  double fSelTrackRecoEndClosestToVertexY;
  double fSelTrackRecoEndClosestToVertexZ;
  double fSelTrackRecoLength;
  int   fSelTrackRecoContained;
  int   fSelTrackRecoMomMethod;
  double fSelTrackRecoCharge;
  double fSelTrackRecoMomMCS;
  double fSelTrackRecoMomContained;
  //int fSelTrackRecoNMatchedVertices;
  double fSelTrackRecoVertexX;
  double fSelTrackRecoVertexY;
  double fSelTrackRecoVertexZ;
  int fSelTrackRecoNChildPFP;
  int fSelTrackRecoNChildTrackPFP;
  int fSelTrackRecoNChildShowerPFP;
  //MVA bits
  double fSelTrackMVAElectron;
  double fSelTrackMVAPion;
  double fSelTrackMVAMuon;
  double fSelTrackMVAProton;
  double fSelTrackMVAPhoton;
  //Pandizzle var
  double fSelTrackMichelNHits;
  double fSelTrackMichelElectronMVA;
  double fSelTrackMichelRecoEnergyPlane2;
  double fSelTrackDeflecAngleSD;
  double fSelTrackLength;
  double fSelTrackEvalRatio;
  double fSelTrackConcentration;
  double fSelTrackCoreHaloRatio;
  double fSelTrackConicalness;
  double fSelTrackdEdxStart;
  double fSelTrackdEdxEnd;
  double fSelTrackdEdxEndRatio;
  double fSelTrackPandizzleVar;
  //DeepPan vars
  double fSelTrackDeepPanMuVar;
  double fSelTrackDeepPanPiVar;
  double fSelTrackDeepPanProtonVar;

  //All tracks
  int fNRecoTracks;
  bool fRecoTrackIsPrimary[kDefMaxNRecoTracks];
  int fRecoTrackTruePDG[kDefMaxNRecoTracks];
  bool fRecoTrackTruePrimary[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomX[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomY[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomZ[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomT[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartX[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartY[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartZ[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartT[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndX[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndY[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndZ[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndT[kDefMaxNRecoTracks];
  //reco bits
  int fRecoTrackRecoNHits[kDefMaxNRecoTracks];
  double fRecoTrackRecoCompleteness[kDefMaxNRecoTracks];
  double fRecoTrackRecoHitPurity[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomX[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomY[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomT[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartX[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartY[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartT[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndX[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndY[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndT[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamX[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamY[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamT[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamX[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamY[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamT[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndClosestToVertexX[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndClosestToVertexY[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndClosestToVertexZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoLength[kDefMaxNRecoTracks];
  bool   fRecoTrackRecoContained[kDefMaxNRecoTracks];
  double fRecoTrackRecoCharge[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomMCS[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomContained[kDefMaxNRecoTracks];
  //int fRecoTrackRecoNMatchedVertices[kDefMaxNRecoTracks];
  double fRecoTrackRecoVertexX[kDefMaxNRecoTracks];
  double fRecoTrackRecoVertexY[kDefMaxNRecoTracks];
  double fRecoTrackRecoVertexZ[kDefMaxNRecoTracks];
  int fRecoTrackRecoNChildPFP[kDefMaxNRecoTracks];
  int fRecoTrackRecoNChildTrackPFP[kDefMaxNRecoTracks];
  int fRecoTrackRecoNChildShowerPFP[kDefMaxNRecoTracks];
  //MVA bits
  double fRecoTrackMVAElectron[kDefMaxNRecoTracks];
  double fRecoTrackMVAPion[kDefMaxNRecoTracks];
  double fRecoTrackMVAMuon[kDefMaxNRecoTracks];
  double fRecoTrackMVAProton[kDefMaxNRecoTracks];
  double fRecoTrackMVAPhoton[kDefMaxNRecoTracks];
  //DeepPan
  double fRecoTrackDeepPanMuVar[kDefMaxNRecoTracks];
  double fRecoTrackDeepPanPiVar[kDefMaxNRecoTracks];
  double fRecoTrackDeepPanProtonVar[kDefMaxNRecoTracks];
  //Pandizzle
  float fRecoTrackMichelNHits[kDefMaxNRecoTracks];
  float fRecoTrackMichelElectronMVA[kDefMaxNRecoTracks];
  float fRecoTrackMichelRecoEnergyPlane2[kDefMaxNRecoTracks];
  float fRecoTrackDeflecAngleSD[kDefMaxNRecoTracks];
  float fRecoTrackLength[kDefMaxNRecoTracks];
  float fRecoTrackEvalRatio[kDefMaxNRecoTracks];
  float fRecoTrackConcentration[kDefMaxNRecoTracks];
  float fRecoTrackCoreHaloRatio[kDefMaxNRecoTracks];
  float fRecoTrackConicalness[kDefMaxNRecoTracks];
  float fRecoTrackdEdxStart[kDefMaxNRecoTracks];
  float fRecoTrackdEdxEnd[kDefMaxNRecoTracks];
  float fRecoTrackdEdxEndRatio[kDefMaxNRecoTracks];
  double fRecoTrackPandizzleVar[kDefMaxNRecoTracks];

  //Event-level bits
  double fRecoEventCharge; //Total collected charge (as measured by the collection planes)
  //reco energy bits
  double fNumuRecoMomLep; //Reco lepton momentum
  double fNumuRecoEHad; //Reco hadronic energy
  double fNumuRecoENu; //Reco neutrino energy

  //Reco neutrino bits
  double fRecoNuVtxX;
  double fRecoNuVtxY;
  double fRecoNuVtxZ;
  int fRecoNuVtxNShowers;
  int fRecoNuVtxNTracks;
  int fRecoNuVtxNChildren;

  //Shower Selection stuff
  //true bits
  int fSelShowerTruePDG;
  bool fSelShowerTruePrimary;
  double fSelShowerTrueMomX;
  double fSelShowerTrueMomY;
  double fSelShowerTrueMomZ;
  double fSelShowerTrueMomT;
  double fSelShowerTrueStartX;
  double fSelShowerTrueStartY;
  double fSelShowerTrueStartZ;
  double fSelShowerTrueStartT;
  double fSelShowerTrueEndX;
  double fSelShowerTrueEndY;
  double fSelShowerTrueEndZ;
  double fSelShowerTrueEndT;
  //reco bits
  int fSelShowerRecoNHits;
  double fSelShowerRecoCompleteness;
  double fSelShowerRecoHitPurity;
  double fSelShowerRecoMomX;
  double fSelShowerRecoMomY;
  double fSelShowerRecoMomZ;
  double fSelShowerRecoMomT;
  double fSelShowerRecoStartX;
  double fSelShowerRecoStartY;
  double fSelShowerRecoStartZ;
  double fSelShowerRecoStartT;
  double fSelShowerRecoEndX;
  double fSelShowerRecoEndY;
  double fSelShowerRecoEndZ;
  double fSelShowerRecoEndT;
  double fSelShowerRecoUpstreamX;
  double fSelShowerRecoUpstreamY;
  double fSelShowerRecoUpstreamZ;
  double fSelShowerRecoUpstreamT;
  double fSelShowerRecoDownstreamX;
  double fSelShowerRecoDownstreamY;
  double fSelShowerRecoDownstreamZ;
  double fSelShowerRecoDownstreamT;
  double fSelShowerRecoEndClosestToVertexX;
  double fSelShowerRecoEndClosestToVertexY;
  double fSelShowerRecoEndClosestToVertexZ;
  bool   fSelShowerRecoContained;
  double fSelShowerRecoCharge;
  double fSelShowerRecoMomMCS;
  double fSelShowerRecoMomContained;
  int fSelShowerRecoNMatchedVertices;
  double fSelShowerRecoVertexX;
  double fSelShowerRecoVertexY;
  double fSelShowerRecoVertexZ;
  int fSelShowerRecoNChildPFP;
  int fSelShowerRecoNChildTrackPFP;
  int fSelShowerRecoNChildShowerPFP;
  double fSelShowerRecoDirX;
  double fSelShowerRecoDirY;
  double fSelShowerRecoDirZ;
  double fSelShowerRecodEdx[3];
  int fSelShowerRecoBestPlane;
  double fSelShowerRecoLength;
  double fSelShowerRecoOpeningAngle;
  bool fSelShowerRecoIsPrimaryPFPDaughter;
  //MVA bits
  double fSelShowerMVAElectron;
  double fSelShowerMVAPion;
  double fSelShowerMVAMuon;
  double fSelShowerMVAProton;
  double fSelShowerMVAPhoton;
  double fSelShowerRecoEnergy[3];
  //Pandrizzle
  double fSelShowerPandrizzleConnectionBDTScore;
  double fSelShowerPandrizzlePathwayLengthMin;
  double fSelShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D;
  double fSelShowerPandrizzleMaxNPostShowerStartHits;
  double fSelShowerPandrizzleMaxPostShowerStartScatterAngle;
  double fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry;
  double fSelShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry;
  double fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
  double fSelShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius;
  double fSelShowerPandrizzleMaxPostShowerStartOpeningAngle;
  double fSelShowerPandrizzleMaxFoundHitRatio;
  double fSelShowerPandrizzleMaxInitialGapSize;
  double fSelShowerPandrizzleMinLargestProjectedGapSize;
  double fSelShowerPandrizzleNViewsWithAmbiguousHits;
  double fSelShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy;
  double fSelShowerPandrizzleEvalRatio;
  double fSelShowerPandrizzleConcentration;
  double fSelShowerPandrizzleCoreHaloRatio;
  double fSelShowerPandrizzleConicalness;
  double fSelShowerPandrizzledEdxBestPlane;
  double fSelShowerPandrizzleDisplacement;
  double fSelShowerPandrizzleDCA;
  double fSelShowerPandrizzleWideness;
  double fSelShowerPandrizzleEnergyDensity;
  double fSelShowerPandrizzleBDTMethod;
  double fSelShowerJamPandrizzleMVAScore;
  double fSelShowerPandrizzleMVAScore;
  bool   fSelShowerPandrizzleIsFilled;
  double fSelShowerDistanceToNuVertexU;
  double fSelShowerDistanceToNuVertexV;
  double fSelShowerDistanceToNuVertexW;

  //All reco showers
  int fNRecoShowers;
  int fRecoShowerTruePDG[kDefMaxNRecoShowers];
  bool fRecoShowerTruePrimary[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomX[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomY[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomZ[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomT[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartX[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartY[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartZ[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartT[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndX[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndY[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndZ[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndT[kDefMaxNRecoShowers];
  int fRecoShowerRecoNHits[kDefMaxNRecoShowers];
  int fRecoShowerActiveCheatElectronAlg[kDefMaxNRecoShowers];
  int fRecoShowerActiveHybridElectronAlg[kDefMaxNRecoShowers];
  int fRecoShowerActiveCheatGammaAlg[kDefMaxNRecoShowers];
  int fRecoShowerActiveHybridGammaAlg[kDefMaxNRecoShowers];
  double fRecoShowerRecoCompleteness[kDefMaxNRecoShowers];
  double fRecoShowerRecoHitPurity[kDefMaxNRecoShowers];
  double fRecoShowerRecoMomX[kDefMaxNRecoShowers];
  double fRecoShowerRecoMomY[kDefMaxNRecoShowers];
  double fRecoShowerRecoMomZ[kDefMaxNRecoShowers];
  double fRecoShowerRecoMomT[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartX[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartY[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartZ[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartT[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndX[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndY[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndZ[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndT[kDefMaxNRecoShowers];
  double fRecoShowerRecoUpstreamX[kDefMaxNRecoShowers];
  double fRecoShowerRecoUpstreamY[kDefMaxNRecoShowers];
  double fRecoShowerRecoUpstreamZ[kDefMaxNRecoShowers];
  double fRecoShowerRecoUpstreamT[kDefMaxNRecoShowers];
  double fRecoShowerRecoDownstreamX[kDefMaxNRecoShowers];
  double fRecoShowerRecoDownstreamY[kDefMaxNRecoShowers];
  double fRecoShowerRecoDownstreamZ[kDefMaxNRecoShowers];
  double fRecoShowerRecoDownstreamT[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndClosestToVertexX[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndClosestToVertexY[kDefMaxNRecoShowers];
  double fRecoShowerRecoEndClosestToVertexZ[kDefMaxNRecoShowers];
  bool   fRecoShowerRecoContained[kDefMaxNRecoShowers];
  double fRecoShowerRecoCharge[kDefMaxNRecoShowers];
  double fRecoShowerRecoMomMCS[kDefMaxNRecoShowers];
  double fRecoShowerRecoMomContained[kDefMaxNRecoShowers];
  int fRecoShowerRecoNMatchedVertices[kDefMaxNRecoShowers];
  double fRecoShowerRecoVertexX[kDefMaxNRecoShowers];
  double fRecoShowerRecoVertexY[kDefMaxNRecoShowers];
  double fRecoShowerRecoVertexZ[kDefMaxNRecoShowers];
  int fRecoShowerRecoNChildPFP[kDefMaxNRecoShowers];
  int fRecoShowerRecoNChildTrackPFP[kDefMaxNRecoShowers];
  int fRecoShowerRecoNChildShowerPFP[kDefMaxNRecoShowers];
  double fRecoShowerRecoDirX[kDefMaxNRecoShowers];
  double fRecoShowerRecoDirY[kDefMaxNRecoShowers];
  double fRecoShowerRecoDirZ[kDefMaxNRecoShowers];
  double fRecoShowerRecodEdx[kDefMaxNRecoShowers][3];
  int fRecoShowerRecoBestPlane[kDefMaxNRecoShowers];
  double fRecoShowerRecoLength[kDefMaxNRecoShowers];
  double fRecoShowerRecoOpeningAngle[kDefMaxNRecoShowers];
  bool fRecoShowerRecoIsPrimaryPFPDaughter[kDefMaxNRecoShowers];
  double fRecoShowerMVAElectron[kDefMaxNRecoShowers];
  double fRecoShowerMVAPion[kDefMaxNRecoShowers];
  double fRecoShowerMVAMuon[kDefMaxNRecoShowers];
  double fRecoShowerMVAProton[kDefMaxNRecoShowers];
  double fRecoShowerMVAPhoton[kDefMaxNRecoShowers];
  double fRecoShowerRecoEnergy[kDefMaxNRecoShowers][3];
  double fRecoShowerPandrizzleConnectionBDTScore[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzlePathwayLengthMin[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxNPostShowerStartHits[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxFoundHitRatio[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMaxInitialGapSize[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMinLargestProjectedGapSize[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleNViewsWithAmbiguousHits[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleEvalRatio[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleConcentration[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleCoreHaloRatio[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleConicalness[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzledEdxBestPlane[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleDisplacement[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleDCA[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleWideness[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleEnergyDensity[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleBDTMethod[kDefMaxNRecoShowers];
  double fRecoShowerJamPandrizzleMVAScore[kDefMaxNRecoShowers];
  double fRecoShowerPandrizzleMVAScore[kDefMaxNRecoShowers];
  bool   fRecoShowerPandrizzleIsFilled[kDefMaxNRecoShowers];
  double fRecoShowerDistanceToNuVertexU[kDefMaxNRecoShowers];
  double fRecoShowerDistanceToNuVertexV[kDefMaxNRecoShowers];
  double fRecoShowerDistanceToNuVertexW[kDefMaxNRecoShowers];

  //reco energy bits
  double fNueRecoMomLep; //Reco lepton momentum
  double fNueRecoEHad; //Reco hadronic energy
  double fNueRecoENu; //Reco neutrino energy



  //POT tree stuff
  TTree* fPOTTree;
  double fPOT;

  //Fhicl pset labels
  std::string fNuGenModuleLabel;
  std::string fLargeantModuleLabel;
  std::string fWireModuleLabel;
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fVertexModuleLabel;
  std::string fPIDModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fHitsModuleLabel;
  std::string fPOTModuleLabel;
  std::string fCVNModuleLabel;
  std::string fNumuEnergyRecoModuleLabel;
  std::string fNueEnergyRecoModuleLabel;
  std::string fPandizzleWeightFileName;

  //Algs
  PandizzleAlg fPandizzleAlg;
  PandrizzleAlg fPandrizzleAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  //shower::ShowerEnergyAlg fShowerEnergyAlg;
  ctp::CTPHelper fConvTrackPID;
  dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;

  //Processing flags
  bool fUsePandoraVertex;

  //Tools
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;
  std::unique_ptr<FDSelectionTools::RecoShowerSelector> fRecoShowerSelector;

  // DUNE CVN stuff for comparison
  double fCVNResultNue;
  double fCVNResultNumu;
  double fCVNResultNutau;
  double fCVNResultNC;
};

FDSelection::CCNuSelection::CCNuSelection(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)   ,
  fNVertexParticles(0),
  fNRecoTracks(0),
  fNRecoShowers(0),
  fNuGenModuleLabel        (pset.get< std::string >("ModuleLabels.NuGenModuleLabel")),
  fLargeantModuleLabel     (pset.get< std::string >("ModuleLabels.LargeantModuleLabel")),
  fWireModuleLabel        (pset.get< std::string >("ModuleLabels.WireModuleLabel")),
  fTrackModuleLabel        (pset.get< std::string >("ModuleLabels.TrackModuleLabel")),
  fShowerModuleLabel        (pset.get< std::string >("ModuleLabels.ShowerModuleLabel")),
  fPFParticleModuleLabel   (pset.get< std::string >("ModuleLabels.PFParticleModuleLabel")),
  fVertexModuleLabel   (pset.get< std::string >("ModuleLabels.VertexModuleLabel")),
  fPIDModuleLabel          (pset.get< std::string >("ModuleLabels.PIDModuleLabel")),
  fParticleIDModuleLabel          (pset.get< std::string >("ModuleLabels.ParticleIDModuleLabel")),
  fHitsModuleLabel         (pset.get< std::string >("ModuleLabels.HitsModuleLabel")),
  fPOTModuleLabel          (pset.get< std::string >("ModuleLabels.POTModuleLabel")),
  fCVNModuleLabel(pset.get<std::string>("ModuleLabels.CVNModuleLabel")),
  fNumuEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.NumuEnergyRecoModuleLabel")),
  fNueEnergyRecoModuleLabel   (pset.get< std::string >("ModuleLabels.NueEnergyRecoModuleLabel")),
  fPandizzleWeightFileName        (pset.get< std::string > ("PandizzleWeightFileName")),
  fPandizzleAlg(pset) ,
  fPandrizzleAlg(pset),
  fCalorimetryAlg          (pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fConvTrackPID(pset.get<fhicl::ParameterSet>("ctpHelper")),
  fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"),fTrackModuleLabel,fShowerModuleLabel,
        fHitsModuleLabel,fWireModuleLabel,fTrackModuleLabel,fShowerModuleLabel,fPFParticleModuleLabel),
  fUsePandoraVertex        (pset.get< bool >("UsePandoraVertex")),
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(pset.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))}
{
}

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
  std::cout << "reset" << std::endl;
  Reset(); //Reset at the start of the event

  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();
  fIsMC = !(evt.isRealData());

  GetEventInfo(evt);
  if (fIsMC) GetTruthInfo(evt);
  FillVertexInformation(evt);
  GetRecoTrackInfo(evt);
  RunTrackSelection(evt);
  GetRecoShowerInfo(evt);
  RunShowerSelection(evt);
  fPandizzleAlg.Run(evt);
  fPandrizzleAlg.Run(evt);

  fTree->Fill();
}

void FDSelection::CCNuSelection::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("ccnusel","CC nu selection");
    fTree->Branch("Run",&fRun);
    fTree->Branch("SubRun",&fSubRun);
    fTree->Branch("Event",&fEvent);
    fTree->Branch("IsMC",&fIsMC);
    fTree->Branch("T0",&fT0);
    fTree->Branch("NuPdg",&fNuPdg);
    fTree->Branch("BeamPdg",&fBeamPdg);
    fTree->Branch("NC",&fNC);
    fTree->Branch("Mode",&fMode);
    fTree->Branch("TargetZ",&fTargetZ);
    fTree->Branch("Q2",&fQ2);
    fTree->Branch("Enu",&fENu);
    fTree->Branch("W",&fW);
    fTree->Branch("X",&fX);
    fTree->Branch("Y",&fY);
    fTree->Branch("NuMomX",&fNuMomX);
    fTree->Branch("NuMomY",&fNuMomY);
    fTree->Branch("NuMomZ",&fNuMomZ);
    fTree->Branch("NuMomT",&fNuMomT);
    fTree->Branch("NuX",&fNuX);
    fTree->Branch("NuY",&fNuY);
    fTree->Branch("NuZ",&fNuZ);
    fTree->Branch("NuT",&fNuT);
    fTree->Branch("NPiP",&fNPiP);
    fTree->Branch("NPim",&fNPim);
    fTree->Branch("NPi0",&fNPi0);
    fTree->Branch("NPhotons",&fNPhotons);
    fTree->Branch("NProton",&fNProtons);
    fTree->Branch("NNeutrons",&fNNeutrons);
    fTree->Branch("NOther",&fNOther);
    fTree->Branch("NVertexParticles",&fNVertexParticles);
    fTree->Branch("VertexParticleIsGHEP",fVertexParticleIsGHEP,"VertexParticleIsGHEP[NVertexParticles]/O");
    fTree->Branch("VertexParticlePDG",fVertexParticlePDG,"VertexParticlePDG[NVertexParticles]/I");
    fTree->Branch("VertexParticleStatus",fVertexParticleStatus,"VertexParticleStatus[NVertexParticles]/I");
    fTree->Branch("VertexParticleNChildren",fVertexParticleNChildren,"VertexParticleNChildren[NVertexParticles]/I");
    fTree->Branch("VertexParticleMomX",fVertexParticleMomX,"VertexParticleMomX[NVertexParticles]/D");
    fTree->Branch("VertexParticleMomY",fVertexParticleMomY,"VertexParticleMomY[NVertexParticles]/D");
    fTree->Branch("VertexParticleMomZ",fVertexParticleMomZ,"VertexParticleMomZ[NVertexParticles]/D");
    fTree->Branch("VertexParticleMomT",fVertexParticleMomT,"VertexParticleMomT[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndX",fVertexParticleEndX,"VertexParticleEndX[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndY",fVertexParticleEndY,"VertexParticleEndY[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndZ",fVertexParticleEndZ,"VertexParticleEndZ[NVertexParticles]/D");
    fTree->Branch("VertexParticleEndT",fVertexParticleEndT,"VertexParticleEndT[NVertexParticles]/D");
    fTree->Branch("LepPDG",&fLepPDG);
    fTree->Branch("MomLepX",&fMomLepX);
    fTree->Branch("MomLepY",&fMomLepY);
    fTree->Branch("MomLepZ",&fMomLepZ);
    fTree->Branch("MomLepT",&fMomLepT);
    fTree->Branch("LepEndX",&fLepEndX);
    fTree->Branch("LepEndY",&fLepEndY);
    fTree->Branch("LepEndZ",&fLepEndZ);
    fTree->Branch("LepEndT",&fLepEndT);
    fTree->Branch("LepNuAngle",&fLepNuAngle);
    fTree->Branch("NuMomTranMag",&fNuMomTranMag);
    fTree->Branch("TargNuclMomTranMag",&fTargNuclMomTranMag);
    fTree->Branch("InitalMomTranMag",&fInitalMomTranMag);
    fTree->Branch("LepMomTranMag",&fLepMomTranMag);
    fTree->Branch("NuclRemMomTranMag",&fNuclRemMomTranMag);
    fTree->Branch("FinalMomTranMagNoLepNoRem",&fFinalMomTranMagNoLepNoRem);
    fTree->Branch("FinalMomTranMagNoLepWithRem",&fFinalMomTranMagNoLepWithRem);
    fTree->Branch("FinalMomTranMagWithLepNoRem",&fFinalMomTranMagWithLepNoRem);
    fTree->Branch("FinalMomTranMagWithLepWithRem",&fFinalMomTranMagWithLepWithRem);
    fTree->Branch("SelTrackTruePDG",&fSelTrackTruePDG);
    fTree->Branch("SelTrackTruePrimary",&fSelTrackTruePrimary);
    fTree->Branch("SelTrackTrueMomX",&fSelTrackTrueMomX);
    fTree->Branch("SelTrackTrueMomY",&fSelTrackTrueMomY);
    fTree->Branch("SelTrackTrueMomZ",&fSelTrackTrueMomZ);
    fTree->Branch("SelTrackTrueMomT",&fSelTrackTrueMomT);
    fTree->Branch("SelTrackTrueStartX",&fSelTrackTrueStartX);
    fTree->Branch("SelTrackTrueStartY",&fSelTrackTrueStartY);
    fTree->Branch("SelTrackTrueStartZ",&fSelTrackTrueStartZ);
    fTree->Branch("SelTrackTrueStartT",&fSelTrackTrueStartT);
    fTree->Branch("SelTrackTrueEndX",&fSelTrackTrueEndX);
    fTree->Branch("SelTrackTrueEndY",&fSelTrackTrueEndY);
    fTree->Branch("SelTrackTrueEndZ",&fSelTrackTrueEndZ);
    fTree->Branch("SelTrackTrueEndT",&fSelTrackTrueEndT);
    fTree->Branch("SelTrackRecoNHits",&fSelTrackRecoNHits);
    fTree->Branch("SelTrackRecoCompleteness",&fSelTrackRecoCompleteness);
    fTree->Branch("SelTrackRecoHitPurity",&fSelTrackRecoHitPurity);
    fTree->Branch("SelTrackRecoMomX",&fSelTrackRecoMomX);
    fTree->Branch("SelTrackRecoMomY",&fSelTrackRecoMomY);
    fTree->Branch("SelTrackRecoMomZ",&fSelTrackRecoMomZ);
    fTree->Branch("SelTrackRecoMomT",&fSelTrackRecoMomT);
    fTree->Branch("SelTrackRecoStartX",&fSelTrackRecoStartX);
    fTree->Branch("SelTrackRecoStartY",&fSelTrackRecoStartY);
    fTree->Branch("SelTrackRecoStartZ",&fSelTrackRecoStartZ);
    fTree->Branch("SelTrackRecoStartT",&fSelTrackRecoStartT);
    fTree->Branch("SelTrackRecoEndX",&fSelTrackRecoEndX);
    fTree->Branch("SelTrackRecoEndY",&fSelTrackRecoEndY);
    fTree->Branch("SelTrackRecoEndZ",&fSelTrackRecoEndZ);
    fTree->Branch("SelTrackRecoEndT",&fSelTrackRecoEndT);
    fTree->Branch("SelTrackRecoUpstreamX",&fSelTrackRecoUpstreamX);
    fTree->Branch("SelTrackRecoUpstreamY",&fSelTrackRecoUpstreamY);
    fTree->Branch("SelTrackRecoUpstreamZ",&fSelTrackRecoUpstreamZ);
    fTree->Branch("SelTrackRecoUpstreamT",&fSelTrackRecoUpstreamT);
    fTree->Branch("SelTrackRecoDownstreamX",&fSelTrackRecoDownstreamX);
    fTree->Branch("SelTrackRecoDownstreamY",&fSelTrackRecoDownstreamY);
    fTree->Branch("SelTrackRecoDownstreamZ",&fSelTrackRecoDownstreamZ);
    fTree->Branch("SelTrackRecoDownstreamT",&fSelTrackRecoDownstreamT);
    fTree->Branch("SelTrackRecoEndClosestToVertexX",&fSelTrackRecoEndClosestToVertexX);
    fTree->Branch("SelTrackRecoEndClosestToVertexY",&fSelTrackRecoEndClosestToVertexY);
    fTree->Branch("SelTrackRecoEndClosestToVertexZ",&fSelTrackRecoEndClosestToVertexZ);
    fTree->Branch("SelTrackRecoLength",&fSelTrackRecoLength);
    fTree->Branch("SelTrackRecoContained",&fSelTrackRecoContained);
    fTree->Branch("SelTrackRecoMomMethod",&fSelTrackRecoMomMethod);
    fTree->Branch("SelTrackRecoCharge",&fSelTrackRecoCharge);
    fTree->Branch("SelTrackRecoMomMCS",&fSelTrackRecoMomMCS);
    fTree->Branch("SelTrackRecoMomContained",&fSelTrackRecoMomContained);
    //fTree->Branch("SelTrackRecoNMatchedVertices",&fSelTrackRecoNMatchedVertices);
    fTree->Branch("SelTrackRecoVertexX",&fSelTrackRecoVertexX);
    fTree->Branch("SelTrackRecoVertexY",&fSelTrackRecoVertexY);
    fTree->Branch("SelTrackRecoVertexZ",&fSelTrackRecoVertexZ);
    fTree->Branch("SelTrackRecoNChildPFP",&fSelTrackRecoNChildPFP);
    fTree->Branch("SelTrackRecoNChildTrackPFP",&fSelTrackRecoNChildTrackPFP);
    fTree->Branch("SelTrackRecoNChildShowerPFP",&fSelTrackRecoNChildShowerPFP);
    fTree->Branch("SelTrackMVAElectron",&fSelTrackMVAElectron);
    fTree->Branch("SelTrackMVAPion",&fSelTrackMVAPion);
    fTree->Branch("SelTrackMVAMuon",&fSelTrackMVAMuon);
    fTree->Branch("SelTrackMVAProton",&fSelTrackMVAProton);
    fTree->Branch("SelTrackMVAPhoton",&fSelTrackMVAPhoton);
    fTree->Branch("SelTrackDeepPanMuVar",&fSelTrackDeepPanMuVar);
    fTree->Branch("SelTrackDeepPanPiVar",&fSelTrackDeepPanPiVar);
    fTree->Branch("SelTrackDeepPanProtonVar",&fSelTrackDeepPanProtonVar);
    fTree->Branch("SelTrackMichelNHits", &fSelTrackMichelNHits);
    fTree->Branch("SelTrackMichelElectronMVA", &fSelTrackMichelElectronMVA);
    fTree->Branch("SelTrackMichelRecoEnergyPlane2", &fSelTrackMichelRecoEnergyPlane2);
    fTree->Branch("SelTrackDeflecAngleSD", &fSelTrackDeflecAngleSD);
    fTree->Branch("SelTrackLength", &fSelTrackLength);
    fTree->Branch("SelTrackEvalRatio", &fSelTrackEvalRatio);
    fTree->Branch("SelTrackConcentration", &fSelTrackConcentration);
    fTree->Branch("SelTrackCoreHaloRatio", &fSelTrackCoreHaloRatio);
    fTree->Branch("SelTrackConicalness", &fSelTrackConicalness);
    fTree->Branch("SelTrackdEdxStart", &fSelTrackdEdxStart);
    fTree->Branch("SelTrackdEdxEnd", &fSelTrackdEdxEnd);
    fTree->Branch("SelTrackdEdxEndRatio", &fSelTrackdEdxEndRatio);
    fTree->Branch("SelTrackPandizzleVar", &fSelTrackPandizzleVar);
    fTree->Branch("NRecoTracks",&fNRecoTracks);
    fTree->Branch("RecoTrackIsPrimary",fRecoTrackIsPrimary,"RecoTrackIsPrimary[NRecoTracks]/O");
    fTree->Branch("RecoTrackTruePDG",fRecoTrackTruePDG,"RecoTrackTruePDG[NRecoTracks]/I");
    fTree->Branch("RecoTrackTruePrimary",fRecoTrackTruePrimary,"RecoTrackTruePrimary[NRecoTracks]/O");
    fTree->Branch("RecoTrackTrueMomX",fRecoTrackTrueMomX,"RecoTrackTrueMomX[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueMomY",fRecoTrackTrueMomY,"RecoTrackTrueMomY[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueMomZ",fRecoTrackTrueMomZ,"RecoTrackTrueMomZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueMomT",fRecoTrackTrueMomT,"RecoTrackTrueMomT[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueStartX",fRecoTrackTrueStartX,"RecoTrackTrueStartX[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueStartY",fRecoTrackTrueStartY,"RecoTrackTrueStartY[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueStartZ",fRecoTrackTrueStartZ,"RecoTrackTrueStartZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueStartT",fRecoTrackTrueStartT,"RecoTrackTrueStartT[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueEndX",fRecoTrackTrueEndX,"RecoTrackTrueEndX[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueEndY",fRecoTrackTrueEndY,"RecoTrackTrueEndY[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueEndZ",fRecoTrackTrueEndZ,"RecoTrackTrueEndZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackTrueEndT",fRecoTrackTrueEndT,"RecoTrackTrueEndT[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoNHits",fRecoTrackRecoNHits,"RecoTrackRecoNHits[NRecoTracks]/I");
    fTree->Branch("RecoTrackRecoCompleteness",fRecoTrackRecoCompleteness,"RecoTrackRecoCompleteness[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoHitPurity",fRecoTrackRecoHitPurity,"RecoTrackRecoHitPurity[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomX",fRecoTrackRecoMomX,"RecoTrackRecoMomX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomY",fRecoTrackRecoMomY,"RecoTrackRecoMomY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomZ",fRecoTrackRecoMomZ,"RecoTrackRecoMomZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomT",fRecoTrackRecoMomT,"RecoTrackRecoMomT[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoStartX",fRecoTrackRecoStartX,"RecoTrackRecoStartX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoStartY",fRecoTrackRecoStartY,"RecoTrackRecoStartY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoStartZ",fRecoTrackRecoStartZ,"RecoTrackRecoStartZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoStartT",fRecoTrackRecoStartT,"RecoTrackRecoStartT[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndX",fRecoTrackRecoEndX,"RecoTrackRecoEndX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndY",fRecoTrackRecoEndY,"RecoTrackRecoEndY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndZ",fRecoTrackRecoEndZ,"RecoTrackRecoEndZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndT",fRecoTrackRecoEndT,"RecoTrackRecoEndT[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamX",fRecoTrackRecoUpstreamX,"RecoTrackRecoUpstreamX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamY",fRecoTrackRecoUpstreamY,"RecoTrackRecoUpstreamY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamZ",fRecoTrackRecoUpstreamZ,"RecoTrackRecoUpstreamZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoUpstreamT",fRecoTrackRecoUpstreamT,"RecoTrackRecoUpstreamT[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamX",fRecoTrackRecoDownstreamX,"RecoTrackRecoDownstreamX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamY",fRecoTrackRecoDownstreamY,"RecoTrackRecoDownstreamY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamZ",fRecoTrackRecoDownstreamZ,"RecoTrackRecoDownstreamZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoDownstreamT",fRecoTrackRecoDownstreamT,"RecoTrackRecoDownstreamT[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexX",fRecoTrackRecoEndClosestToVertexX,"RecoTrackRecoEndClosestToVertexX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexY",fRecoTrackRecoEndClosestToVertexY,"RecoTrackRecoEndClosestToVertexY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoEndClosestToVertexZ",fRecoTrackRecoEndClosestToVertexZ,"RecoTrackRecoEndClosestToVertexZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoLength",fRecoTrackRecoLength,"RecoTrackRecoLength[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoContained",fRecoTrackRecoContained,"RecoTrackRecoContained[NRecoTracks]/O");
    fTree->Branch("RecoTrackRecoCharge",fRecoTrackRecoCharge,"RecoTrackRecoCharge[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomMCS",fRecoTrackRecoMomMCS,"RecoTrackRecoMomMCS[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoMomContained",fRecoTrackRecoMomContained,"RecoTrackRecoMomContained[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoVertexX",fRecoTrackRecoVertexX,"RecoTrackRecoVertexX[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoVertexY",fRecoTrackRecoVertexY,"RecoTrackRecoVertexY[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoVertexZ",fRecoTrackRecoVertexZ,"RecoTrackRecoVertexZ[NRecoTracks]/D");
    fTree->Branch("RecoTrackRecoNChildPFP",fRecoTrackRecoNChildPFP,"RecoTrackRecoNChildPFP[NRecoTracks]/I");
    fTree->Branch("RecoTrackRecoNChildTrackPFP",fRecoTrackRecoNChildTrackPFP,"RecoTrackRecoNChildTrackPFP[NRecoTracks]/I");
    fTree->Branch("RecoTrackRecoNChildShowerPFP",fRecoTrackRecoNChildShowerPFP,"RecoTrackRecoNChildShowerPFP[NRecoTracks]/I");
    fTree->Branch("RecoTrackMVAElectron",fRecoTrackMVAElectron,"RecoTrackMVAElectron[NRecoTracks]/D");
    fTree->Branch("RecoTrackMVAPion",fRecoTrackMVAPion,"RecoTrackMVAPion[NRecoTracks]/D");
    fTree->Branch("RecoTrackMVAMuon",fRecoTrackMVAMuon,"RecoTrackMVAMuon[NRecoTracks]/D");
    fTree->Branch("RecoTrackMVAProton",fRecoTrackMVAProton,"RecoTrackMVAProton[NRecoTracks]/D");
    fTree->Branch("RecoTrackMVAPhoton",fRecoTrackMVAPhoton,"RecoTrackMVAPhoton[NRecoTracks]/D");
    fTree->Branch("RecoTrackDeepPanMuVar",fRecoTrackDeepPanMuVar,"RecoTrackDeepPanMuVar[NRecoTracks]/D");
    fTree->Branch("RecoTrackDeepPanPiVar",fRecoTrackDeepPanPiVar,"RecoTrackDeepPanPiVar[NRecoTracks]/D");
    fTree->Branch("RecoTrackDeepPanProtonVar",fRecoTrackDeepPanProtonVar,"RecoTrackDeepPanProtonVar[NRecoTracks]/D");
    fTree->Branch("RecoTrackMichelNHits", fRecoTrackMichelNHits, "RecoTrackMichelNHits/D");
    fTree->Branch("RecoTrackMichelElectronMVA", fRecoTrackMichelElectronMVA, "RecoTrackMichelElectronMVA/D");
    fTree->Branch("RecoTrackMichelRecoEnergyPlane2", fRecoTrackMichelRecoEnergyPlane2, "RecoTrackMichelRecoEnergyPlane2/D");
    fTree->Branch("RecoTrackDeflecAngleSD", fRecoTrackDeflecAngleSD, "RecoTrackDeflecAngleSD/D");
    fTree->Branch("RecoTrackLength", fRecoTrackLength, "RecoTrackLength/D");
    fTree->Branch("RecoTrackEvalRatio", fRecoTrackEvalRatio, "RecoTrackEvalRatio/D");
    fTree->Branch("RecoTrackConcentration", fRecoTrackConcentration, "RecoTrackConcentration/D");
    fTree->Branch("RecoTrackCoreHaloRatio", fRecoTrackCoreHaloRatio, "RecoTrackCoreHaloRatio/D");
    fTree->Branch("RecoTrackConicalness", fRecoTrackConicalness, "RecoTrackConicalness/D");
    fTree->Branch("RecoTrackdEdxStart", fRecoTrackdEdxStart, "RecoTrackdEdxStart/D");
    fTree->Branch("RecoTrackdEdxEnd", fRecoTrackdEdxEnd, "RecoTrackdEdxEnd/D");
    fTree->Branch("RecoTrackdEdxEndRatio", fRecoTrackdEdxEndRatio, "RecoTrackdEdxEndRatio/D");
    fTree->Branch("RecoTrackPandizzleVar",fRecoTrackPandizzleVar,"RecoTrackPandizzleVar[NRecoTracks]/D");

    fTree->Branch("SelShowerTruePDG",&fSelShowerTruePDG);
    fTree->Branch("SelShowerTruePrimary",&fSelShowerTruePrimary);
    fTree->Branch("SelShowerTrueMomX",&fSelShowerTrueMomX);
    fTree->Branch("SelShowerTrueMomY",&fSelShowerTrueMomY);
    fTree->Branch("SelShowerTrueMomZ",&fSelShowerTrueMomZ);
    fTree->Branch("SelShowerTrueMomT",&fSelShowerTrueMomT);
    fTree->Branch("SelShowerTrueStartX",&fSelShowerTrueStartX);
    fTree->Branch("SelShowerTrueStartY",&fSelShowerTrueStartY);
    fTree->Branch("SelShowerTrueStartZ",&fSelShowerTrueStartZ);
    fTree->Branch("SelShowerTrueStartT",&fSelShowerTrueStartT);
    fTree->Branch("SelShowerTrueEndX",&fSelShowerTrueEndX);
    fTree->Branch("SelShowerTrueEndY",&fSelShowerTrueEndY);
    fTree->Branch("SelShowerTrueEndZ",&fSelShowerTrueEndZ);
    fTree->Branch("SelShowerTrueEndT",&fSelShowerTrueEndT);
    fTree->Branch("SelShowerRecoNHits",&fSelShowerRecoNHits);
    fTree->Branch("SelShowerRecoCompleteness",&fSelShowerRecoCompleteness);
    fTree->Branch("SelShowerRecoHitPurity",&fSelShowerRecoHitPurity);
    fTree->Branch("SelShowerRecoMomX",&fSelShowerRecoMomX);
    fTree->Branch("SelShowerRecoMomY",&fSelShowerRecoMomY);
    fTree->Branch("SelShowerRecoMomZ",&fSelShowerRecoMomZ);
    fTree->Branch("SelShowerRecoMomT",&fSelShowerRecoMomT);
    fTree->Branch("SelShowerRecoStartX",&fSelShowerRecoStartX);
    fTree->Branch("SelShowerRecoStartY",&fSelShowerRecoStartY);
    fTree->Branch("SelShowerRecoStartZ",&fSelShowerRecoStartZ);
    fTree->Branch("SelShowerRecoStartT",&fSelShowerRecoStartT);
    fTree->Branch("SelShowerRecoEndX",&fSelShowerRecoEndX);
    fTree->Branch("SelShowerRecoEndY",&fSelShowerRecoEndY);
    fTree->Branch("SelShowerRecoEndZ",&fSelShowerRecoEndZ);
    fTree->Branch("SelShowerRecoEndT",&fSelShowerRecoEndT);
    fTree->Branch("SelShowerRecoUpstreamX",&fSelShowerRecoUpstreamX);
    fTree->Branch("SelShowerRecoUpstreamY",&fSelShowerRecoUpstreamY);
    fTree->Branch("SelShowerRecoUpstreamZ",&fSelShowerRecoUpstreamZ);
    fTree->Branch("SelShowerRecoUpstreamT",&fSelShowerRecoUpstreamT);
    fTree->Branch("SelShowerRecoDownstreamX",&fSelShowerRecoDownstreamX);
    fTree->Branch("SelShowerRecoDownstreamY",&fSelShowerRecoDownstreamY);
    fTree->Branch("SelShowerRecoDownstreamZ",&fSelShowerRecoDownstreamZ);
    fTree->Branch("SelShowerRecoDownstreamT",&fSelShowerRecoDownstreamT);
    fTree->Branch("SelShowerRecoEndClosestToVertexX",&fSelShowerRecoEndClosestToVertexX);
    fTree->Branch("SelShowerRecoEndClosestToVertexY",&fSelShowerRecoEndClosestToVertexY);
    fTree->Branch("SelShowerRecoEndClosestToVertexZ",&fSelShowerRecoEndClosestToVertexZ);
    fTree->Branch("SelShowerRecoContained",&fSelShowerRecoContained);
    fTree->Branch("SelShowerRecoCharge",&fSelShowerRecoCharge);
    fTree->Branch("SelShowerRecoMomMCS",&fSelShowerRecoMomMCS);
    fTree->Branch("SelShowerRecoMomContained",&fSelShowerRecoMomContained);
    fTree->Branch("SelShowerRecoNMatchedVertices",&fSelShowerRecoNMatchedVertices);
    fTree->Branch("SelShowerRecoVertexX",&fSelShowerRecoVertexX);
    fTree->Branch("SelShowerRecoVertexY",&fSelShowerRecoVertexY);
    fTree->Branch("SelShowerRecoVertexZ",&fSelShowerRecoVertexZ);
    fTree->Branch("SelShowerRecoNChildPFP",&fSelShowerRecoNChildPFP);
    fTree->Branch("SelShowerRecoNChildTrackPFP",&fSelShowerRecoNChildTrackPFP);
    fTree->Branch("SelShowerRecoNChildShowerPFP",&fSelShowerRecoNChildShowerPFP);
    fTree->Branch("SelShowerRecoDirX",&fSelShowerRecoDirX);
    fTree->Branch("SelShowerRecoDirY",&fSelShowerRecoDirY);
    fTree->Branch("SelShowerRecoDirZ",&fSelShowerRecoDirZ);
    fTree->Branch("SelShowerRecodEdx",&fSelShowerRecodEdx,"SelShowerRecodEdx[3]/D");
    fTree->Branch("SelShowerRecoBestPlane",&fSelShowerRecoBestPlane);
    fTree->Branch("SelShowerRecoLength",&fSelShowerRecoLength);
    fTree->Branch("SelShowerRecoOpeningAngle",&fSelShowerRecoOpeningAngle);
    fTree->Branch("SelShowerRecoIsPrimaryPFPDaughter",&fSelShowerRecoIsPrimaryPFPDaughter);
    fTree->Branch("SelShowerMVAElectron",&fSelShowerMVAElectron);
    fTree->Branch("SelShowerMVAPion",&fSelShowerMVAPion);
    fTree->Branch("SelShowerMVAMuon",&fSelShowerMVAMuon);
    fTree->Branch("SelShowerMVAProton",&fSelShowerMVAProton);
    fTree->Branch("SelShowerMVAPhoton",&fSelShowerMVAPhoton);
    fTree->Branch("SelShowerRecoEnergy",&fSelShowerRecoEnergy,"SelShowerRecoEnergy[3]/D");
    fTree->Branch("SelShowerPandrizzleConnectionBDTScore", &fSelShowerPandrizzleConnectionBDTScore);
    fTree->Branch("SelShowerPandrizzlePathwayLengthMin", &fSelShowerPandrizzlePathwayLengthMin);
    fTree->Branch("SelShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D", &fSelShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D);
    fTree->Branch("SelShowerPandrizzleMaxNPostShowerStartHits", &fSelShowerPandrizzleMaxNPostShowerStartHits);
    fTree->Branch("SelShowerPandrizzleMaxPostShowerStartScatterAngle", &fSelShowerPandrizzleMaxPostShowerStartScatterAngle);
    fTree->Branch("SelShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry", &fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry);
    fTree->Branch("SelShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry", &fSelShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry);
    fTree->Branch("SelShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", &fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
    fTree->Branch("SelShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius", &fSelShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius);
    fTree->Branch("SelShowerPandrizzleMaxPostShowerStartOpeningAngle", &fSelShowerPandrizzleMaxPostShowerStartOpeningAngle);
    fTree->Branch("SelShowerPandrizzleMaxFoundHitRatio", &fSelShowerPandrizzleMaxFoundHitRatio);
    fTree->Branch("SelShowerPandrizzleMaxInitialGapSize", &fSelShowerPandrizzleMaxInitialGapSize);
    fTree->Branch("SelShowerPandrizzleMinLargestProjectedGapSize", &fSelShowerPandrizzleMinLargestProjectedGapSize);
    fTree->Branch("SelShowerPandrizzleNViewsWithAmbiguousHits", &fSelShowerPandrizzleNViewsWithAmbiguousHits);
    fTree->Branch("SelShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy", &fSelShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy);
    fTree->Branch("SelShowerPandrizzleEvalRatio",&fSelShowerPandrizzleEvalRatio);
    fTree->Branch("SelShowerPandrizzleConcentration",&fSelShowerPandrizzleConcentration);
    fTree->Branch("SelShowerPandrizzleCoreHaloRatio",&fSelShowerPandrizzleCoreHaloRatio);
    fTree->Branch("SelShowerPandrizzleConicalness",&fSelShowerPandrizzleConicalness);
    fTree->Branch("SelShowerPandrizzledEdxBestPlane",&fSelShowerPandrizzledEdxBestPlane);
    fTree->Branch("SelShowerPandrizzleDisplacement",&fSelShowerPandrizzleDisplacement);
    fTree->Branch("SelShowerPandrizzleDCA",&fSelShowerPandrizzleDCA);
    fTree->Branch("SelShowerPandrizzleWideness",&fSelShowerPandrizzleWideness);
    fTree->Branch("SelShowerPandrizzleEnergyDensity",&fSelShowerPandrizzleEnergyDensity);
    fTree->Branch("SelShowerPandrizzleBDTMethod", &fSelShowerPandrizzleBDTMethod);
    fTree->Branch("SelShowerJamPandrizzleMVAScore",&fSelShowerJamPandrizzleMVAScore);
    fTree->Branch("SelShowerPandrizzleMVAScore",&fSelShowerPandrizzleMVAScore);
    fTree->Branch("SelShowerPandrizzleIsFilled",&fSelShowerPandrizzleIsFilled);
    fTree->Branch("SelShowerDistanceToNuVertexU", &fSelShowerDistanceToNuVertexU);
    fTree->Branch("SelShowerDistanceToNuVertexV", &fSelShowerDistanceToNuVertexV);
    fTree->Branch("SelShowerDistanceToNuVertexW", &fSelShowerDistanceToNuVertexW);

    fTree->Branch("NRecoShowers",&fNRecoShowers);
    fTree->Branch("RecoShowerTruePDG",fRecoShowerTruePDG,"RecoShowerTruePDG[NRecoShowers]/I");
    fTree->Branch("RecoShowerTruePrimary",fRecoShowerTruePrimary,"RecoShowerTruePrimary[NRecoShowers]/O");
    fTree->Branch("RecoShowerTrueMomX",fRecoShowerTrueMomX,"RecoShowerTrueMomX[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueMomY",fRecoShowerTrueMomY,"RecoShowerTrueMomY[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueMomZ",fRecoShowerTrueMomZ,"RecoShowerTrueMomZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueMomT",fRecoShowerTrueMomT,"RecoShowerTrueMomT[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueStartX",fRecoShowerTrueStartX,"RecoShowerTrueStartX[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueStartY",fRecoShowerTrueStartY,"RecoShowerTrueStartY[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueStartZ",fRecoShowerTrueStartZ,"RecoShowerTrueStartZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueStartT",fRecoShowerTrueStartT,"RecoShowerTrueStartT[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueEndX",fRecoShowerTrueEndX,"RecoShowerTrueEndX[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueEndY",fRecoShowerTrueEndY,"RecoShowerTrueEndY[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueEndZ",fRecoShowerTrueEndZ,"RecoShowerTrueEndZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerTrueEndT",fRecoShowerTrueEndT,"RecoShowerTrueEndT[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoNHits",fRecoShowerRecoNHits,"RecoShowerRecoNHits[NRecoShowers]/I");
    fTree->Branch("RecoShowerActiveCheatElectronAlg", fRecoShowerActiveCheatElectronAlg, "RecoShowerActiveCheatElectronAlg[NRecoShowers]/I");
    fTree->Branch("RecoShowerActiveHybridElectronAlg", fRecoShowerActiveHybridElectronAlg, "RecoShowerActiveHybridElectronAlg[NRecoShowers]/I");
    fTree->Branch("RecoShowerActiveCheatGammaAlg", fRecoShowerActiveCheatGammaAlg, "RecoShowerActiveCheatGammaAlg[NRecoShowers]/I");
    fTree->Branch("RecoShowerActiveHybridGammaAlg", fRecoShowerActiveHybridGammaAlg, "RecoShowerActiveHybridGammaAlg[NRecoShowers]/I");
    fTree->Branch("RecoShowerRecoCompleteness",fRecoShowerRecoCompleteness,"RecoShowerRecoCompleteness[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoHitPurity",fRecoShowerRecoHitPurity,"RecoShowerRecoHitPurity[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoMomX",fRecoShowerRecoMomX,"RecoShowerRecoMomX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoMomY",fRecoShowerRecoMomY,"RecoShowerRecoMomY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoMomZ",fRecoShowerRecoMomZ,"RecoShowerRecoMomZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoMomT",fRecoShowerRecoMomT,"RecoShowerRecoMomT[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoStartX",fRecoShowerRecoStartX,"RecoShowerRecoStartX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoStartY",fRecoShowerRecoStartY,"RecoShowerRecoStartY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoStartZ",fRecoShowerRecoStartZ,"RecoShowerRecoStartZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoStartT",fRecoShowerRecoStartT,"RecoShowerRecoStartT[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndX",fRecoShowerRecoEndX,"RecoShowerRecoEndX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndY",fRecoShowerRecoEndY,"RecoShowerRecoEndY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndZ",fRecoShowerRecoEndZ,"RecoShowerRecoEndZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndT",fRecoShowerRecoEndT,"RecoShowerRecoEndT[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoUpstreamX",fRecoShowerRecoUpstreamX,"RecoShowerRecoUpstreamX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoUpstreamY",fRecoShowerRecoUpstreamY,"RecoShowerRecoUpstreamY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoUpstreamZ",fRecoShowerRecoUpstreamZ,"RecoShowerRecoUpstreamZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoUpstreamT",fRecoShowerRecoUpstreamT,"RecoShowerRecoUpstreamT[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDownstreamX",fRecoShowerRecoDownstreamX,"RecoShowerRecoDownstreamX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDownstreamY",fRecoShowerRecoDownstreamY,"RecoShowerRecoDownstreamY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDownstreamZ",fRecoShowerRecoDownstreamZ,"RecoShowerRecoDownstreamZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDownstreamT",fRecoShowerRecoDownstreamT,"RecoShowerRecoDownstreamT[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndClosestToVertexX",fRecoShowerRecoEndClosestToVertexX,"RecoShowerRecoEndClosestToVertexX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndClosestToVertexY",fRecoShowerRecoEndClosestToVertexY,"RecoShowerRecoEndClosestToVertexY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEndClosestToVertexZ",fRecoShowerRecoEndClosestToVertexZ,"RecoShowerRecoEndClosestToVertexZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoContained",fRecoShowerRecoContained,"RecoShowerRecoContained[NRecoShowers]/O");
    fTree->Branch("RecoShowerRecoCharge",fRecoShowerRecoCharge,"RecoShowerRecoCharge[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoMomMCS",fRecoShowerRecoMomMCS,"RecoShowerRecoMomMCS[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoMomContained",fRecoShowerRecoMomContained,"RecoShowerRecoMomContained[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoNMatchedVertices",fRecoShowerRecoNMatchedVertices,"RecoShowerRecoNMatchedVertices[NRecoShowers]/I");
    fTree->Branch("RecoShowerRecoVertexX",fRecoShowerRecoVertexX,"RecoShowerRecoVertexX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoVertexY",fRecoShowerRecoVertexY,"RecoShowerRecoVertexY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoVertexZ",fRecoShowerRecoVertexZ,"RecoShowerRecoVertexZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoNChildPFP",fRecoShowerRecoNChildPFP,"RecoShowerRecoNChildPFP[NRecoShowers]/I");
    fTree->Branch("RecoShowerRecoNChildTrackPFP",fRecoShowerRecoNChildTrackPFP,"RecoShowerRecoNChildTrackPFP[NRecoShowers]/I");
    fTree->Branch("RecoShowerRecoNChildShowerPFP",fRecoShowerRecoNChildShowerPFP,"RecoShowerRecoNChildShowerPFP[NRecoShowers]/I");
    fTree->Branch("RecoShowerRecoDirX",fRecoShowerRecoDirX,"RecoShowerRecoDirX[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDirY",fRecoShowerRecoDirY,"RecoShowerRecoDirY[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoDirZ",fRecoShowerRecoDirZ,"RecoShowerRecoDirZ[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecodEdx",fRecoShowerRecodEdx,"RecoShowerRecodEdx[NRecoShowers][3]/D");
    fTree->Branch("RecoShowerRecoBestPlane",&fRecoShowerRecoBestPlane,"RecoShowerRecoBestPlane[3]/I");
    fTree->Branch("RecoShowerRecoLength",&fRecoShowerRecoLength,"RecoShowerRecoLength[3]/D");
    fTree->Branch("RecoShowerRecoOpeningAngle",&fRecoShowerRecoOpeningAngle,"RecoShowerRecoOpeningAngle[3]/D");
    fTree->Branch("RecoShowerRecoIsPrimaryPFPDaughter",fRecoShowerRecoIsPrimaryPFPDaughter,"RecoShowerRecoIsPrimaryPFPDaughter[NRecoShowers]/O");
    fTree->Branch("RecoShowerMVAElectron",fRecoShowerMVAElectron,"RecoShowerMVAElectron[NRecoShowers]/D");
    fTree->Branch("RecoShowerMVAPion",fRecoShowerMVAPion,"RecoShowerMVAPion[NRecoShowers]/D");
    fTree->Branch("RecoShowerMVAMuon",fRecoShowerMVAMuon,"RecoShowerMVAMuon[NRecoShowers]/D");
    fTree->Branch("RecoShowerMVAProton",fRecoShowerMVAProton,"RecoShowerMVAProton[NRecoShowers]/D");
    fTree->Branch("RecoShowerMVAPhoton",fRecoShowerMVAPhoton,"RecoShowerMVAPhoton[NRecoShowers]/D");
    fTree->Branch("RecoShowerRecoEnergy",fRecoShowerRecoEnergy,"RecoShowerRecoEnergy[NRecoShowers][3]/D");

    fTree->Branch("RecoShowerPandrizzleConnectionBDTScore", &fRecoShowerPandrizzleConnectionBDTScore, "RecoShowerPandrizzleConnectionBDTScore[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzlePathwayLengthMin", &fRecoShowerPandrizzlePathwayLengthMin, "RecoShowerPandrizzlePathwayLengthMin[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D", &fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D, 
      "RecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxNPostShowerStartHits", &fRecoShowerPandrizzleMaxNPostShowerStartHits, "RecoShowerPandrizzleMaxNPostShowerStartHits[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartScatterAngle", &fRecoShowerPandrizzleMaxPostShowerStartScatterAngle, "RecoShowerPandrizzleMaxPostShowerStartScatterAngle[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry", &fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry, 
      "RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry", &fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry, 
      "RecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", &fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance, 
      "RecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius", &fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius, 
      "RecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxPostShowerStartOpeningAngle", &fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle, "RecoShowerPandrizzleMaxPostShowerStartOpeningAngle[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxFoundHitRatio", &fRecoShowerPandrizzleMaxFoundHitRatio, "RecoShowerPandrizzleMaxFoundHitRatio[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMaxInitialGapSize", &fRecoShowerPandrizzleMaxInitialGapSize, "RecoShowerPandrizzleMaxInitialGapSize[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleMinLargestProjectedGapSize", &fRecoShowerPandrizzleMinLargestProjectedGapSize, "RecoShowerPandrizzleMinLargestProjectedGapSize[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleNViewsWithAmbiguousHits", &fRecoShowerPandrizzleNViewsWithAmbiguousHits, "RecoShowerPandrizzleNViewsWithAmbiguousHits[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy", &fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy, "RecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleEvalRatio",&fRecoShowerPandrizzleEvalRatio,"RecoShowerPandrizzleEvalRatio[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleConcentration",&fRecoShowerPandrizzleConcentration,"RecoShowerPandrizzleConcentration[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleCoreHaloRatio",&fRecoShowerPandrizzleCoreHaloRatio, "RecoShowerPandrizzleCoreHaloRatio[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleConicalness",&fRecoShowerPandrizzleConicalness, "RecoShowerPandrizzleConicalness[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzledEdxBestPlane",&fRecoShowerPandrizzledEdxBestPlane, "RecoShowerPandrizzledEdxBestPlane[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleDisplacement",&fRecoShowerPandrizzleDisplacement, "RecoShowerPandrizzleDisplacement[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleDCA",&fRecoShowerPandrizzleDCA, "RecoShowerPandrizzleDCA[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleWideness",&fRecoShowerPandrizzleWideness, "RecoShowerPandrizzleWideness[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleEnergyDensity",&fRecoShowerPandrizzleEnergyDensity, "RecoShowerPandrizzleEnergyDensity[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleBDTMethod", &fRecoShowerPandrizzleBDTMethod, "RecoShowerPandrizzleBDTMethod[NRecoShowers]/D");
    fTree->Branch("RecoShowerJamPandrizzleMVAScore",&fRecoShowerJamPandrizzleMVAScore, "RecoShowerJamPandrizzleMVAScore[NRecoShowers]/D"); 
    fTree->Branch("RecoShowerPandrizzleMVAScore",&fRecoShowerPandrizzleMVAScore, "RecoShowerPandrizzleMVAScore[NRecoShowers]/D");
    fTree->Branch("RecoShowerPandrizzleIsFilled",&fRecoShowerPandrizzleIsFilled, "RecoShowerPandrizzleIsFilled[NRecoShowers]/O");
    fTree->Branch("RecoShowerDistanceToNuVertexU", &fRecoShowerDistanceToNuVertexU, "RecoShowerDistanceToNuVertexU[NRecoShowers]/D");
    fTree->Branch("RecoShowerDistanceToNuVertexV", &fRecoShowerDistanceToNuVertexV, "RecoShowerDistanceToNuVertexV[NRecoShowers]/D");
    fTree->Branch("RecoShowerDistanceToNuVertexW", &fRecoShowerDistanceToNuVertexW, "RecoShowerDistanceToNuVertexW[NRecoShowers]/D");

    fTree->Branch("CVNResultNue", &fCVNResultNue);
    fTree->Branch("CVNResultNumu", &fCVNResultNumu);
    fTree->Branch("CVNResultNutau", &fCVNResultNutau);
    fTree->Branch("CVNResultNC", &fCVNResultNC);

    fTree->Branch("RecoNuVtxX",&fRecoNuVtxX);
    fTree->Branch("RecoNuVtxY",&fRecoNuVtxY);
    fTree->Branch("RecoNuVtxZ",&fRecoNuVtxZ);
    fTree->Branch("RecoNuVtxNShowers",&fRecoNuVtxNShowers);
    fTree->Branch("RecoNuVtxNTracks",&fRecoNuVtxNTracks);
    fTree->Branch("RecoNuVtxNChildren",&fRecoNuVtxNChildren);

    fTree->Branch("RecoEventCharge",&fRecoEventCharge);
    fTree->Branch("NumuRecoMomLep",&fNumuRecoMomLep);
    fTree->Branch("NumuRecoEHad",&fNumuRecoEHad);
    fTree->Branch("NumuRecoENu",&fNumuRecoENu);
    fTree->Branch("NueRecoMomLep",&fNueRecoMomLep);
    fTree->Branch("NueRecoEHad",&fNueRecoEHad);
    fTree->Branch("NueRecoENu",&fNueRecoENu);



    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("POT",&fPOT);
    fPOTTree->Branch("Run",&fRun);
    fPOTTree->Branch("SubRun",&fSubRun);

    Reset();  //Default value all variables now
}

void FDSelection::CCNuSelection::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}
void FDSelection::CCNuSelection::endSubRun(const art::SubRun& sr){
  //Need the run and subrun
  fRun = sr.run();
  fSubRun = sr.subRun();
  //Need the POT (obvs)
  art::Handle< sumdata::POTSummary > potListHandle;

  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    fPOT = potListHandle->totpot;
  else
    fPOT = 0.;
  if (fPOTTree) fPOTTree->Fill();
}

void FDSelection::CCNuSelection::endJob()
{
  // Implementation of optional member function here.
}

void FDSelection::CCNuSelection::Reset()
{
  //Generic stuff
  fRun = kDefInt;
  fSubRun = kDefInt;
  fEvent = kDefInt;
  fIsMC = kDefInt;
  //Detector stuff
  fT0 = kDefDoub;
  //Neutrino stuff
  fNuPdg = kDefInt; 
  fBeamPdg = kDefInt; 
  fNC = kDefInt;    
  fMode = kDefInt; 
  fTargetZ = kDefInt;
  fQ2 = kDefDoub; 
  fENu = kDefDoub; 
  fW = kDefDoub; 
  fX = kDefDoub;
  fY = kDefDoub;
  fNuMomX = kDefDoub; 
  fNuMomY = kDefDoub;
  fNuMomZ = kDefDoub;
  fNuMomT = kDefDoub;
  fNuX = kDefDoub; 
  fNuY = kDefDoub;
  fNuZ = kDefDoub;
  fNuT = kDefDoub;
  fNPiP = 0;
  fNPim = 0;
  fNPi0 = 0;
  fNPhotons = 0;
  fNProtons = 0;
  fNNeutrons = 0;
  fNOther = 0;
  for (int i_vertexparticle = 0; i_vertexparticle < (fNVertexParticles == 0 ? kDefMaxNTrueVertexParticles : fNVertexParticles); i_vertexparticle++){
    fVertexParticleIsGHEP[i_vertexparticle] = 0;
    fVertexParticlePDG[i_vertexparticle] = kDefInt;
    fVertexParticleStatus[i_vertexparticle] = kDefInt;
    fVertexParticleNChildren[i_vertexparticle] = kDefInt;
    fVertexParticleMomX[i_vertexparticle] = kDefDoub;
    fVertexParticleMomY[i_vertexparticle] = kDefDoub;
    fVertexParticleMomZ[i_vertexparticle] = kDefDoub;
    fVertexParticleMomT[i_vertexparticle] = kDefDoub;
    fVertexParticleEndX[i_vertexparticle] = kDefDoub;
    fVertexParticleEndY[i_vertexparticle] = kDefDoub;
    fVertexParticleEndZ[i_vertexparticle] = kDefDoub;
    fVertexParticleEndT[i_vertexparticle] = kDefDoub;
  }
  fNVertexParticles = 0;
  //Outgoing lepton stuff
  fLepPDG = kDefInt;
  fMomLepX = kDefDoub;
  fMomLepY = kDefDoub;
  fMomLepZ = kDefDoub;
  fMomLepT = kDefDoub;
  fLepEndX = kDefDoub;
  fLepEndY = kDefDoub;
  fLepEndZ = kDefDoub;
  fLepEndT = kDefDoub;
  fLepNuAngle = kDefDoub;
  fNuMomTranMag = kDefDoub;
  fTargNuclMomTranMag = kDefDoub;
  fInitalMomTranMag = kDefDoub;
  fLepMomTranMag = kDefDoub;
  fNuclRemMomTranMag = kDefDoub;
  fFinalMomTranMagNoLepNoRem = kDefDoub;
  fFinalMomTranMagNoLepWithRem = kDefDoub;
  fFinalMomTranMagWithLepNoRem = kDefDoub;
  fFinalMomTranMagWithLepWithRem = kDefDoub;
  //Selection stuff
  //Track stuff
  //true bits
  fSelTrackTruePDG = kDefInt;
  fSelTrackTruePrimary = kDefInt;
  fSelTrackTrueMomX = kDefDoub;
  fSelTrackTrueMomY = kDefDoub;
  fSelTrackTrueMomZ = kDefDoub;
  fSelTrackTrueMomT = kDefDoub;
  fSelTrackTrueStartX = kDefDoub;
  fSelTrackTrueStartY = kDefDoub;
  fSelTrackTrueStartZ = kDefDoub;
  fSelTrackTrueStartT = kDefDoub;
  fSelTrackTrueEndX = kDefDoub;
  fSelTrackTrueEndY = kDefDoub;
  fSelTrackTrueEndZ = kDefDoub;
  fSelTrackTrueEndT = kDefDoub;
  //reco bits
  fSelTrackRecoNHits = kDefInt;
  fSelTrackRecoCompleteness = kDefDoub;
  fSelTrackRecoHitPurity = kDefDoub;
  fSelTrackRecoMomX = kDefDoub;
  fSelTrackRecoMomY = kDefDoub;
  fSelTrackRecoMomZ = kDefDoub;
  fSelTrackRecoMomT = kDefDoub;
  fSelTrackRecoStartX = kDefDoub;
  fSelTrackRecoStartY = kDefDoub;
  fSelTrackRecoStartZ = kDefDoub;
  fSelTrackRecoStartT = kDefDoub;
  fSelTrackRecoEndX = kDefDoub;
  fSelTrackRecoEndY = kDefDoub;
  fSelTrackRecoEndZ = kDefDoub;
  fSelTrackRecoEndT = kDefDoub;
  fSelTrackRecoUpstreamX = kDefDoub;
  fSelTrackRecoUpstreamY = kDefDoub;
  fSelTrackRecoUpstreamZ = kDefDoub;
  fSelTrackRecoUpstreamT = kDefDoub;
  fSelTrackRecoDownstreamX = kDefDoub;
  fSelTrackRecoDownstreamY = kDefDoub;
  fSelTrackRecoDownstreamZ = kDefDoub;
  fSelTrackRecoDownstreamT = kDefDoub;
  fSelTrackRecoEndClosestToVertexX = kDefDoub;
  fSelTrackRecoEndClosestToVertexY = kDefDoub;
  fSelTrackRecoEndClosestToVertexZ = kDefDoub;
  fSelTrackRecoLength = kDefDoub;
  fSelTrackRecoContained = kDefInt;
  fSelTrackRecoMomMethod = kDefInt;
  fSelTrackRecoCharge = kDefDoub;
  fSelTrackRecoMomMCS = kDefDoub;
  fSelTrackRecoMomContained = kDefDoub;
  //fSelTrackRecoNMatchedVertices = kDefInt;
  fSelTrackRecoVertexX = kDefDoub;
  fSelTrackRecoVertexY = kDefDoub;
  fSelTrackRecoVertexZ = kDefDoub;
  fSelTrackRecoNChildPFP = kDefInt;
  fSelTrackRecoNChildTrackPFP = kDefInt;
  fSelTrackRecoNChildShowerPFP = kDefInt;
  //MVA bits
  fSelTrackMVAElectron = kDefDoub;
  fSelTrackMVAPion = kDefDoub;
  fSelTrackMVAMuon = kDefDoub;
  fSelTrackMVAProton = kDefDoub;
  fSelTrackMVAPhoton = kDefDoub;
  //Pandizzle
  fSelTrackMichelNHits = kDefDoub;
  fSelTrackMichelElectronMVA = kDefDoub;
  fSelTrackMichelRecoEnergyPlane2 = kDefDoub;
  fSelTrackDeflecAngleSD = kDefDoub;
  fSelTrackLength = kDefDoub;
  fSelTrackEvalRatio = kDefDoub;
  fSelTrackConcentration = kDefDoub;
  fSelTrackCoreHaloRatio = kDefDoub;
  fSelTrackConicalness = kDefDoub;
  fSelTrackdEdxStart = kDefDoub;
  fSelTrackdEdxEnd = kDefDoub;
  fSelTrackdEdxEndRatio = kDefDoub;
  fSelTrackPandizzleVar = kDefDoub;
  //DeepPan
  fSelTrackDeepPanMuVar = kDefDoub;
  fSelTrackDeepPanPiVar = kDefDoub;
  fSelTrackDeepPanProtonVar = kDefDoub;
  //CVN
  fCVNResultNue = kDefDoub;
  fCVNResultNumu = kDefDoub;
  fCVNResultNutau = kDefDoub;
  fCVNResultNC = kDefDoub;

  for (int i_recotrack = 0; i_recotrack < (fNRecoTracks == 0 ? kDefMaxNRecoTracks : fNRecoTracks); i_recotrack++){
    fRecoTrackIsPrimary[i_recotrack] = false;
    fRecoTrackTruePDG[i_recotrack] = kDefInt;
    fRecoTrackTruePrimary[i_recotrack] = false;
    fRecoTrackTrueMomX[i_recotrack] = kDefDoub;
    fRecoTrackTrueMomY[i_recotrack] = kDefDoub;
    fRecoTrackTrueMomZ[i_recotrack] = kDefDoub;
    fRecoTrackTrueMomT[i_recotrack] = kDefDoub;
    fRecoTrackTrueStartX[i_recotrack] = kDefDoub;
    fRecoTrackTrueStartY[i_recotrack] = kDefDoub;
    fRecoTrackTrueStartZ[i_recotrack] = kDefDoub;
    fRecoTrackTrueStartT[i_recotrack] = kDefDoub;
    fRecoTrackTrueEndX[i_recotrack] = kDefDoub;
    fRecoTrackTrueEndY[i_recotrack] = kDefDoub;
    fRecoTrackTrueEndZ[i_recotrack] = kDefDoub;
    fRecoTrackTrueEndT[i_recotrack] = kDefDoub;
    fRecoTrackRecoNHits[i_recotrack] = kDefInt;
    fRecoTrackRecoCompleteness[i_recotrack] = kDefDoub;
    fRecoTrackRecoHitPurity[i_recotrack] = kDefDoub;
    fRecoTrackRecoMomX[i_recotrack] = kDefDoub;
    fRecoTrackRecoMomY[i_recotrack] = kDefDoub;
    fRecoTrackRecoMomZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoMomT[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartX[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartY[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartT[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndX[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndY[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndT[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamX[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamY[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamT[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamX[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamY[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamT[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndClosestToVertexX[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndClosestToVertexY[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndClosestToVertexZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoLength[i_recotrack] = kDefDoub;
    fRecoTrackRecoContained[i_recotrack] = 0;
    fRecoTrackRecoCharge[i_recotrack] = kDefDoub;
    fRecoTrackRecoMomMCS[i_recotrack] = kDefDoub;
    fRecoTrackRecoMomContained[i_recotrack] = kDefDoub;
    fRecoTrackRecoVertexX[i_recotrack] = kDefDoub;
    fRecoTrackRecoVertexY[i_recotrack] = kDefDoub;
    fRecoTrackRecoVertexZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoNChildPFP[i_recotrack] = kDefInt;
    fRecoTrackRecoNChildTrackPFP[i_recotrack] = kDefInt;
    fRecoTrackRecoNChildShowerPFP[i_recotrack] = kDefInt;

    fRecoTrackMVAElectron[i_recotrack] = kDefDoub;
    fRecoTrackMVAPion[i_recotrack] = kDefDoub;
    fRecoTrackMVAMuon[i_recotrack] = kDefDoub;
    fRecoTrackMVAProton[i_recotrack] = kDefDoub;
    fRecoTrackMVAPhoton[i_recotrack] = kDefDoub;
    fRecoTrackPandizzleVar[i_recotrack] = kDefDoub;
    fRecoTrackDeepPanMuVar[i_recotrack] = kDefDoub;
    fRecoTrackDeepPanPiVar[i_recotrack] = kDefDoub;
    fRecoTrackDeepPanProtonVar[i_recotrack] = kDefDoub;

    fRecoTrackMichelNHits[i_recotrack] = kDefDoub;
    fRecoTrackMichelElectronMVA[i_recotrack] = kDefDoub;
    fRecoTrackMichelRecoEnergyPlane2[i_recotrack] = kDefDoub;
    fRecoTrackDeflecAngleSD[i_recotrack] = kDefDoub;
    fRecoTrackLength[i_recotrack] = kDefDoub;
    fRecoTrackEvalRatio[i_recotrack] = kDefDoub;
    fRecoTrackConcentration[i_recotrack] = kDefDoub;
    fRecoTrackCoreHaloRatio[i_recotrack] = kDefDoub;
    fRecoTrackConicalness[i_recotrack] = kDefDoub;
    fRecoTrackdEdxStart[i_recotrack] = kDefDoub;
    fRecoTrackdEdxEnd[i_recotrack] = kDefDoub;
    fRecoTrackdEdxEndRatio[i_recotrack] = kDefDoub;
  }
  fNRecoTracks=0;

  //Reco energy bits
  fNumuRecoEHad = kDefDoub;
  fNumuRecoMomLep = kDefDoub;
  fNumuRecoENu = kDefDoub; //Neutrino reco energy

  //Shower bits
  //true bits
  fSelShowerTruePDG = kDefInt;
  fSelShowerTruePrimary = false;
  fSelShowerTrueMomX = kDefDoub;
  fSelShowerTrueMomY = kDefDoub;
  fSelShowerTrueMomZ = kDefDoub;
  fSelShowerTrueMomT = kDefDoub;
  fSelShowerTrueStartX = kDefDoub;
  fSelShowerTrueStartY = kDefDoub;
  fSelShowerTrueStartZ = kDefDoub;
  fSelShowerTrueStartT = kDefDoub;
  fSelShowerTrueEndX = kDefDoub;
  fSelShowerTrueEndY = kDefDoub;
  fSelShowerTrueEndZ = kDefDoub;
  fSelShowerTrueEndT = kDefDoub;
  //reco bits
  fSelShowerRecoNHits = kDefInt;
  fSelShowerRecoCompleteness = kDefDoub;
  fSelShowerRecoHitPurity = kDefDoub;
  fSelShowerRecoMomX = kDefDoub;
  fSelShowerRecoMomY = kDefDoub;
  fSelShowerRecoMomZ = kDefDoub;
  fSelShowerRecoMomT = kDefDoub;
  fSelShowerRecoStartX = kDefDoub;
  fSelShowerRecoStartY = kDefDoub;
  fSelShowerRecoStartZ = kDefDoub;
  fSelShowerRecoStartT = kDefDoub;
  fSelShowerRecoEndX = kDefDoub;
  fSelShowerRecoEndY = kDefDoub;
  fSelShowerRecoEndZ = kDefDoub;
  fSelShowerRecoEndT = kDefDoub;
  fSelShowerRecoUpstreamX = kDefDoub;
  fSelShowerRecoUpstreamY = kDefDoub;
  fSelShowerRecoUpstreamZ = kDefDoub;
  fSelShowerRecoUpstreamT = kDefDoub;
  fSelShowerRecoDownstreamX = kDefDoub;
  fSelShowerRecoDownstreamY = kDefDoub;
  fSelShowerRecoDownstreamZ = kDefDoub;
  fSelShowerRecoDownstreamT = kDefDoub;
  fSelShowerRecoEndClosestToVertexX = kDefDoub;
  fSelShowerRecoEndClosestToVertexY = kDefDoub;
  fSelShowerRecoEndClosestToVertexZ = kDefDoub;
  fSelShowerRecoContained = 0;
  fSelShowerRecoCharge = kDefDoub;
  fSelShowerRecoMomMCS = kDefDoub;
  fSelShowerRecoMomContained = kDefDoub;
  fSelShowerRecoNMatchedVertices = kDefInt;
  fSelShowerRecoVertexX = kDefDoub;
  fSelShowerRecoVertexY = kDefDoub;
  fSelShowerRecoVertexZ = kDefDoub;
  fSelShowerRecoNChildPFP = kDefInt;
  fSelShowerRecoNChildTrackPFP = kDefInt;
  fSelShowerRecoNChildShowerPFP = kDefInt;
  fSelShowerRecoDirX = kDefDoub;
  fSelShowerRecoDirY = kDefDoub;
  fSelShowerRecoDirZ = kDefDoub;
  fSelShowerRecoBestPlane = kDefInt;
  fSelShowerRecoLength = kDefDoub;
  fSelShowerRecoOpeningAngle = kDefDoub;
  fSelShowerRecoIsPrimaryPFPDaughter = 0;
  //MVA bits
  fSelShowerMVAElectron = kDefDoub;
  fSelShowerMVAPion = kDefDoub;
  fSelShowerMVAMuon = kDefDoub;
  fSelShowerMVAProton = kDefDoub;
  fSelShowerMVAPhoton = kDefDoub;
  fSelShowerPandrizzleEvalRatio     = kDefDoub;
  fSelShowerPandrizzleConcentration = kDefDoub;
  fSelShowerPandrizzleCoreHaloRatio = kDefDoub;
  fSelShowerPandrizzleConicalness   = kDefDoub;
  fSelShowerPandrizzledEdxBestPlane = kDefDoub;
  fSelShowerPandrizzleDisplacement  = kDefDoub;
  fSelShowerPandrizzleDCA           = kDefDoub;
  fSelShowerPandrizzleWideness      = kDefDoub;
  fSelShowerPandrizzleEnergyDensity = kDefDoub;
  fSelShowerPandrizzleBDTMethod = kDefDoub;
  fSelShowerJamPandrizzleMVAScore = kDefDoub;
  fSelShowerPandrizzleMVAScore      = kDefDoub;

  fSelShowerPandrizzleIsFilled      = 0;
  fSelShowerDistanceToNuVertexU     = kDefDoub;
  fSelShowerDistanceToNuVertexV     = kDefDoub;
  fSelShowerDistanceToNuVertexW     = kDefDoub;

  fSelShowerPandrizzleConnectionBDTScore = kDefDoub;
  fSelShowerPandrizzlePathwayLengthMin = kDefDoub;
  fSelShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D = kDefDoub;
  fSelShowerPandrizzleMaxNPostShowerStartHits = kDefDoub;
  fSelShowerPandrizzleMaxPostShowerStartScatterAngle = kDefDoub;
  fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry = kDefDoub;
  fSelShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry = kDefDoub;
  fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = kDefDoub;
  fSelShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius = kDefDoub;
  fSelShowerPandrizzleMaxPostShowerStartOpeningAngle = kDefDoub;
  fSelShowerPandrizzleMaxFoundHitRatio = kDefDoub;
  fSelShowerPandrizzleMaxInitialGapSize = kDefDoub;
  fSelShowerPandrizzleMinLargestProjectedGapSize = kDefDoub;
  fSelShowerPandrizzleNViewsWithAmbiguousHits = kDefDoub;
  fSelShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy = kDefDoub;

  for (int i_plane = 0; i_plane < 3; i_plane++) 
  {
      fSelShowerRecodEdx[i_plane] = kDefDoub;
      fSelShowerRecoEnergy[i_plane] = kDefDoub;
  }
  //All shower bits
  for (int i_shower = 0; i_shower < (fNRecoShowers == 0 ? kDefMaxNRecoShowers : fNRecoShowers); i_shower++){
    fRecoShowerTruePDG[i_shower] = kDefInt;
    fRecoShowerTruePrimary[i_shower] = false;
    fRecoShowerTrueMomX[i_shower] = kDefDoub;
    fRecoShowerTrueMomY[i_shower] = kDefDoub;
    fRecoShowerTrueMomZ[i_shower] = kDefDoub;
    fRecoShowerTrueMomT[i_shower] = kDefDoub;
    fRecoShowerTrueStartX[i_shower] = kDefDoub;
    fRecoShowerTrueStartY[i_shower] = kDefDoub;
    fRecoShowerTrueStartZ[i_shower] = kDefDoub;
    fRecoShowerTrueStartT[i_shower] = kDefDoub;
    fRecoShowerTrueEndX[i_shower] = kDefDoub;
    fRecoShowerTrueEndY[i_shower] = kDefDoub;
    fRecoShowerTrueEndZ[i_shower] = kDefDoub;
    fRecoShowerTrueEndT[i_shower] = kDefDoub;
    fRecoShowerRecoNHits[i_shower] = kDefInt;
    fRecoShowerActiveCheatElectronAlg[i_shower] = kDefInt;
    fRecoShowerActiveHybridElectronAlg[i_shower] = kDefInt;
    fRecoShowerActiveCheatGammaAlg[i_shower] = kDefInt;
    fRecoShowerActiveHybridGammaAlg[i_shower] = kDefInt;
    fRecoShowerRecoCompleteness[i_shower] = kDefDoub;
    fRecoShowerRecoHitPurity[i_shower] = kDefDoub;
    fRecoShowerRecoMomX[i_shower] = kDefDoub;
    fRecoShowerRecoMomY[i_shower] = kDefDoub;
    fRecoShowerRecoMomZ[i_shower] = kDefDoub;
    fRecoShowerRecoMomT[i_shower] = kDefDoub;
    fRecoShowerRecoStartX[i_shower] = kDefDoub;
    fRecoShowerRecoStartY[i_shower] = kDefDoub;
    fRecoShowerRecoStartZ[i_shower] = kDefDoub;
    fRecoShowerRecoStartT[i_shower] = kDefDoub;
    fRecoShowerRecoEndX[i_shower] = kDefDoub;
    fRecoShowerRecoEndY[i_shower] = kDefDoub;
    fRecoShowerRecoEndZ[i_shower] = kDefDoub;
    fRecoShowerRecoEndT[i_shower] = kDefDoub;
    fRecoShowerRecoUpstreamX[i_shower] = kDefDoub;
    fRecoShowerRecoUpstreamY[i_shower] = kDefDoub;
    fRecoShowerRecoUpstreamZ[i_shower] = kDefDoub;
    fRecoShowerRecoUpstreamT[i_shower] = kDefDoub;
    fRecoShowerRecoDownstreamX[i_shower] = kDefDoub;
    fRecoShowerRecoDownstreamY[i_shower] = kDefDoub;
    fRecoShowerRecoDownstreamZ[i_shower] = kDefDoub;
    fRecoShowerRecoDownstreamT[i_shower] = kDefDoub;
    fRecoShowerRecoEndClosestToVertexX[i_shower] = kDefDoub;
    fRecoShowerRecoEndClosestToVertexY[i_shower] = kDefDoub;
    fRecoShowerRecoEndClosestToVertexZ[i_shower] = kDefDoub;
    fRecoShowerRecoContained[i_shower] = kDefDoub;
    fRecoShowerRecoCharge[i_shower] = kDefDoub;
    fRecoShowerRecoMomMCS[i_shower] = kDefDoub;
    fRecoShowerRecoMomContained[i_shower] = kDefDoub;
    fRecoShowerRecoNMatchedVertices[i_shower] = kDefInt;
    fRecoShowerRecoVertexX[i_shower] = kDefDoub;
    fRecoShowerRecoVertexY[i_shower] = kDefDoub;
    fRecoShowerRecoVertexZ[i_shower] = kDefDoub;
    fRecoShowerRecoNChildPFP[i_shower] = kDefInt;
    fRecoShowerRecoNChildTrackPFP[i_shower] = kDefInt;
    fRecoShowerRecoNChildShowerPFP[i_shower] = kDefInt;
    fRecoShowerRecoDirX[i_shower] = kDefDoub;
    fRecoShowerRecoDirY[i_shower] = kDefDoub;
    fRecoShowerRecoDirZ[i_shower] = kDefDoub;
    fRecoShowerRecoBestPlane[i_shower] = kDefInt;
    fRecoShowerRecoLength[i_shower] = kDefDoub;
    fRecoShowerRecoOpeningAngle[i_shower] = kDefDoub;
    fRecoShowerRecoIsPrimaryPFPDaughter[i_shower] = false;
    fRecoShowerMVAElectron[i_shower] = kDefDoub;
    fRecoShowerMVAPion[i_shower] = kDefDoub;
    fRecoShowerMVAMuon[i_shower] = kDefDoub;
    fRecoShowerMVAProton[i_shower] = kDefDoub;
    fRecoShowerMVAPhoton[i_shower] = kDefDoub;
    for (int i_plane = 0; i_plane < 3; i_plane++){
        fRecoShowerRecodEdx[i_shower][i_plane] = kDefDoub;
        fRecoShowerRecoEnergy[i_shower][i_plane] = kDefDoub;
    }

    fRecoShowerPandrizzleConnectionBDTScore[i_shower] = kDefDoub;
    fRecoShowerPandrizzlePathwayLengthMin[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxNPostShowerStartHits[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxFoundHitRatio[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMaxInitialGapSize[i_shower] = kDefDoub;
    fRecoShowerPandrizzleMinLargestProjectedGapSize[i_shower] = kDefDoub;
    fRecoShowerPandrizzleNViewsWithAmbiguousHits[i_shower] = kDefDoub;
    fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[i_shower] = kDefDoub;

    fRecoShowerPandrizzleEvalRatio[i_shower]     = kDefDoub;
    fRecoShowerPandrizzleConcentration[i_shower] = kDefDoub;
    fRecoShowerPandrizzleCoreHaloRatio[i_shower] = kDefDoub;
    fRecoShowerPandrizzleConicalness[i_shower]   = kDefDoub;
    fRecoShowerPandrizzledEdxBestPlane[i_shower] = kDefDoub;
    fRecoShowerPandrizzleDisplacement[i_shower]  = kDefDoub;
    fRecoShowerPandrizzleDCA[i_shower]           = kDefDoub;
    fRecoShowerPandrizzleWideness[i_shower]      = kDefDoub;
    fRecoShowerPandrizzleEnergyDensity[i_shower] = kDefDoub;
    fRecoShowerPandrizzleBDTMethod[i_shower]     = kDefDoub;
    fRecoShowerJamPandrizzleMVAScore[i_shower]   = kDefDoub;
    fRecoShowerPandrizzleMVAScore[i_shower]      = kDefDoub;
    fRecoShowerPandrizzleIsFilled[i_shower]      = 0;
    fRecoShowerDistanceToNuVertexU[i_shower]     = kDefDoub;
    fRecoShowerDistanceToNuVertexV[i_shower]     = kDefDoub;
    fRecoShowerDistanceToNuVertexW[i_shower]     = kDefDoub;
  }
  fNRecoShowers = 0;
  //Reco nu bits
  fRecoNuVtxX = kDefDoub;
  fRecoNuVtxY = kDefDoub;
  fRecoNuVtxZ = kDefDoub;
  fRecoNuVtxNShowers = 0;
  fRecoNuVtxNTracks = 0;
  fRecoNuVtxNChildren = 0;



  //Reco energy bits
  fNueRecoEHad = kDefDoub;
  fNueRecoMomLep = kDefDoub;
  fNueRecoENu = kDefDoub; //Neutrino reco energy

  //Event level stuff
  fRecoEventCharge = kDefDoub;

}

//////////////////////////////////////////////

void FDSelection::CCNuSelection::GetEventInfo(art::Event const & evt)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

  // T0
  fT0 = trigger_offset(clockData);

  // Get total event charge
  std::vector<art::Ptr<recob::Hit>> hitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
  fRecoEventCharge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, hitList); 

  // Get CVN results
  art::Handle<std::vector<cvn::Result>> cvnResult;
  evt.getByLabel(fCVNModuleLabel, cvnResult);

  if(!cvnResult->empty()) 
  {
      fCVNResultNue = (*cvnResult)[0].GetNueProbability();
      fCVNResultNumu = (*cvnResult)[0].GetNumuProbability();
      fCVNResultNutau = (*cvnResult)[0].GetNutauProbability();
      fCVNResultNC = (*cvnResult)[0].GetNCProbability();
  }

  return;
}

//////////////////////////////////////////////

void FDSelection::CCNuSelection::GetTruthInfo(art::Event const & evt)
{
  //Get the generator record
  art::Handle< std::vector<simb::MCTruth> > mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle)){
    art::fill_ptr_vector(mcList, mcTruthListHandle);
  }

  if (mcList.size() > 1)
    mf::LogWarning("CCNuSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one.";

  //Get the flux record
  art::Handle< std::vector<simb::MCFlux> > mcFluxListHandle;
  std::vector<art::Ptr<simb::MCFlux> > mcFlux;
  if (evt.getByLabel(fNuGenModuleLabel, mcFluxListHandle)){
    art::fill_ptr_vector(mcFlux, mcFluxListHandle);
  }

  //need the assns for later
  art::FindManyP<simb::MCParticle> fmpt(mcTruthListHandle, evt, fLargeantModuleLabel);

  art::Ptr<simb::MCTruth> mcTruth = mcList.at(0);

  if (mcTruth->Origin() != simb::kBeamNeutrino) 
  {
    mf::LogWarning("CCNuSelection") << "Origin for this event is " << mcTruth->Origin() << " and not simb::kBeamNeutrino (" << simb::kBeamNeutrino<<")";
    return;
  }
 
  // Neutrino
  fNuPdg = mcTruth->GetNeutrino().Nu().PdgCode();

  if (mcFluxListHandle.isValid()) 
    fBeamPdg  = mcFlux[0]->fntype;

  fNC = mcTruth->GetNeutrino().CCNC();
  fMode = mcTruth->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  fTargetZ = mcTruth->GetNeutrino().Target()%100000000/10000;
  fENu = mcTruth->GetNeutrino().Nu().E();
  fQ2 = mcTruth->GetNeutrino().QSqr();
  fW = mcTruth->GetNeutrino().W();
  fX = mcTruth->GetNeutrino().X();
  fY = mcTruth->GetNeutrino().Y();
  fNuX = mcTruth->GetNeutrino().Nu().Vx();
  fNuY = mcTruth->GetNeutrino().Nu().Vy();
  fNuZ = mcTruth->GetNeutrino().Nu().Vz();
  fNuT = mcTruth->GetNeutrino().Nu().T();
  fNuMomX = mcTruth->GetNeutrino().Nu().Momentum().X();
  fNuMomY = mcTruth->GetNeutrino().Nu().Momentum().Y();
  fNuMomZ = mcTruth->GetNeutrino().Nu().Momentum().Z();
  fNuMomT = mcTruth->GetNeutrino().Nu().Momentum().T();

  //Leading lepton
  fLepPDG = mcTruth->GetNeutrino().Lepton().PdgCode();
  fMomLepX = mcTruth->GetNeutrino().Lepton().Momentum().X();
  fMomLepY = mcTruth->GetNeutrino().Lepton().Momentum().Y();
  fMomLepZ = mcTruth->GetNeutrino().Lepton().Momentum().Z();
  fMomLepT = mcTruth->GetNeutrino().Lepton().Momentum().T();
  fLepEndX = mcTruth->GetNeutrino().Lepton().EndPosition().X();
  fLepEndY = mcTruth->GetNeutrino().Lepton().EndPosition().Y();
  fLepEndZ = mcTruth->GetNeutrino().Lepton().EndPosition().Z();
  fLepEndY = mcTruth->GetNeutrino().Lepton().EndPosition().T();
  fLepNuAngle = mcTruth->GetNeutrino().Nu().Momentum().Vect().Angle(mcTruth->GetNeutrino().Lepton().Momentum().Vect());

  //Vertex particles
  const std::vector<art::Ptr<simb::MCParticle> > associated_particles = fmpt.at(mcTruth.key());

  fNVertexParticles = 0;

  for (unsigned int i_part = 0; i_part < associated_particles.size(); i_part++)
  {
    art::Ptr<simb::MCParticle> particle = associated_particles[i_part];
    if (particle->StatusCode() != 1) continue; //count tracked particles
    if (particle->Mother() != 0) continue; //count primary particles
    fVertexParticleIsGHEP[fNVertexParticles] = 0;
    fVertexParticlePDG[fNVertexParticles] = particle->PdgCode();;
    fVertexParticleStatus[fNVertexParticles] = particle->StatusCode();
    fVertexParticleNChildren[fNVertexParticles] = particle->NumberDaughters();
    fVertexParticleMomX[fNVertexParticles] = particle->Momentum(0).X();
    fVertexParticleMomY[fNVertexParticles] = particle->Momentum(0).Y();
    fVertexParticleMomZ[fNVertexParticles] = particle->Momentum(0).Z();
    fVertexParticleMomT[fNVertexParticles] = particle->Momentum(0).T();
    fVertexParticleEndX[fNVertexParticles] = particle->EndPosition().X();
    fVertexParticleEndY[fNVertexParticles] = particle->EndPosition().Y();
    fVertexParticleEndZ[fNVertexParticles] = particle->EndPosition().Z();
    fVertexParticleEndT[fNVertexParticles] = particle->EndPosition().T();
    fNVertexParticles++;
  }

  //Loop over the final state particles from the ghep vertex
  for (int i_part = 0; i_part < mcTruth->NParticles(); i_part++)
  {
    const simb::MCParticle& vertex_particle = mcTruth->GetParticle(i_part);
    int pdg = vertex_particle.PdgCode();
    fVertexParticleIsGHEP[fNVertexParticles] = 1;
    fVertexParticlePDG[fNVertexParticles] = pdg;
    fVertexParticleStatus[fNVertexParticles] = vertex_particle.StatusCode();
    fVertexParticleNChildren[fNVertexParticles] = vertex_particle.NumberDaughters();
    fVertexParticleMomX[fNVertexParticles] = vertex_particle.Momentum(0).X();
    fVertexParticleMomY[fNVertexParticles] = vertex_particle.Momentum(0).Y();
    fVertexParticleMomZ[fNVertexParticles] = vertex_particle.Momentum(0).Z();
    fVertexParticleMomT[fNVertexParticles] = vertex_particle.Momentum(0).T();
    fVertexParticleEndX[fNVertexParticles] = vertex_particle.EndPosition().X();
    fVertexParticleEndY[fNVertexParticles] = vertex_particle.EndPosition().Y();
    fVertexParticleEndZ[fNVertexParticles] = vertex_particle.EndPosition().Z();
    fVertexParticleEndT[fNVertexParticles] = vertex_particle.EndPosition().T();
    fNVertexParticles++;
    if (!(vertex_particle.StatusCode() == 1)) continue;
    if (pdg >= 2000000000) continue;
    if (std::abs(pdg) >= 11 && std::abs(pdg) <= 16) continue;
    if (pdg==211) fNPiP++;
    else if (pdg==-211) fNPim++;
    else if (pdg==111) fNPi0++;
    else if (pdg==22) fNPhotons++;
    else if (pdg==2212) fNProtons++;
    else if (pdg==2112) fNNeutrons++;
    else fNOther++;
  }

  if (fNVertexParticles > kDefMaxNTrueVertexParticles)
  {
    std::cout << "CCNuSelection::GetTruthInfo VERTEX ARRAY IS GOING TO BE OVERFILLED - VERY BAD" << std::endl;
    throw;
  }

  //do some transverse momentum stuff
  TVector3 beam_axis(0, 0, 1);
  beam_axis.RotateX(-0.101);
  TVector3 nu_mom_vect = mcTruth->GetNeutrino().Nu().Momentum().Vect();
  TVector3 nu_mom_tran_vect = ProjectVectorOntoPlane(nu_mom_vect, beam_axis);
  fNuMomTranMag = nu_mom_tran_vect.Mag(); 

  TVector3 total_initial_mom_tran_vect;
  total_initial_mom_tran_vect += nu_mom_tran_vect;

  TVector3 total_final_mom_tran_vect_nolep_norem;
  TVector3 total_final_mom_tran_vect_nolep_withrem;
  TVector3 total_final_mom_tran_vect_withlep_norem;
  TVector3 total_final_mom_tran_vect_withlep_withrem;
  
  //Loop over the particles
  for (int i_part = 0; i_part < mcTruth->NParticles(); i_part++){
    const simb::MCParticle& particle = mcTruth->GetParticle(i_part);
    int status = particle.StatusCode();
    if (status == 11){ //Nucleon target
      TVector3 target_nucleon_mom_vect = particle.Momentum(0).Vect();
      TVector3 target_nucleon_mom_tran_vect = ProjectVectorOntoPlane(target_nucleon_mom_vect, beam_axis); 
      total_initial_mom_tran_vect += target_nucleon_mom_tran_vect;
      fTargNuclMomTranMag = target_nucleon_mom_tran_vect.Mag();
    }
    else if (status == 15){ //final nuclear remnant
      TVector3 nuclear_rem_mom_vect = particle.Momentum(0).Vect();
      TVector3 nuclear_rem_mom_tran_vect = ProjectVectorOntoPlane(nuclear_rem_mom_vect, beam_axis);
      total_final_mom_tran_vect_nolep_withrem += nuclear_rem_mom_tran_vect;
      total_final_mom_tran_vect_withlep_withrem += nuclear_rem_mom_tran_vect;
      fNuclRemMomTranMag = nuclear_rem_mom_tran_vect.Mag();
    }
    else if (status == 1){ //final state particle
      TVector3 fin_state_mom_vect = particle.Momentum(0).Vect();
      TVector3 fin_state_mom_tran_vect = ProjectVectorOntoPlane(fin_state_mom_vect, beam_axis);
      int pdg = particle.PdgCode();
      if (std::abs(pdg) >= 11 && std::abs(pdg) <= 16){
        total_final_mom_tran_vect_withlep_withrem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_withlep_norem += fin_state_mom_tran_vect;
        fLepMomTranMag = fin_state_mom_tran_vect.Mag();
      }
      else {
        total_final_mom_tran_vect_withlep_withrem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_withlep_norem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_nolep_withrem += fin_state_mom_tran_vect;
        total_final_mom_tran_vect_nolep_norem += fin_state_mom_tran_vect;
      }
    }
  }

  //Get all of the mom mags now
  fInitalMomTranMag = total_initial_mom_tran_vect.Mag();
  fFinalMomTranMagNoLepNoRem = total_final_mom_tran_vect_nolep_norem.Mag();
  fFinalMomTranMagNoLepWithRem = total_final_mom_tran_vect_nolep_withrem.Mag();
  fFinalMomTranMagWithLepNoRem = total_final_mom_tran_vect_withlep_norem.Mag();
  fFinalMomTranMagWithLepWithRem = total_final_mom_tran_vect_withlep_withrem.Mag();
}

////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::GetRecoTrackInfo(art::Event const & evt)
{
  /*
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > trackList;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle)){
    art::fill_ptr_vector(trackList, trackListHandle);
  }

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  //Now make the entire PFPPArticleMap
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);
  

  //Get the track -> hits assns
  art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
  */

  // Get PFP list from event, loop through PFPs and get tracks.. sigh.. 


  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel);

  std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fPFParticleModuleLabel);

  int trackCounter = 0;

  for (art::Ptr<recob::PFParticle> pfp : pfps)
  {
    if (trackCounter == kDefMaxNRecoTracks)
      break;

    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
      continue;

    for (art::Ptr<recob::PFParticle> nuChild : nuChildren)
    {
      if (nuChild->Self() == pfp->Self())
      {
        fRecoTrackIsPrimary[trackCounter] = true;
        break;
      }
    }

    const std::vector<art::Ptr<recob::Hit>> current_track_hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPFParticleModuleLabel);
    fRecoTrackRecoNHits[trackCounter] = current_track_hits.size();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    fRecoTrackRecoCharge[trackCounter]  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, current_track_hits); 

    art::Ptr<recob::Track> current_track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);

    // TODO
    fRecoTrackRecoMomX[trackCounter] = kDefDoub; //temp
    fRecoTrackRecoMomY[trackCounter] = kDefDoub; //temp
    fRecoTrackRecoMomZ[trackCounter] = kDefDoub; //temp
    fRecoTrackRecoMomT[trackCounter] = kDefDoub; //temp

    // General track variables
    recob::Track::Point_t trackStart, trackEnd;
    std::tie(trackStart, trackEnd) = current_track->Extent(); 
    fRecoTrackRecoStartX[trackCounter] = trackStart.X();
    fRecoTrackRecoStartY[trackCounter] = trackStart.Y();
    fRecoTrackRecoStartZ[trackCounter] = trackStart.Z();
    fRecoTrackRecoEndX[trackCounter] = trackEnd.X();
    fRecoTrackRecoEndY[trackCounter] = trackEnd.Y();
    fRecoTrackRecoEndZ[trackCounter] = trackEnd.Z();
    if (fRecoTrackRecoEndZ[trackCounter] > fRecoTrackRecoStartZ[trackCounter]){
      fRecoTrackRecoUpstreamX[trackCounter] = fRecoTrackRecoStartX[trackCounter];
      fRecoTrackRecoUpstreamY[trackCounter] = fRecoTrackRecoStartY[trackCounter];
      fRecoTrackRecoUpstreamZ[trackCounter] = fRecoTrackRecoStartZ[trackCounter];
      fRecoTrackRecoDownstreamX[trackCounter] = fRecoTrackRecoEndX[trackCounter];
      fRecoTrackRecoDownstreamY[trackCounter] = fRecoTrackRecoEndY[trackCounter];
      fRecoTrackRecoDownstreamZ[trackCounter] = fRecoTrackRecoEndZ[trackCounter];
    }
    else{
      fRecoTrackRecoDownstreamX[trackCounter] = fRecoTrackRecoStartX[trackCounter];
      fRecoTrackRecoDownstreamY[trackCounter] = fRecoTrackRecoStartY[trackCounter];
      fRecoTrackRecoDownstreamZ[trackCounter] = fRecoTrackRecoStartZ[trackCounter];
      fRecoTrackRecoUpstreamX[trackCounter] = fRecoTrackRecoEndX[trackCounter];
      fRecoTrackRecoUpstreamY[trackCounter] = fRecoTrackRecoEndY[trackCounter];
      fRecoTrackRecoUpstreamZ[trackCounter] = fRecoTrackRecoEndZ[trackCounter];
    }
    fRecoTrackRecoLength[trackCounter] = current_track->Length();

    TVector3 upstream_end(fRecoTrackRecoUpstreamX[trackCounter],fRecoTrackRecoUpstreamY[trackCounter],fRecoTrackRecoUpstreamZ[trackCounter]);
    TVector3 downstream_end(fRecoTrackRecoDownstreamX[trackCounter],fRecoTrackRecoDownstreamY[trackCounter],fRecoTrackRecoDownstreamZ[trackCounter]);
    TVector3 vertex_pos(fRecoTrackRecoVertexX[trackCounter],fRecoTrackRecoVertexY[trackCounter],fRecoTrackRecoVertexZ[trackCounter]);
    if ((vertex_pos-upstream_end).Mag() < (vertex_pos-downstream_end).Mag()){
      fRecoTrackRecoEndClosestToVertexX[trackCounter] = fRecoTrackRecoUpstreamX[trackCounter];
      fRecoTrackRecoEndClosestToVertexY[trackCounter] = fRecoTrackRecoUpstreamY[trackCounter];
      fRecoTrackRecoEndClosestToVertexZ[trackCounter] = fRecoTrackRecoUpstreamZ[trackCounter];
    }
    else{
      fRecoTrackRecoEndClosestToVertexX[trackCounter] = fRecoTrackRecoDownstreamX[trackCounter];
      fRecoTrackRecoEndClosestToVertexY[trackCounter] = fRecoTrackRecoDownstreamY[trackCounter];
      fRecoTrackRecoEndClosestToVertexZ[trackCounter] = fRecoTrackRecoDownstreamZ[trackCounter];
    }

    // Child particle info
    FillChildPFPInformation(pfp, evt, fRecoTrackRecoNChildPFP[trackCounter], fRecoTrackRecoNChildTrackPFP[trackCounter], fRecoTrackRecoNChildShowerPFP[trackCounter]);

    std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, current_track_hits, 1);
    fRecoTrackRecoCompleteness[trackCounter] = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, current_track_hits, eventHitList, g4id);
    fRecoTrackRecoHitPurity[trackCounter] = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, current_track_hits, g4id);

    if (TruthMatchUtils::Valid(g4id)){
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
        if (matched_mcparticle){
          fRecoTrackTruePDG[trackCounter] = matched_mcparticle->PdgCode();
          if (matched_mcparticle->Mother() == 0) fRecoTrackTruePrimary[trackCounter] = true;
          else fRecoTrackTruePrimary[trackCounter] = false;
          fRecoTrackTrueMomX[trackCounter] = matched_mcparticle->Momentum().X();
          fRecoTrackTrueMomY[trackCounter] = matched_mcparticle->Momentum().Y();
          fRecoTrackTrueMomZ[trackCounter] = matched_mcparticle->Momentum().Z();
          fRecoTrackTrueMomT[trackCounter] = matched_mcparticle->Momentum().T();
          fRecoTrackTrueStartX[trackCounter] = matched_mcparticle->Position(0).X();
          fRecoTrackTrueStartY[trackCounter] = matched_mcparticle->Position(0).Y();
          fRecoTrackTrueStartZ[trackCounter] = matched_mcparticle->Position(0).Z();
          fRecoTrackTrueStartT[trackCounter] = matched_mcparticle->Position(0).T();
          fRecoTrackTrueEndX[trackCounter] = matched_mcparticle->EndPosition().X();
          fRecoTrackTrueEndY[trackCounter] = matched_mcparticle->EndPosition().Y();
          fRecoTrackTrueEndZ[trackCounter] = matched_mcparticle->EndPosition().Z();
          fRecoTrackTrueEndT[trackCounter] = matched_mcparticle->EndPosition().T();
        }
    }

    // MVA PID
    art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Track> >{current_track}, evt, fPIDModuleLabel);
    art::Ptr<anab::MVAPIDResult> pid(findPIDResult.at(0));

    if (pid.isAvailable())
    {
      std::map<std::string,double> mvaOutMap = pid->mvaOutput;
      if (!mvaOutMap.empty())
      {
        fRecoTrackMVAElectron[trackCounter] = mvaOutMap["electron"];
        fRecoTrackMVAPion[trackCounter] = mvaOutMap["pich"];
        fRecoTrackMVAMuon[trackCounter] = mvaOutMap["muon"];
        fRecoTrackMVAProton[trackCounter] = mvaOutMap["proton"];
        fRecoTrackMVAPhoton[trackCounter] = mvaOutMap["photon"];
      }
    }

    // Pandizzle variables
    FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(current_track, evt));
    fRecoTrackMichelNHits[trackCounter] = (float)pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelNHits);
    fRecoTrackMichelElectronMVA[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelElectronMVA);
    fRecoTrackMichelRecoEnergyPlane2[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelRecoEnergyPlane2);
    fRecoTrackDeflecAngleSD[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD);
    fRecoTrackLength[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength);
    fRecoTrackEvalRatio[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio);
    fRecoTrackConcentration[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration);
    fRecoTrackCoreHaloRatio[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio);
    fRecoTrackConicalness[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness);
    fRecoTrackdEdxStart[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart);
    fRecoTrackdEdxEnd[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd);
    fRecoTrackdEdxEndRatio[trackCounter] = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio);
    fRecoTrackPandizzleVar[trackCounter] = pandizzleRecord.GetMVAScore();

    trackCounter++;
  }

  fNRecoTracks = trackCounter;

  return;
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunTrackSelection(art::Event const & evt){
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(evt.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return;
  }

  //Get the selected track
  art::Ptr<recob::Track> sel_track = fRecoTrackSelector->FindSelectedTrack(evt);

  //If we didn't find a selected track then what's the point?
  if (!(sel_track.isAvailable())) {
    std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no track returned from selection" << std::endl; 
    return;
  }

  //24/07/18 DBrailsford Get the reco energy data product for neutrinos
  //art::Handle<dune::EnergyRecoOutput> energyRecoHandle;
  //if (!evt.getByLabel(fNumuEnergyRecoModuleLabel, energyRecoHandle)) {
  //  std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no energy reconstruction found with label " << fNumuEnergyRecoModuleLabel << std::endl;
  //  return;
  //}

  //Get the hits for said track
  art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > sel_track_hits = fmht.at(sel_track.key());
  fSelTrackRecoNHits = sel_track_hits.size();

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  //Start filling some variables
  recob::Track::Point_t trackStart, trackEnd;
  std::tie(trackStart, trackEnd) = sel_track->Extent(); 
  fSelTrackRecoMomX = kDefDoub; //temp
  fSelTrackRecoMomY = kDefDoub; //temp
  fSelTrackRecoMomZ = kDefDoub; //temp
  fSelTrackRecoMomT = kDefDoub; //temp
  fSelTrackRecoStartX = trackStart.X();
  fSelTrackRecoStartY = trackStart.Y();
  fSelTrackRecoStartZ = trackStart.Z();
  fSelTrackRecoEndX = trackEnd.X();
  fSelTrackRecoEndY = trackEnd.Y();
  fSelTrackRecoEndZ = trackEnd.Z();
  if (fSelTrackRecoEndZ > fSelTrackRecoStartZ){
    fSelTrackRecoUpstreamX = fSelTrackRecoStartX;
    fSelTrackRecoUpstreamY = fSelTrackRecoStartY;
    fSelTrackRecoUpstreamZ = fSelTrackRecoStartZ;
    fSelTrackRecoDownstreamX = fSelTrackRecoEndX;
    fSelTrackRecoDownstreamY = fSelTrackRecoEndY;
    fSelTrackRecoDownstreamZ = fSelTrackRecoEndZ;
  }
  else{
    fSelTrackRecoDownstreamX = fSelTrackRecoStartX;
    fSelTrackRecoDownstreamY = fSelTrackRecoStartY;
    fSelTrackRecoDownstreamZ = fSelTrackRecoStartZ;
    fSelTrackRecoUpstreamX = fSelTrackRecoEndX;
    fSelTrackRecoUpstreamY = fSelTrackRecoEndY;
    fSelTrackRecoUpstreamZ = fSelTrackRecoEndZ;
  }
  fSelTrackRecoLength = sel_track->Length();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  fSelTrackRecoCharge  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, sel_track_hits); 

  //Now that the vertex information has been filled.  Calculate which end of th etrack is closest
  TVector3 upstream_end(fSelTrackRecoUpstreamX,fSelTrackRecoUpstreamY,fSelTrackRecoUpstreamZ);
  TVector3 downstream_end(fSelTrackRecoDownstreamX,fSelTrackRecoDownstreamY,fSelTrackRecoDownstreamZ);
  TVector3 vertex_pos(fSelTrackRecoVertexX,fSelTrackRecoVertexY,fSelTrackRecoVertexZ);
  if ((vertex_pos-upstream_end).Mag() < (vertex_pos-downstream_end).Mag()){
    fSelTrackRecoEndClosestToVertexX = fSelTrackRecoUpstreamX;
    fSelTrackRecoEndClosestToVertexY = fSelTrackRecoUpstreamY;
    fSelTrackRecoEndClosestToVertexZ = fSelTrackRecoUpstreamZ;
  }
  else{
    fSelTrackRecoEndClosestToVertexX = fSelTrackRecoDownstreamX;
    fSelTrackRecoEndClosestToVertexY = fSelTrackRecoDownstreamY;
    fSelTrackRecoEndClosestToVertexZ = fSelTrackRecoDownstreamZ;
  }

  //28/11/18 DBrailsford Fill child PFP info
  FillChildPFPInformation(sel_track, evt, fSelTrackRecoNChildPFP, fSelTrackRecoNChildTrackPFP, fSelTrackRecoNChildShowerPFP);
  //24/07/18 DBrailsford Use the data product to get the neutrino energy
  //fSelTrackRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);
  std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(sel_track, evt)));
  fNumuRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNumuRecoEHad = energyRecoHandle->fHadLorentzVector.E();

  fSelTrackRecoContained = energyRecoHandle->longestTrackContained; 
  fSelTrackRecoMomMethod = energyRecoHandle->trackMomMethod;

  if (energyRecoHandle->trackMomMethod==1){ //momentum by range was used to calculate ENu
    fSelTrackRecoMomContained = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fNumuRecoMomLep = fSelTrackRecoMomContained;
  }
  else if (energyRecoHandle->trackMomMethod==0){//momentum by MCS
    fSelTrackRecoMomMCS = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fNumuRecoMomLep = fSelTrackRecoMomMCS;
  }

  int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sel_track_hits, 1);
  fSelTrackRecoCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, sel_track_hits, hitList, g4id);
  fSelTrackRecoHitPurity = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, sel_track_hits, g4id);

  if (TruthMatchUtils::Valid(g4id))
  {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
      if (matched_mcparticle){
        //Fill variables
        fSelTrackTruePDG = matched_mcparticle->PdgCode();
        if (matched_mcparticle->Mother()==0) fSelTrackTruePrimary = 1;
        else fSelTrackTruePrimary = 0;
        fSelTrackTrueMomX = matched_mcparticle->Momentum().X();
        fSelTrackTrueMomY = matched_mcparticle->Momentum().Y();
        fSelTrackTrueMomZ = matched_mcparticle->Momentum().Z();
        fSelTrackTrueMomT = matched_mcparticle->Momentum().T();
        fSelTrackTrueStartX = matched_mcparticle->Position(0).X();
        fSelTrackTrueStartY = matched_mcparticle->Position(0).Y();
        fSelTrackTrueStartZ = matched_mcparticle->Position(0).Z();
        fSelTrackTrueStartT = matched_mcparticle->Position(0).T();
        fSelTrackTrueEndX = matched_mcparticle->EndPosition().X();
        fSelTrackTrueEndY = matched_mcparticle->EndPosition().Y();
        fSelTrackTrueEndZ = matched_mcparticle->EndPosition().Z();
        fSelTrackTrueEndT = matched_mcparticle->EndPosition().T();
      }
  }

  //Now get the pid stuff
  art::FindManyP<anab::MVAPIDResult> fmpidt(trackListHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(sel_track.key());
  if (pids.at(0).isAvailable()){
      std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
      if (!(mvaOutMap.empty())){
          //Get the PIDs
          fSelTrackMVAElectron = mvaOutMap["electron"];
          fSelTrackMVAPion = mvaOutMap["pich"];
          fSelTrackMVAMuon = mvaOutMap["muon"];
          fSelTrackMVAProton = mvaOutMap["proton"];
          fSelTrackMVAPhoton = mvaOutMap["photon"];
      }
  }

  //Get the pandizzle variables
  art::Ptr<recob::PFParticle> track_pfp = GetPFParticleMatchedToTrack(sel_track, evt);
  FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(sel_track, evt));

  fSelTrackMichelNHits = (float)pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelNHits);
  fSelTrackMichelElectronMVA = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelElectronMVA);
  fSelTrackMichelRecoEnergyPlane2 = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelRecoEnergyPlane2);
  fSelTrackDeflecAngleSD = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD);
  fSelTrackLength = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength);
  fSelTrackEvalRatio = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio);
  fSelTrackConcentration = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration);
  fSelTrackCoreHaloRatio = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio);
  fSelTrackConicalness = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness);
  fSelTrackdEdxStart = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart);
  fSelTrackdEdxEnd = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd);
  fSelTrackdEdxEndRatio = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio);
  fSelTrackPandizzleVar = pandizzleRecord.GetMVAScore();

  //20/04/20 DBrailsford
  //Get the DeepPan variables
  ctp::CTPResult deepPanPIDResult = fConvTrackPID.RunConvolutionalTrackPID(track_pfp, evt);
  if (deepPanPIDResult.IsValid()){
      fSelTrackDeepPanMuVar = deepPanPIDResult.GetMuonScore();
      fSelTrackDeepPanPiVar = deepPanPIDResult.GetPionScore();
      fSelTrackDeepPanProtonVar = deepPanPIDResult.GetProtonScore();
  }
}

//double FDSelection::CCNuSelection::CalculateTrackCharge(art::Ptr<recob::Track> const track, std::vector< art::Ptr< recob::Hit> > const track_hits){
//  double charge = 0;
//  for (unsigned int i_hit = 0; i_hit < track_hits.size(); i_hit++){
//    if (track_hits[i_hit]->WireID().Plane != 2) continue; //Collection only
//    charge += track_hits[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(track_hits[i_hit]->PeakTime(), fT0);
//  }
//  return charge;
//}


bool FDSelection::CCNuSelection::IsTrackContained(art::Ptr<recob::Track> const track, std::vector< art::Ptr<recob::Hit > > const track_hits, art::Event const & evt){
  //Get the space points for each of the hits
  //Annoyingly we have to go from the start of the handle for the hits...
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }
  art::FindManyP<recob::SpacePoint> fmhs(hitListHandle, evt, fTrackModuleLabel);
  for (unsigned int i_hit = 0; i_hit < track_hits.size(); i_hit++){
    if (track_hits[i_hit]->WireID().Plane != 2) continue;
    std::vector<art::Ptr<recob::SpacePoint> > space_points = fmhs.at(track_hits[i_hit].key());
    if (space_points.size()){
      //Make a TVector3
      TVector3 position(space_points[0]->XYZ()[0],space_points[0]->XYZ()[1],space_points[0]->XYZ()[2]);
      bool is_in_tpc = FDSelectionUtils::IsInsideTPC(position,20); //20cm buffer from the wall
      if (!is_in_tpc){
        return false;
      }
    }
  }
  return true;
}

void FDSelection::CCNuSelection::FillVertexInformation(art::Event const & evt){
  //Have to get the PFP and then get the VTX
  //Need the PFP vector handle
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);



  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);
  //Loop over the neutrinos
  if (nu_pfps.size() != 1){
    std::cout<<"FDSelection::CCNuSelection: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return; //do nothing
  }
  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  //Count the number of children
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);
  for (int i_childpfp = 0; i_childpfp < nu_pfp->NumDaughters(); i_childpfp++){
    int child_id = nu_pfp->Daughter(i_childpfp);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[child_id];
    int pdg = child_pfp->PdgCode();
    if (pdg == 11) fRecoNuVtxNShowers++;
    else if (pdg == 13) fRecoNuVtxNTracks++;
    fRecoNuVtxNChildren++;
  }

  art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());
  /*
  fSelTrackRecoNMatchedVertices = sel_pfp_vertices.size();
  if (fSelTrackRecoNMatchedVertices == 0){ //Nothing to do
    return;
  }
  else if (fSelTrackRecoNMatchedVertices > 1){
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << fSelTrackRecoNMatchedVertices << std::endl;
  }
  */
  if (sel_pfp_vertices.size() == 0){ //Nothing to do
    return;
  }
  else if (sel_pfp_vertices.size() > 1){
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  return;
}

void FDSelection::CCNuSelection::FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp)
{
  n_child_pfp = 0;
  n_child_track_pfp = 0;
  n_child_shower_pfp = 0;

  std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfp, evt, fPFParticleModuleLabel);

  for (art::Ptr<recob::PFParticle> childPFP : childPFPs)
  {
    int pdg = childPFP->PdgCode();
    if (pdg == 13) n_child_track_pfp++;
    else if (pdg == 11) n_child_shower_pfp++;
    else std::cout << "FillChildPFPInformation: found a child PFP with an unexpected pdg code: " << pdg << std::endl;
  }

  return;
}


void FDSelection::CCNuSelection::FillChildPFPInformation(art::Ptr<recob::Track> const track, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp){
  art::Ptr<recob::PFParticle> pfp = GetPFParticleMatchedToTrack(track, evt);
  FillChildPFPInformation(pfp, evt, n_child_pfp, n_child_track_pfp, n_child_shower_pfp);
  return;
}

void FDSelection::CCNuSelection::FillChildPFPInformation(art::Ptr<recob::Shower> const shower, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp){
  art::Ptr<recob::PFParticle> pfp = GetPFParticleMatchedToShower(shower, evt);
  FillChildPFPInformation(pfp, evt, n_child_pfp, n_child_track_pfp, n_child_shower_pfp);
  return;
}

art::Ptr<recob::PFParticle> FDSelection::CCNuSelection::GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt){
  art::Ptr<recob::PFParticle> matched_pfp;
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(evt.getByLabel(fTrackModuleLabel, trackListHandle))){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToTrack Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel << std::endl;
    return matched_pfp;
  }
  art::FindManyP<recob::PFParticle> fmpfpt(trackListHandle, evt, fTrackModuleLabel);
  const std::vector<art::Ptr<recob::PFParticle> > sel_track_pfps = fmpfpt.at(track.key());
  if (sel_track_pfps.size() != 1){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToTrack NUMBER OF PFP MATCHED TO A TRACK DOES NOT EQUAL 1: " << sel_track_pfps.size() << std::endl;
    return matched_pfp;
  }
  matched_pfp = sel_track_pfps[0];
  return matched_pfp;
}

art::Ptr<recob::PFParticle> FDSelection::CCNuSelection::GetPFParticleMatchedToShower(art::Ptr<recob::Shower> const shower, art::Event const & evt){
  art::Ptr<recob::PFParticle> matched_pfp;
  art::Handle< std::vector<recob::Shower> > showerListHandle;

  if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle))){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToShower Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
    return matched_pfp;
  }

  art::FindManyP<recob::PFParticle> fmpfpt(showerListHandle, evt, fShowerModuleLabel);
  const std::vector<art::Ptr<recob::PFParticle> > sel_shower_pfps = fmpfpt.at(shower.key());
  if (sel_shower_pfps.size() != 1){
    std::cout<<"CCNuSelection::GetPFParticleMatchedToShower NUMBER OF PFP MATCHED TO A SHOWER DOES NOT EQUAL 1: " << sel_shower_pfps.size() << std::endl;
    return matched_pfp;
  }

  matched_pfp = sel_shower_pfps[0];

  return matched_pfp;
}

art::Ptr<recob::Track> FDSelection::CCNuSelection::GetTrackMatchedPFParticle(art::Ptr<recob::PFParticle> const pfparticle, art::Event const & evt){
  art::Ptr<recob::Track> matched_track;
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return matched_track; //empty
  }
  art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Track> > sel_pfp_tracks = fmtpfp.at(pfparticle.key());
  if (sel_pfp_tracks.size() != 1){
    std::cout<<"CCNuSelection::GetTrackMatchedPFParticle number of tracks matched to PFP does not equal 1: " << sel_pfp_tracks.size() << std::endl;
    return matched_track; //empty
  }
  matched_track = sel_pfp_tracks[0];
  return matched_track;
}




void FDSelection::CCNuSelection::GetRecoShowerInfo(art::Event const & evt){
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerList;
  if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle))){
    std::cout<<"Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
    return;
  }
  else art::fill_ptr_vector(showerList, showerListHandle);

  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle > > pfparticleList;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))){
    std::cout<<"Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  else art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  //const std::vector<art::Ptr<recob::Hit> > sel_shower_hits = fmhs.at(sel_shower.key());
  //fSelShowerRecoNHits = sel_shower_hits.size();

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn(pfparticleListHandle, evt, "pandoraSel");

  fNRecoShowers = (showerList.size() > kDefMaxNRecoShowers ? kDefMaxNRecoShowers : showerList.size());

  const unsigned int loopLimit = fNRecoShowers;

  //Loop
  for (unsigned int i_shower = 0; i_shower < loopLimit; i_shower++){
    art::Ptr<recob::Shower> current_shower = showerList[i_shower];
    const std::vector<art::Ptr<recob::Hit> > current_shower_hits = fmhs.at(current_shower.key());

    fRecoShowerRecoNHits[i_shower] = current_shower_hits.size();

    art::Ptr<recob::PFParticle> current_shower_pfp = GetPFParticleMatchedToShower(current_shower, evt);

    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(current_shower_pfp.key());
    larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetadata[0]->GetPropertiesMap());

    fRecoShowerActiveHybridElectronAlg[i_shower] = (propertiesMap.find("ActiveHybridElectronAlg") != propertiesMap.end()) ? 1 : 0;
    fRecoShowerActiveCheatElectronAlg[i_shower] = (propertiesMap.find("ActiveCheatingElectronAlg") != propertiesMap.end()) ? 1 : 0;
    fRecoShowerActiveHybridGammaAlg[i_shower] = (propertiesMap.find("ActiveHybridGammaAlg") != propertiesMap.end()) ? 1 : 0;
    fRecoShowerActiveCheatGammaAlg[i_shower] = (propertiesMap.find("ActiveCheatingGammaAlg") != propertiesMap.end()) ? 1 : 0;

    //Is this PFParticle a primary daughter of the neutrino
    fRecoShowerRecoIsPrimaryPFPDaughter[i_shower] = IsPFParticlePrimaryDaughter(current_shower_pfp, pfparticleList);

    fRecoShowerRecoDirX[i_shower]   = current_shower->Direction().X();
    fRecoShowerRecoDirY[i_shower]   = current_shower->Direction().Y();
    fRecoShowerRecoDirZ[i_shower]   = current_shower->Direction().Z();
    fRecoShowerRecoStartX[i_shower] = current_shower->ShowerStart().X();
    fRecoShowerRecoStartY[i_shower] = current_shower->ShowerStart().Y();
    fRecoShowerRecoStartZ[i_shower] = current_shower->ShowerStart().Z();

    fRecoShowerRecoBestPlane[i_shower] = current_shower->best_plane();
    fRecoShowerRecoLength[i_shower] = current_shower->Length();
    fRecoShowerRecoOpeningAngle[i_shower] = current_shower->OpenAngle();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    fRecoShowerRecoCharge[i_shower]  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, current_shower_hits); 

    //24/07/18 DBrailsford Use the data product to get the neutrino energy
    //fRecoShowerRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);
    //std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(current_shower, evt)));
    //fNueRecoENu = energyRecoHandle->fNuLorentzVector.E();
    //fNueRecoEHad = energyRecoHandle->fHadLorentzVector.E();
    //fNueRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    FillChildPFPInformation(current_shower, evt, fRecoShowerRecoNChildPFP[i_shower], fRecoShowerRecoNChildTrackPFP[i_shower], fRecoShowerRecoNChildShowerPFP[i_shower]);

    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, current_shower_hits, 1);
    fRecoShowerRecoCompleteness[i_shower] = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, current_shower_hits, hitList, g4id);
    fRecoShowerRecoHitPurity[i_shower] = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, current_shower_hits, g4id);

    if (TruthMatchUtils::Valid(g4id))
    {
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->TrackIdToParticle_P(g4id);
        if (matched_mcparticle){
          //Fill variables
          fRecoShowerTruePDG[i_shower] = matched_mcparticle->PdgCode();
          if (matched_mcparticle->Mother()==0) fRecoShowerTruePrimary[i_shower] = 1;
          else fRecoShowerTruePrimary[i_shower] = 0;
          fRecoShowerTrueMomX[i_shower] = matched_mcparticle->Momentum().X();
          fRecoShowerTrueMomY[i_shower] = matched_mcparticle->Momentum().Y();
          fRecoShowerTrueMomZ[i_shower] = matched_mcparticle->Momentum().Z();
          fRecoShowerTrueMomT[i_shower] = matched_mcparticle->Momentum().T();
          fRecoShowerTrueStartX[i_shower] = matched_mcparticle->Position(0).X();
          fRecoShowerTrueStartY[i_shower] = matched_mcparticle->Position(0).Y();
          fRecoShowerTrueStartZ[i_shower] = matched_mcparticle->Position(0).Z();
          fRecoShowerTrueStartT[i_shower] = matched_mcparticle->Position(0).T();
          fRecoShowerTrueEndX[i_shower] = matched_mcparticle->EndPosition().X();
          fRecoShowerTrueEndY[i_shower] = matched_mcparticle->EndPosition().Y();
          fRecoShowerTrueEndZ[i_shower] = matched_mcparticle->EndPosition().Z();
          fRecoShowerTrueEndT[i_shower] = matched_mcparticle->EndPosition().T();
        }
    }

    //Now get the pid stuff
    art::FindManyP<anab::MVAPIDResult> fmpidt(showerListHandle, evt, fPIDModuleLabel);
    std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(current_shower.key());
    if (pids.at(0).isAvailable()){
        std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
        if (!(mvaOutMap.empty())){
            //Get the PIDs
            fRecoShowerMVAElectron[i_shower] = mvaOutMap["electron"];
            fRecoShowerMVAPion[i_shower] = mvaOutMap["pich"];
            fRecoShowerMVAMuon[i_shower] = mvaOutMap["muon"];
            fRecoShowerMVAProton[i_shower] = mvaOutMap["proton"];
            fRecoShowerMVAPhoton[i_shower] = mvaOutMap["photon"];
        }
    }

    if (current_shower->dEdx().size() > 0)
    {
        for (int i_plane = 0; i_plane < 3; i_plane++){
          //fRecoShowerRecoEnergy[i_shower][i_plane] = fShowerEnergyAlg.ShowerEnergy(current_shower_hits, i_plane);
          fRecoShowerRecodEdx[i_shower][i_plane] = current_shower->dEdx()[i_plane];
          fRecoShowerRecoEnergy[i_shower][i_plane] = current_shower->Energy()[i_plane];
        }
    }

    // fill pandrizzle variables
    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(current_shower, evt));//fPandrizzleAlg.RunPID(current_shower, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), evt));

    fRecoShowerPandrizzleEvalRatio[i_shower]     = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
    fRecoShowerPandrizzleConcentration[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
    fRecoShowerPandrizzleCoreHaloRatio[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
    fRecoShowerPandrizzleConicalness[i_shower]   = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
    fRecoShowerPandrizzledEdxBestPlane[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
    fRecoShowerPandrizzleDisplacement[i_shower]  = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
    fRecoShowerPandrizzleDCA[i_shower]           = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
    fRecoShowerPandrizzleWideness[i_shower]      = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
    fRecoShowerPandrizzleEnergyDensity[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);

    fRecoShowerPandrizzlePathwayLengthMin[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kPathwayLengthMin);
    fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxShowerStartPathwayScatteringAngle2D);
    fRecoShowerPandrizzleMaxNPostShowerStartHits[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxNPostShowerStartHits);
    fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartScatterAngle);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartShowerStartEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[i_shower] = 
      pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
    fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinPostShowerStartShowerStartMoliereRadius);
    fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartOpeningAngle);
    fRecoShowerPandrizzleMaxFoundHitRatio[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxFoundHitRatio);
    fRecoShowerPandrizzleMaxInitialGapSize[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxInitialGapSize);
    fRecoShowerPandrizzleMinLargestProjectedGapSize[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinLargestProjectedGapSize);
    fRecoShowerPandrizzleNViewsWithAmbiguousHits[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kNViewsWithAmbiguousHits);
    fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kAmbiguousHitMaxUnaccountedEnergy);

    fRecoShowerPandrizzleBDTMethod[i_shower] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);
    float pandrizzleScore(pandrizzleRecord.GetMVAScore());
    fRecoShowerPandrizzleMVAScore[i_shower] = (std::fabs(fRecoShowerPandrizzleBDTMethod[i_shower] - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerJamPandrizzleMVAScore[i_shower] = (std::fabs(fRecoShowerPandrizzleBDTMethod[i_shower] - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerPandrizzleIsFilled[i_shower] = pandrizzleRecord.IsFilled();

    std::vector<art::Ptr<recob::Hit>> hitListU = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(current_shower_hits, 0);
    std::vector<art::Ptr<recob::Hit>> hitListV = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(current_shower_hits, 1);
    std::vector<art::Ptr<recob::Hit>> hitListW = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(current_shower_hits, 2);

    // Fill vertex separation info
    fRecoShowerDistanceToNuVertexU[i_shower] = DistanceToNuVertex(evt, hitListU, TVector3(fNuX, fNuY, fNuZ));
    fRecoShowerDistanceToNuVertexV[i_shower] = DistanceToNuVertex(evt, hitListV, TVector3(fNuX, fNuY, fNuZ));
    fRecoShowerDistanceToNuVertexW[i_shower] = DistanceToNuVertex(evt, hitListW, TVector3(fNuX, fNuY, fNuZ));
  }

    

  return;
}


void FDSelection::CCNuSelection::RunShowerSelection(art::Event const & evt)
{
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector <art::Ptr<recob::Shower>> showerList;
  if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle)))
  {
    std::cout << "Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
    return;
  } 
  else art::fill_ptr_vector(showerList, showerListHandle);

  // Get the selected shower
  art::Ptr<recob::Shower> sel_shower = fRecoShowerSelector->FindSelectedShower(evt);

  // If we didn't find a selected shower then what's the point?
  if (!(sel_shower.isAvailable())) 
  {
    std::cout << "FDSelection::CCNuSelection::RunShowerSelection - no shower selected by tool" << std::endl;
    return;
  }

  //10/10/18 DBrailsford start assessing PFParticle stuff
  //Now get the PFParticles
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle > > pfparticleList;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)))
  {
    std::cout << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel << std::endl;
    return;
  }
  else art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  art::FindManyP<recob::PFParticle> fmpfps(showerListHandle, evt, fShowerModuleLabel);
  const std::vector<art::Ptr<recob::PFParticle> > sel_shower_pfps = fmpfps.at(sel_shower.key());
  if (sel_shower_pfps.size() != 1)
  {
    std::cout << "CCNuSelection::GetPFParticleMatchedToShower NUMBER OF PFP MATCHED TO A SHOWER DOES NOT EQUAL 1: " << sel_shower_pfps.size() << std::endl;
    return;
  }

  art::Ptr<recob::PFParticle> sel_pfp = sel_shower_pfps[0];

  //Get the hits for said shower
  art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);
  const std::vector<art::Ptr<recob::Hit> > sel_shower_hits = fmhs.at(sel_shower.key());
  fSelShowerRecoNHits = sel_shower_hits.size();

  //Get the hit list handle
  art::Handle< std::vector<recob::Hit> > hitListHandle; 
  std::vector<art::Ptr<recob::Hit> > hitList; 
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitList, hitListHandle);

  //Is this PFParticle a primary daughter of the neutrino
  fSelShowerRecoIsPrimaryPFPDaughter = IsPFParticlePrimaryDaughter(sel_pfp, pfparticleList);

  fSelShowerRecoDirX   = sel_shower->Direction().X();
  fSelShowerRecoDirY   = sel_shower->Direction().Y();
  fSelShowerRecoDirZ   = sel_shower->Direction().Z();
  fSelShowerRecoStartX = sel_shower->ShowerStart().X();
  fSelShowerRecoStartY = sel_shower->ShowerStart().Y();
  fSelShowerRecoStartZ = sel_shower->ShowerStart().Z();

  fSelShowerRecoBestPlane = sel_shower->best_plane();
  fSelShowerRecoLength = sel_shower->Length();
  fSelShowerRecoOpeningAngle = sel_shower->OpenAngle();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  fSelShowerRecoCharge  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, sel_shower_hits); 

  //24/07/18 DBrailsford Use the data product to get the neutrino energy
  //fSelShowerRecoContained = IsTrackContained(sel_track, sel_track_hits, evt);

  //24/07/18 DBrailsford Get the reco energy data product for neutrinos
  //art::Handle<dune::EnergyRecoOutput> energyRecoHandle;
  //if (!evt.getByLabel(fNueEnergyRecoModuleLabel, energyRecoHandle)) {
  //std::cout<<"FDSelection::CCNuSelection::RunShowerSelection - Not able to find energy reconstruction container with name " << fNueEnergyRecoModuleLabel << std::endl;
  //return;
  //}

  std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(sel_shower, evt)));
  fNueRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNueRecoEHad = energyRecoHandle->fHadLorentzVector.E();
  fNueRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

  FillChildPFPInformation(sel_shower, evt, fSelShowerRecoNChildPFP, fSelShowerRecoNChildTrackPFP, fSelShowerRecoNChildShowerPFP);

  int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sel_shower_hits, 1);
  fSelShowerRecoCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, sel_shower_hits, hitList, g4id);
  fSelShowerRecoHitPurity = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, sel_shower_hits, g4id);

  if (TruthMatchUtils::Valid(g4id))
  {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      const simb::MCParticle* matched_mcparticle = pi_serv->TrackIdToParticle_P(g4id);
      if (matched_mcparticle){
        //Fill variables
        fSelShowerTruePDG = matched_mcparticle->PdgCode();
        if (matched_mcparticle->Mother()==0) fSelShowerTruePrimary = 1;
        else fSelShowerTruePrimary = 0;
        fSelShowerTrueMomX = matched_mcparticle->Momentum().X();
        fSelShowerTrueMomY = matched_mcparticle->Momentum().Y();
        fSelShowerTrueMomZ = matched_mcparticle->Momentum().Z();
        fSelShowerTrueMomT = matched_mcparticle->Momentum().T();
        fSelShowerTrueStartX = matched_mcparticle->Position(0).X();
        fSelShowerTrueStartY = matched_mcparticle->Position(0).Y();
        fSelShowerTrueStartZ = matched_mcparticle->Position(0).Z();
        fSelShowerTrueStartT = matched_mcparticle->Position(0).T();
        fSelShowerTrueEndX = matched_mcparticle->EndPosition().X();
        fSelShowerTrueEndY = matched_mcparticle->EndPosition().Y();
        fSelShowerTrueEndZ = matched_mcparticle->EndPosition().Z();
        fSelShowerTrueEndT = matched_mcparticle->EndPosition().T();
      }
  }
  //Now get the pid stuff
  art::FindManyP<anab::MVAPIDResult> fmpidt(showerListHandle, evt, fPIDModuleLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(sel_shower.key());
  if (pids.at(0).isAvailable()){
      std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
      if (!(mvaOutMap.empty())){
          //Get the PIDs
          fSelShowerMVAElectron = mvaOutMap["electron"];
          fSelShowerMVAPion = mvaOutMap["pich"];
          fSelShowerMVAMuon = mvaOutMap["muon"];
          fSelShowerMVAProton = mvaOutMap["proton"];
          fSelShowerMVAPhoton = mvaOutMap["photon"];
      }
  }

  //25/07/19 DBrailsford
  //Calculate shower energy
  if (sel_shower->dEdx().size() > 0)
  {
      for (unsigned int i_plane = 0; i_plane < 3; i_plane++){
        //fSelShowerRecoEnergy[i_plane] = fShowerEnergyAlg.ShowerEnergy(sel_shower_hits, i_plane);
        fSelShowerRecodEdx[i_plane] = sel_shower->dEdx()[i_plane];
        fSelShowerRecoEnergy[i_plane] = sel_shower->Energy()[i_plane];
      }
  }

  //Pandrizzle
  FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(sel_shower, evt));//fPandrizzleAlg.RunPID(sel_shower, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), evt));
  fSelShowerPandrizzleEvalRatio     = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
  fSelShowerPandrizzleConcentration = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
  fSelShowerPandrizzleCoreHaloRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
  fSelShowerPandrizzleConicalness   = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
  fSelShowerPandrizzledEdxBestPlane = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
  fSelShowerPandrizzleDisplacement  = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
  fSelShowerPandrizzleDCA           = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
  fSelShowerPandrizzleWideness      = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
  fSelShowerPandrizzleEnergyDensity = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);
  fSelShowerPandrizzleEvalRatio     = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
  fSelShowerPandrizzleConcentration = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
  fSelShowerPandrizzleCoreHaloRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
  fSelShowerPandrizzleConicalness   = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
  fSelShowerPandrizzledEdxBestPlane = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
  fSelShowerPandrizzleDisplacement  = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
  fSelShowerPandrizzleDCA           = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
  fSelShowerPandrizzleWideness      = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
  fSelShowerPandrizzleEnergyDensity = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);

  fSelShowerPandrizzlePathwayLengthMin = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kPathwayLengthMin);
  fSelShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxShowerStartPathwayScatteringAngle2D);
  fSelShowerPandrizzleMaxNPostShowerStartHits = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxNPostShowerStartHits);
  fSelShowerPandrizzleMaxPostShowerStartScatterAngle = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartScatterAngle);
  fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyAsymmetry);
  fSelShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartShowerStartEnergyAsymmetry);
  fSelShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = 
      pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
  fSelShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinPostShowerStartShowerStartMoliereRadius);
  fSelShowerPandrizzleMaxPostShowerStartOpeningAngle = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartOpeningAngle);
  fSelShowerPandrizzleMaxFoundHitRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxFoundHitRatio);
  fSelShowerPandrizzleMaxInitialGapSize = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxInitialGapSize);
  fSelShowerPandrizzleMinLargestProjectedGapSize = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinLargestProjectedGapSize);
  fSelShowerPandrizzleNViewsWithAmbiguousHits = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kNViewsWithAmbiguousHits);
  fSelShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kAmbiguousHitMaxUnaccountedEnergy);

  fSelShowerPandrizzleBDTMethod = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);
  float pandrizzleScore(pandrizzleRecord.GetMVAScore());
  fSelShowerPandrizzleMVAScore = (std::fabs(fSelShowerPandrizzleBDTMethod - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
  fSelShowerJamPandrizzleMVAScore = (std::fabs(fSelShowerPandrizzleBDTMethod - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
  fSelShowerPandrizzleIsFilled = pandrizzleRecord.IsFilled();
  fSelShowerPandrizzleIsFilled      = pandrizzleRecord.IsFilled();

  std::vector<art::Ptr<recob::Hit>> hitListU = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(sel_shower_hits, 0);
  std::vector<art::Ptr<recob::Hit>> hitListV = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(sel_shower_hits, 1);
  std::vector<art::Ptr<recob::Hit>> hitListW = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(sel_shower_hits, 2);

  // Fill vertex separation info
  fSelShowerDistanceToNuVertexU = DistanceToNuVertex(evt, hitListU, TVector3(fNuX, fNuY, fNuZ));
  fSelShowerDistanceToNuVertexV = DistanceToNuVertex(evt, hitListV, TVector3(fNuX, fNuY, fNuZ));
  fSelShowerDistanceToNuVertexW = DistanceToNuVertex(evt, hitListW, TVector3(fNuX, fNuY, fNuZ));

}

//double FDSelection::CCNuSelection::CalculateShowerCharge(art::Ptr<recob::Shower> const shower, std::vector< art::Ptr< recob::Hit> > const shower_hits){
//  double charge = 0;
//  for (unsigned int i_hit = 0; i_hit < shower_hits.size(); i_hit++){
//    if (shower_hits[i_hit]->WireID().Plane != 2) continue; //Collection only
//    charge += shower_hits[i_hit]->Integral() * fCalorimetryAlg.LifetimeCorrection(shower_hits[i_hit]->PeakTime(), fT0);
//  }
//  return charge;
//}

bool FDSelection::CCNuSelection::IsPFParticlePrimaryDaughter(art::Ptr<recob::PFParticle> const & pfparticle, std::vector<art::Ptr<recob::PFParticle > > const & pfparticles){
  size_t parent_id = pfparticle->Parent();
  for (unsigned int i_pfp = 0; i_pfp < pfparticles.size(); i_pfp++){
    art::Ptr<recob::PFParticle> parent_candidate = pfparticles[i_pfp];
    if (parent_candidate->Self() == parent_id && parent_candidate->IsPrimary()){
      return true;
    }
  }
  return false;
}

TVector3 FDSelection::CCNuSelection::ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector){
  TVector3 projected_vector = vector_to_project - (vector_to_project.Dot(plane_norm_vector)/(plane_norm_vector.Mag()*plane_norm_vector.Mag()))*plane_norm_vector;
  return projected_vector;
}

//////////////////////////////////////////////////////////////////

double FDSelection::CCNuSelection::DistanceToNuVertex(art::Event const & evt, std::vector<art::Ptr<recob::Hit>> const artHitList, const TVector3 &nuVertex)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    double closestDistanceSquared(std::numeric_limits<double>::max());

    for (art::Ptr<recob::Hit> hit : artHitList)
    {
        const geo::WireID hit_WireID(hit->WireID());
        const double hit_Time(hit->PeakTime());
        const geo::View_t hit_View(hit->View());
        const geo::CryostatID cryostatID(hit_WireID.Cryostat);

        const geo::View_t pandora_View(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));

        TVector3 vertexXYZ = nuVertex;
        geo::Point_t hitXYZ = theGeometry->Cryostat(cryostatID).TPC(hit_WireID.TPC).Plane(hit_WireID.Plane).Wire(hit_WireID.Wire).GetCenter();
        const double hitY(hitXYZ.Y());
        const double hitZ(hitXYZ.Z());

        double hitPosition[3];
        hitPosition[0] = detProp.ConvertTicksToX(hit_Time, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat);
        hitPosition[1] = 0.0;

        double nuVertexPosition[3];
        nuVertexPosition[0] = nuVertex[0];
        nuVertexPosition[1] = 0.0;

        if (pandora_View == geo::kW)
        {
            hitPosition[2] = YZtoW(hit_WireID, hitY, hitZ);
            nuVertexPosition[2] = YZtoW(hit_WireID, vertexXYZ[1], vertexXYZ[2]);
        }
        else if (pandora_View == geo::kU)
        {
            hitPosition[2] = YZtoU(hit_WireID, hitY, hitZ);
            nuVertexPosition[2] = YZtoU(hit_WireID, vertexXYZ[1], vertexXYZ[2]);
        }
        else if (pandora_View == geo::kV) 
        {
            hitPosition[2] = YZtoV(hit_WireID, hitY, hitZ);
            nuVertexPosition[2] = YZtoV(hit_WireID, vertexXYZ[1], vertexXYZ[2]);
        }
        else 
        {
            throw cet::exception("LArPandora")
                << "CreatePandoraHits2D - this wire view not recognised (View=" << hit_View << ") ";
        }

        double deltaX = hitPosition[0] - nuVertexPosition[0];
        double deltaY = hitPosition[1] - nuVertexPosition[1];
        double deltaZ = hitPosition[2] - nuVertexPosition[2];
        double separationSquared = (deltaX * deltaX) + (deltaY * deltaY) + (deltaZ * deltaZ);

        if (separationSquared < closestDistanceSquared)
            closestDistanceSquared = separationSquared;
    }

    return std::sqrt(closestDistanceSquared);
}

////////////////////////////////////////////////////////////////////////////////////////////////

double FDSelection::CCNuSelection::YZtoU(const geo::WireID &hit_WireID, double y, double z)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    //const double wireAngleU(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kU, hit_WireID.TPC, hit_WireID.Cryostat));
    const double wireAngleU(0.623257);
    return (z * std::cos(wireAngleU) - y * std::sin(wireAngleU));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double FDSelection::CCNuSelection::YZtoV(const geo::WireID &hit_WireID, double y, double z)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    //const double wireAngleV(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kV, hit_WireID.TPC, hit_WireID.Cryostat));
    const double wireAngleV(-0.623257);
    return (z * std::cos(wireAngleV) - y * std::sin(wireAngleV));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double FDSelection::CCNuSelection::YZtoW(const geo::WireID &hit_WireID, double y, double z)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    //const double wireAngleW(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kW, hit_WireID.TPC, hit_WireID.Cryostat));
    const double wireAngleW(0.0);
    return (z * std::cos(wireAngleW) - y * std::sin(wireAngleW));
}

////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(FDSelection::CCNuSelection)
