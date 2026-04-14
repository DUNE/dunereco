///////////////////////////////////////////////
// CCNuSelection analyzer
//
// Creates tree on which to perform Pandora-based numu/nue selection
// I Mawby, D Brailsford May 2023
///////////////////////////////////////////////

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
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
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
#include "larsim/Utils/TruthMatchUtils.h"

//DUNE
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/TrackPID/algorithms/CTPHelper.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
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

constexpr int kDefMaxNTrueVertexParticles = 150;
constexpr int kDefMaxNRecoTracks = 50;
constexpr int kDefMaxNRecoShowers = 50;

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
  void Reset();
  void GetEventInfo(art::Event const & evt);
  void GetTruthInfo(art::Event const & evt);
  void FillVertexInfo(art::Event const & evt);
  void GetRecoTrackInfo(art::Event const & evt);
  void RunTrackSelection(art::Event const & evt);
  void FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp);
  void GetRecoShowerInfo(art::Event const & evt);
  void RunShowerSelection(art::Event const & evt);

  TVector3 ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector);
  double DistanceToNuVertex(art::Event const & evt, std::vector<art::Ptr<recob::Hit>> const artHitList, const TVector3 &nuVertex);
  double YZtoU(const geo::WireID &hit_WireID, double y, double x);
  double YZtoV(const geo::WireID &hit_WireID, double y, double x);
  double YZtoW(const geo::WireID &hit_WireID, double y, double x);

  // Trees
  TTree *fTree; // the selection tree
  TTree* fPOTTree;
  double fPOT;
  ////////////////////////////////////////
  // Event info
  // Event identification
  int fRun;
  int fSubRun;
  int fEvent;
  int fIsMC;
  //Detector 
  double fT0;
  //Neutrino 
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
  double fNuMomX; //Neutrino momentum
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
  int fVertexParticlePDG[kDefMaxNTrueVertexParticles];
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
  //Outgoing Lepton 
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
  //Transverse Mom
  double  fNuMomTranMag;
  double  fTargNuclMomTranMag;
  double  fInitalMomTranMag;
  double  fLepMomTranMag;
  double  fNuclRemMomTranMag;
  double  fFinalMomTranMagNoLepNoRem;
  double  fFinalMomTranMagNoLepWithRem;
  double  fFinalMomTranMagWithLepNoRem;
  double  fFinalMomTranMagWithLepWithRem;
  ////////////////////////////////////////
  // Selected Track Info
  // Truth
  int fSelTrackTruePDG;
  bool fSelTrackTruePrimary;
  double fSelTrackTrueMomX;
  double fSelTrackTrueMomY;
  double fSelTrackTrueMomZ;
  double fSelTrackTrueMomT;
  double fSelTrackTrueStartX;
  double fSelTrackTrueStartY;
  double fSelTrackTrueStartZ;
  double fSelTrackTrueEndX;
  double fSelTrackTrueEndY;
  double fSelTrackTrueEndZ;
  // Reco
  int fSelTrackRecoNHits;
  double fSelTrackRecoCompleteness;
  double fSelTrackRecoHitPurity;
  double fSelTrackRecoStartX;
  double fSelTrackRecoStartY;
  double fSelTrackRecoStartZ;
  double fSelTrackRecoEndX;
  double fSelTrackRecoEndY;
  double fSelTrackRecoEndZ;
  double fSelTrackRecoUpstreamX;
  double fSelTrackRecoUpstreamY;
  double fSelTrackRecoUpstreamZ;
  double fSelTrackRecoDownstreamX;
  double fSelTrackRecoDownstreamY;
  double fSelTrackRecoDownstreamZ;
  double fSelTrackRecoEndClosestToVertexX;
  double fSelTrackRecoEndClosestToVertexY;
  double fSelTrackRecoEndClosestToVertexZ;
  double fSelTrackRecoLength;
  int   fSelTrackRecoContained;
  int   fSelTrackRecoMomMethod;
  double fSelTrackRecoCharge;
  double fSelTrackRecoMomMCS;
  double fSelTrackRecoMomContained;
  double fSelTrackRecoVertexX;
  double fSelTrackRecoVertexY;
  double fSelTrackRecoVertexZ;
  int fSelTrackRecoNChildPFP;
  int fSelTrackRecoNChildTrackPFP;
  int fSelTrackRecoNChildShowerPFP;
  // MVA PID
  double fSelTrackMVAElectron;
  double fSelTrackMVAPion;
  double fSelTrackMVAMuon;
  double fSelTrackMVAProton;
  double fSelTrackMVAPhoton;
  // Pandizzle PID
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
  // DeepPan PID
  double fSelTrackDeepPanMuVar;
  double fSelTrackDeepPanPiVar;
  double fSelTrackDeepPanProtonVar;
  ////////////////////////////////////////
  // All tracks info
  // Truth
  int fRecoTrackTruePDG[kDefMaxNRecoTracks];
  bool fRecoTrackTruePrimary[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomX[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomY[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomZ[kDefMaxNRecoTracks];
  double fRecoTrackTrueMomT[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartX[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartY[kDefMaxNRecoTracks];
  double fRecoTrackTrueStartZ[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndX[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndY[kDefMaxNRecoTracks];
  double fRecoTrackTrueEndZ[kDefMaxNRecoTracks];
  // Reco
  int fNRecoTracks;
  bool fRecoTrackIsPrimary[kDefMaxNRecoTracks];
  int fRecoTrackRecoNHits[kDefMaxNRecoTracks];
  double fRecoTrackRecoCompleteness[kDefMaxNRecoTracks];
  double fRecoTrackRecoHitPurity[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartX[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartY[kDefMaxNRecoTracks];
  double fRecoTrackRecoStartZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndX[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndY[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamX[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamY[kDefMaxNRecoTracks];
  double fRecoTrackRecoUpstreamZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamX[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamY[kDefMaxNRecoTracks];
  double fRecoTrackRecoDownstreamZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndClosestToVertexX[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndClosestToVertexY[kDefMaxNRecoTracks];
  double fRecoTrackRecoEndClosestToVertexZ[kDefMaxNRecoTracks];
  double fRecoTrackRecoLength[kDefMaxNRecoTracks];
  int   fRecoTrackRecoContained[kDefMaxNRecoTracks];
  int   fRecoTrackRecoMomMethod[kDefMaxNRecoTracks];
  double fRecoTrackRecoCharge[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomMCS[kDefMaxNRecoTracks];
  double fRecoTrackRecoMomContained[kDefMaxNRecoTracks];
  double fRecoTrackRecoVertexX[kDefMaxNRecoTracks];
  double fRecoTrackRecoVertexY[kDefMaxNRecoTracks];
  double fRecoTrackRecoVertexZ[kDefMaxNRecoTracks];
  int fRecoTrackRecoNChildPFP[kDefMaxNRecoTracks];
  int fRecoTrackRecoNChildTrackPFP[kDefMaxNRecoTracks];
  int fRecoTrackRecoNChildShowerPFP[kDefMaxNRecoTracks];
  // MVA PID
  double fRecoTrackMVAElectron[kDefMaxNRecoTracks];
  double fRecoTrackMVAPion[kDefMaxNRecoTracks];
  double fRecoTrackMVAMuon[kDefMaxNRecoTracks];
  double fRecoTrackMVAProton[kDefMaxNRecoTracks];
  double fRecoTrackMVAPhoton[kDefMaxNRecoTracks];
  // DeepPan PID
  double fRecoTrackDeepPanMuVar[kDefMaxNRecoTracks];
  double fRecoTrackDeepPanPiVar[kDefMaxNRecoTracks];
  double fRecoTrackDeepPanProtonVar[kDefMaxNRecoTracks];
  // Pandizzle PID
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
  ////////////////////////////////////////
  // Selected shower info
  // Truth
  int fSelShowerTruePDG;
  bool fSelShowerTruePrimary;
  double fSelShowerTrueMomX;
  double fSelShowerTrueMomY;
  double fSelShowerTrueMomZ;
  double fSelShowerTrueMomT;
  double fSelShowerTrueStartX;
  double fSelShowerTrueStartY;
  double fSelShowerTrueStartZ;
  double fSelShowerTrueEndX;
  double fSelShowerTrueEndY;
  double fSelShowerTrueEndZ;
  // Reco
  int fSelShowerRecoNHits;
  double fSelShowerRecoCompleteness;
  double fSelShowerRecoHitPurity;
  double fSelShowerRecoStartX;
  double fSelShowerRecoStartY;
  double fSelShowerRecoStartZ;
  // trj -- comment out unused variables to make clang succeed
  //double fSelShowerRecoEndX;
  //double fSelShowerRecoEndY;
  //double fSelShowerRecoEndZ;
  double fSelShowerRecoVertexX;
  double fSelShowerRecoVertexY;
  double fSelShowerRecoVertexZ;
  double fSelShowerDistanceToNuVertexU;
  double fSelShowerDistanceToNuVertexV;
  double fSelShowerDistanceToNuVertexW;
  double fSelShowerRecoCharge;
  double fSelShowerRecoMom;
  int fSelShowerRecoNChildPFP;
  int fSelShowerRecoNChildTrackPFP;
  int fSelShowerRecoNChildShowerPFP;
  double fSelShowerRecoDirX;
  double fSelShowerRecoDirY;
  double fSelShowerRecoDirZ;
  double fSelShowerRecodEdx[3];
  double fSelShowerRecoEnergy[3];
  int fSelShowerRecoBestPlane;
  double fSelShowerRecoLength;
  double fSelShowerRecoOpeningAngle;
  // MVA PID
  double fSelShowerMVAElectron;
  double fSelShowerMVAPion;
  double fSelShowerMVAMuon;
  double fSelShowerMVAProton;
  double fSelShowerMVAPhoton;
  // Pandrizzle PID
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
  double fSelShowerEnhancedPandrizzleScore;
  double fSelShowerBackupPandrizzleScore;
  bool   fSelShowerPandrizzleIsFilled;
  ////////////////////////////////////////
  // All showers info
  // True
  int fNRecoShowers;
  bool fRecoShowerTruePrimary[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomX[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomY[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomZ[kDefMaxNRecoShowers];
  double fRecoShowerTrueMomT[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartX[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartY[kDefMaxNRecoShowers];
  double fRecoShowerTrueStartZ[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndX[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndY[kDefMaxNRecoShowers];
  double fRecoShowerTrueEndZ[kDefMaxNRecoShowers];
  // Reco
  int fRecoShowerTruePDG[kDefMaxNRecoShowers];
  int fRecoShowerRecoNHits[kDefMaxNRecoShowers];
  double fRecoShowerRecoCompleteness[kDefMaxNRecoShowers];
  double fRecoShowerRecoHitPurity[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartX[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartY[kDefMaxNRecoShowers];
  double fRecoShowerRecoStartZ[kDefMaxNRecoShowers];
  double fRecoShowerRecoCharge[kDefMaxNRecoShowers];
  double fRecoShowerRecoMom[kDefMaxNRecoShowers];
  double fRecoShowerRecoVertexX[kDefMaxNRecoShowers];
  double fRecoShowerRecoVertexY[kDefMaxNRecoShowers];
  double fRecoShowerRecoVertexZ[kDefMaxNRecoShowers];
  double fRecoShowerDistanceToNuVertexU[kDefMaxNRecoShowers];
  double fRecoShowerDistanceToNuVertexV[kDefMaxNRecoShowers];
  double fRecoShowerDistanceToNuVertexW[kDefMaxNRecoShowers];
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
  double fRecoShowerRecoEnergy[kDefMaxNRecoShowers][3];
  // MVA PID
  double fRecoShowerMVAElectron[kDefMaxNRecoShowers];
  double fRecoShowerMVAPion[kDefMaxNRecoShowers];
  double fRecoShowerMVAMuon[kDefMaxNRecoShowers];
  double fRecoShowerMVAProton[kDefMaxNRecoShowers];
  double fRecoShowerMVAPhoton[kDefMaxNRecoShowers];
  // Pandrizzle MVA
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
  double fRecoShowerEnhancedPandrizzleScore[kDefMaxNRecoShowers];
  double fRecoShowerBackupPandrizzleScore[kDefMaxNRecoShowers];
  bool   fRecoShowerPandrizzleIsFilled[kDefMaxNRecoShowers];
  ////////////////////////////////////////
  //Event-level
  double fRecoNuVtxX;
  double fRecoNuVtxY;
  double fRecoNuVtxZ;
  int fRecoNuVtxNShowers;
  int fRecoNuVtxNTracks;
  int fRecoNuVtxNChildren;
  double fRecoEventCharge;
  double fNumuRecoMomLep;
  double fNumuRecoEHad;
  double fNumuRecoENu;
  double fNueRecoMomLep;
  double fNueRecoEHad;
  double fNueRecoENu;
  ////////////////////////////////////////
  //Module labels
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
  std::string fCVNProductInstance;
  std::string fNumuEnergyRecoModuleLabel;
  std::string fNueEnergyRecoModuleLabel;
  ////////////////////////////////////////
  //Algs
  bool fMakeSelectionTrainingTrees;
  double fReducedTreeMode;
  PandizzleAlg fPandizzleAlg;
  PandrizzleAlg fPandrizzleAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  ctp::CTPHelper fConvTrackPID;
  dune::NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;

  //Tools
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;
  std::unique_ptr<FDSelectionTools::RecoShowerSelector> fRecoShowerSelector;

  // DUNE CVN Scores
  double fCVNResultNue;
  double fCVNResultNumu;
  double fCVNResultNutau;
  double fCVNResultNC;
};

FDSelection::CCNuSelection::CCNuSelection(fhicl::ParameterSet const & pset) :
  EDAnalyzer(pset),
  fNVertexParticles(0),
  fNRecoTracks(0),
  fNRecoShowers(0),
  fNuGenModuleLabel(pset.get<std::string>("ModuleLabels.NuGenModuleLabel")),
  fLargeantModuleLabel(pset.get<std::string>("ModuleLabels.LargeantModuleLabel")),
  fWireModuleLabel(pset.get<std::string>("ModuleLabels.WireModuleLabel")),
  fTrackModuleLabel(pset.get<std::string>("ModuleLabels.TrackModuleLabel")),
  fShowerModuleLabel(pset.get<std::string>("ModuleLabels.ShowerModuleLabel")),
  fPFParticleModuleLabel(pset.get<std::string>("ModuleLabels.PFParticleModuleLabel")),
  fVertexModuleLabel(pset.get<std::string>("ModuleLabels.VertexModuleLabel")),
  fPIDModuleLabel(pset.get<std::string>("ModuleLabels.PIDModuleLabel")),
  fParticleIDModuleLabel(pset.get<std::string>("ModuleLabels.ParticleIDModuleLabel")),
  fHitsModuleLabel(pset.get<std::string>("ModuleLabels.HitsModuleLabel")),
  fPOTModuleLabel(pset.get<std::string>("ModuleLabels.POTModuleLabel")),
  fCVNModuleLabel(pset.get<std::string>("ModuleLabels.CVNModuleLabel")),
  fCVNProductInstance(pset.get<std::string>("ModuleLabels.CVNProductInstance")),
  fNumuEnergyRecoModuleLabel(pset.get<std::string>("ModuleLabels.NumuEnergyRecoModuleLabel")),
  fNueEnergyRecoModuleLabel(pset.get<std::string>("ModuleLabels.NueEnergyRecoModuleLabel")),
  fMakeSelectionTrainingTrees(pset.get<bool>("MakeSelectionTrainingTrees")),
  fReducedTreeMode(pset.get<bool>("ReducedTreeMode", true)),
  fPandizzleAlg(pset),
  fPandrizzleAlg(pset),
  fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fConvTrackPID(pset.get<fhicl::ParameterSet>("ctpHelper")),
  fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"),fTrackModuleLabel,fShowerModuleLabel,
  fHitsModuleLabel,fWireModuleLabel,fTrackModuleLabel,fShowerModuleLabel,fPFParticleModuleLabel),
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(pset.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))}
{
}

void FDSelection::CCNuSelection::analyze(art::Event const & evt)
{
  Reset(); // reset tree variables

  fRun = evt.run();
  fSubRun = evt.subRun();
  fEvent = evt.event();
  fIsMC = !evt.isRealData();

  GetEventInfo(evt);

  if (fIsMC) 
    GetTruthInfo(evt);

  FillVertexInfo(evt);
  GetRecoTrackInfo(evt);
  RunTrackSelection(evt);
  GetRecoShowerInfo(evt);
  RunShowerSelection(evt);

  if (fMakeSelectionTrainingTrees)
  {
    fPandizzleAlg.Run(evt);
    fPandrizzleAlg.Run(evt);
  }

  fTree->Fill();
}

void FDSelection::CCNuSelection::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("ccnusel","CC nu selection");

    // Event info
    fTree->Branch("Run",&fRun);
    fTree->Branch("SubRun",&fSubRun);
    fTree->Branch("Event",&fEvent);

    // True info
    fTree->Branch("NuPdg",&fNuPdg);
    fTree->Branch("BeamPdg",&fBeamPdg);
    fTree->Branch("NC",&fNC);
    fTree->Branch("Mode",&fMode);
    fTree->Branch("TargetZ",&fTargetZ);
    fTree->Branch("Enu",&fENu);
    fTree->Branch("NuX",&fNuX);
    fTree->Branch("NuY",&fNuY);
    fTree->Branch("NuZ",&fNuZ);

    // Reco info 
    fTree->Branch("RecoNuVtxX",&fRecoNuVtxX);
    fTree->Branch("RecoNuVtxY",&fRecoNuVtxY);
    fTree->Branch("RecoNuVtxZ",&fRecoNuVtxZ);
    fTree->Branch("NueRecoENu",&fNueRecoENu);
    fTree->Branch("NumuRecoENu",&fNumuRecoENu);

    // CVN info
    fTree->Branch("CVNResultNue", &fCVNResultNue);
    fTree->Branch("CVNResultNumu", &fCVNResultNumu);
    fTree->Branch("CVNResultNutau", &fCVNResultNutau);
    fTree->Branch("CVNResultNC", &fCVNResultNC);

    // Selected track things...
    fTree->Branch("SelTrackTruePDG",&fSelTrackTruePDG);
    fTree->Branch("SelTrackRecoContained",&fSelTrackRecoContained);
    fTree->Branch("SelTrackRecoMomMethod",&fSelTrackRecoMomMethod);
    fTree->Branch("SelTrackPandizzleVar", &fSelTrackPandizzleVar);

    // Selected shower things
    fTree->Branch("SelShowerTruePDG",&fSelShowerTruePDG);
    fTree->Branch("SelShowerEnhancedPandrizzleScore",&fSelShowerEnhancedPandrizzleScore);
    fTree->Branch("SelShowerBackupPandrizzleScore",&fSelShowerBackupPandrizzleScore);

    // POT tree
    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("POT",&fPOT);
    fPOTTree->Branch("Run",&fRun);
    fPOTTree->Branch("SubRun",&fSubRun);

    if (!fReducedTreeMode)
    {
        fTree->Branch("IsMC",&fIsMC);
        fTree->Branch("T0",&fT0);
        fTree->Branch("Q2",&fQ2);
        fTree->Branch("W",&fW);
        fTree->Branch("X",&fX);
        fTree->Branch("Y",&fY);
        fTree->Branch("NuMomX",&fNuMomX);
        fTree->Branch("NuMomY",&fNuMomY);
        fTree->Branch("NuMomZ",&fNuMomZ);
        fTree->Branch("NuMomT",&fNuMomT);
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

        fTree->Branch("SelTrackTruePrimary",&fSelTrackTruePrimary);
        fTree->Branch("SelTrackTrueMomX",&fSelTrackTrueMomX);
        fTree->Branch("SelTrackTrueMomY",&fSelTrackTrueMomY);
        fTree->Branch("SelTrackTrueMomZ",&fSelTrackTrueMomZ);
        fTree->Branch("SelTrackTrueMomT",&fSelTrackTrueMomT);
        fTree->Branch("SelTrackTrueStartX",&fSelTrackTrueStartX);
        fTree->Branch("SelTrackTrueStartY",&fSelTrackTrueStartY);
        fTree->Branch("SelTrackTrueStartZ",&fSelTrackTrueStartZ);
        fTree->Branch("SelTrackTrueEndX",&fSelTrackTrueEndX);
        fTree->Branch("SelTrackTrueEndY",&fSelTrackTrueEndY);
        fTree->Branch("SelTrackTrueEndZ",&fSelTrackTrueEndZ);
        fTree->Branch("SelTrackRecoNHits",&fSelTrackRecoNHits);
        fTree->Branch("SelTrackRecoCompleteness",&fSelTrackRecoCompleteness);
        fTree->Branch("SelTrackRecoHitPurity",&fSelTrackRecoHitPurity);
        fTree->Branch("SelTrackRecoStartX",&fSelTrackRecoStartX);
        fTree->Branch("SelTrackRecoStartY",&fSelTrackRecoStartY);
        fTree->Branch("SelTrackRecoStartZ",&fSelTrackRecoStartZ);
        fTree->Branch("SelTrackRecoEndX",&fSelTrackRecoEndX);
        fTree->Branch("SelTrackRecoEndY",&fSelTrackRecoEndY);
        fTree->Branch("SelTrackRecoEndZ",&fSelTrackRecoEndZ);
        fTree->Branch("SelTrackRecoUpstreamX",&fSelTrackRecoUpstreamX);
        fTree->Branch("SelTrackRecoUpstreamY",&fSelTrackRecoUpstreamY);
        fTree->Branch("SelTrackRecoUpstreamZ",&fSelTrackRecoUpstreamZ);
        fTree->Branch("SelTrackRecoDownstreamX",&fSelTrackRecoDownstreamX);
        fTree->Branch("SelTrackRecoDownstreamY",&fSelTrackRecoDownstreamY);
        fTree->Branch("SelTrackRecoDownstreamZ",&fSelTrackRecoDownstreamZ);
        fTree->Branch("SelTrackRecoEndClosestToVertexX",&fSelTrackRecoEndClosestToVertexX);
        fTree->Branch("SelTrackRecoEndClosestToVertexY",&fSelTrackRecoEndClosestToVertexY);
        fTree->Branch("SelTrackRecoEndClosestToVertexZ",&fSelTrackRecoEndClosestToVertexZ);
        fTree->Branch("SelTrackRecoLength",&fSelTrackRecoLength);
        fTree->Branch("SelTrackRecoCharge",&fSelTrackRecoCharge);
        fTree->Branch("SelTrackRecoMomMCS",&fSelTrackRecoMomMCS);
        fTree->Branch("SelTrackRecoMomContained",&fSelTrackRecoMomContained);
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
        fTree->Branch("RecoTrackTrueEndX",fRecoTrackTrueEndX,"RecoTrackTrueEndX[NRecoTracks]/D");
        fTree->Branch("RecoTrackTrueEndY",fRecoTrackTrueEndY,"RecoTrackTrueEndY[NRecoTracks]/D");
        fTree->Branch("RecoTrackTrueEndZ",fRecoTrackTrueEndZ,"RecoTrackTrueEndZ[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoNHits",fRecoTrackRecoNHits,"RecoTrackRecoNHits[NRecoTracks]/I");
        fTree->Branch("RecoTrackRecoCompleteness",fRecoTrackRecoCompleteness,"RecoTrackRecoCompleteness[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoHitPurity",fRecoTrackRecoHitPurity,"RecoTrackRecoHitPurity[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoStartX",fRecoTrackRecoStartX,"RecoTrackRecoStartX[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoStartY",fRecoTrackRecoStartY,"RecoTrackRecoStartY[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoStartZ",fRecoTrackRecoStartZ,"RecoTrackRecoStartZ[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoEndX",fRecoTrackRecoEndX,"RecoTrackRecoEndX[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoEndY",fRecoTrackRecoEndY,"RecoTrackRecoEndY[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoEndZ",fRecoTrackRecoEndZ,"RecoTrackRecoEndZ[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoUpstreamX",fRecoTrackRecoUpstreamX,"RecoTrackRecoUpstreamX[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoUpstreamY",fRecoTrackRecoUpstreamY,"RecoTrackRecoUpstreamY[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoUpstreamZ",fRecoTrackRecoUpstreamZ,"RecoTrackRecoUpstreamZ[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoDownstreamX",fRecoTrackRecoDownstreamX,"RecoTrackRecoDownstreamX[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoDownstreamY",fRecoTrackRecoDownstreamY,"RecoTrackRecoDownstreamY[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoDownstreamZ",fRecoTrackRecoDownstreamZ,"RecoTrackRecoDownstreamZ[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoEndClosestToVertexX",fRecoTrackRecoEndClosestToVertexX,"RecoTrackRecoEndClosestToVertexX[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoEndClosestToVertexY",fRecoTrackRecoEndClosestToVertexY,"RecoTrackRecoEndClosestToVertexY[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoEndClosestToVertexZ",fRecoTrackRecoEndClosestToVertexZ,"RecoTrackRecoEndClosestToVertexZ[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoLength",fRecoTrackRecoLength,"RecoTrackRecoLength[NRecoTracks]/D");
        fTree->Branch("RecoTrackRecoContained",fRecoTrackRecoContained,"RecoTrackRecoContained[NRecoTracks]/I");
        fTree->Branch("RecoTrackRecoMomMethod",fRecoTrackRecoMomMethod,"RecoTrackRecoMomMethod[NRecoTracks]/I");
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
        fTree->Branch("SelShowerTruePrimary",&fSelShowerTruePrimary);
        fTree->Branch("SelShowerTrueMomX",&fSelShowerTrueMomX);
        fTree->Branch("SelShowerTrueMomY",&fSelShowerTrueMomY);
        fTree->Branch("SelShowerTrueMomZ",&fSelShowerTrueMomZ);
        fTree->Branch("SelShowerTrueMomT",&fSelShowerTrueMomT);
        fTree->Branch("SelShowerTrueStartX",&fSelShowerTrueStartX);
        fTree->Branch("SelShowerTrueStartY",&fSelShowerTrueStartY);
        fTree->Branch("SelShowerTrueStartZ",&fSelShowerTrueStartZ);
        fTree->Branch("SelShowerTrueEndX",&fSelShowerTrueEndX);
        fTree->Branch("SelShowerTrueEndY",&fSelShowerTrueEndY);
        fTree->Branch("SelShowerTrueEndZ",&fSelShowerTrueEndZ);
        fTree->Branch("SelShowerRecoNHits",&fSelShowerRecoNHits);
        fTree->Branch("SelShowerRecoCompleteness",&fSelShowerRecoCompleteness);
        fTree->Branch("SelShowerRecoHitPurity",&fSelShowerRecoHitPurity);
        fTree->Branch("SelShowerRecoStartX",&fSelShowerRecoStartX);
        fTree->Branch("SelShowerRecoStartY",&fSelShowerRecoStartY);
        fTree->Branch("SelShowerRecoStartZ",&fSelShowerRecoStartZ);
        fTree->Branch("SelShowerRecoCharge",&fSelShowerRecoCharge);
        fTree->Branch("SelShowerRecoMom",&fSelShowerRecoMom);
        fTree->Branch("SelShowerRecoVertexX",&fSelShowerRecoVertexX);
        fTree->Branch("SelShowerRecoVertexY",&fSelShowerRecoVertexY);
        fTree->Branch("SelShowerRecoVertexZ",&fSelShowerRecoVertexZ);
        fTree->Branch("SelShowerDistanceToNuVertexU", &fSelShowerDistanceToNuVertexU);
        fTree->Branch("SelShowerDistanceToNuVertexV", &fSelShowerDistanceToNuVertexV);
        fTree->Branch("SelShowerDistanceToNuVertexW", &fSelShowerDistanceToNuVertexW);
        fTree->Branch("SelShowerRecoNChildPFP",&fSelShowerRecoNChildPFP);
        fTree->Branch("SelShowerRecoNChildTrackPFP",&fSelShowerRecoNChildTrackPFP);
        fTree->Branch("SelShowerRecoNChildShowerPFP",&fSelShowerRecoNChildShowerPFP);
        fTree->Branch("SelShowerRecoDirX",&fSelShowerRecoDirX);
        fTree->Branch("SelShowerRecoDirY",&fSelShowerRecoDirY);
        fTree->Branch("SelShowerRecoDirZ",&fSelShowerRecoDirZ);
        fTree->Branch("SelShowerRecodEdx",&fSelShowerRecodEdx,"SelShowerRecodEdx[3]/D");
        fTree->Branch("SelShowerRecoEnergy",&fSelShowerRecoEnergy,"SelShowerRecoEnergy[3]/D");
        fTree->Branch("SelShowerRecoBestPlane",&fSelShowerRecoBestPlane);
        fTree->Branch("SelShowerRecoLength",&fSelShowerRecoLength);
        fTree->Branch("SelShowerRecoOpeningAngle",&fSelShowerRecoOpeningAngle);
        fTree->Branch("SelShowerMVAElectron",&fSelShowerMVAElectron);
        fTree->Branch("SelShowerMVAPion",&fSelShowerMVAPion);
        fTree->Branch("SelShowerMVAMuon",&fSelShowerMVAMuon);
        fTree->Branch("SelShowerMVAProton",&fSelShowerMVAProton);
        fTree->Branch("SelShowerMVAPhoton",&fSelShowerMVAPhoton);
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
        fTree->Branch("SelShowerPandrizzleIsFilled",&fSelShowerPandrizzleIsFilled);

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
        fTree->Branch("RecoShowerTrueEndX",fRecoShowerTrueEndX,"RecoShowerTrueEndX[NRecoShowers]/D");
        fTree->Branch("RecoShowerTrueEndY",fRecoShowerTrueEndY,"RecoShowerTrueEndY[NRecoShowers]/D");
        fTree->Branch("RecoShowerTrueEndZ",fRecoShowerTrueEndZ,"RecoShowerTrueEndZ[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoNHits",fRecoShowerRecoNHits,"RecoShowerRecoNHits[NRecoShowers]/I");
        fTree->Branch("RecoShowerRecoCompleteness",fRecoShowerRecoCompleteness,"RecoShowerRecoCompleteness[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoHitPurity",fRecoShowerRecoHitPurity,"RecoShowerRecoHitPurity[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoStartX",fRecoShowerRecoStartX,"RecoShowerRecoStartX[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoStartY",fRecoShowerRecoStartY,"RecoShowerRecoStartY[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoStartZ",fRecoShowerRecoStartZ,"RecoShowerRecoStartZ[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoCharge",fRecoShowerRecoCharge,"RecoShowerRecoCharge[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoMom",fRecoShowerRecoMom,"RecoShowerRecoMom[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoVertexX",fRecoShowerRecoVertexX,"RecoShowerRecoVertexX[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoVertexY",fRecoShowerRecoVertexY,"RecoShowerRecoVertexY[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoVertexZ",fRecoShowerRecoVertexZ,"RecoShowerRecoVertexZ[NRecoShowers]/D");
        fTree->Branch("RecoShowerDistanceToNuVertexU", &fRecoShowerDistanceToNuVertexU, "RecoShowerDistanceToNuVertexU[NRecoShowers]/D");
        fTree->Branch("RecoShowerDistanceToNuVertexV", &fRecoShowerDistanceToNuVertexV, "RecoShowerDistanceToNuVertexV[NRecoShowers]/D");
        fTree->Branch("RecoShowerDistanceToNuVertexW", &fRecoShowerDistanceToNuVertexW, "RecoShowerDistanceToNuVertexW[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoNChildPFP",fRecoShowerRecoNChildPFP,"RecoShowerRecoNChildPFP[NRecoShowers]/I");
        fTree->Branch("RecoShowerRecoNChildTrackPFP",fRecoShowerRecoNChildTrackPFP,"RecoShowerRecoNChildTrackPFP[NRecoShowers]/I");
        fTree->Branch("RecoShowerRecoNChildShowerPFP",fRecoShowerRecoNChildShowerPFP,"RecoShowerRecoNChildShowerPFP[NRecoShowers]/I");
        fTree->Branch("RecoShowerRecoDirX",fRecoShowerRecoDirX,"RecoShowerRecoDirX[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoDirY",fRecoShowerRecoDirY,"RecoShowerRecoDirY[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecoDirZ",fRecoShowerRecoDirZ,"RecoShowerRecoDirZ[NRecoShowers]/D");
        fTree->Branch("RecoShowerRecodEdx",fRecoShowerRecodEdx,"RecoShowerRecodEdx[NRecoShowers][3]/D");
        fTree->Branch("RecoShowerRecoEnergy",fRecoShowerRecoEnergy,"RecoShowerRecoEnergy[NRecoShowers][3]/D");
        fTree->Branch("RecoShowerRecoBestPlane",&fRecoShowerRecoBestPlane,"RecoShowerRecoBestPlane[3]/I");
        fTree->Branch("RecoShowerRecoLength",&fRecoShowerRecoLength,"RecoShowerRecoLength[3]/D");
        fTree->Branch("RecoShowerRecoOpeningAngle",&fRecoShowerRecoOpeningAngle,"RecoShowerRecoOpeningAngle[3]/D");
        fTree->Branch("RecoShowerRecoIsPrimaryPFPDaughter",fRecoShowerRecoIsPrimaryPFPDaughter,"RecoShowerRecoIsPrimaryPFPDaughter[NRecoShowers]/O");
        fTree->Branch("RecoShowerMVAElectron",fRecoShowerMVAElectron,"RecoShowerMVAElectron[NRecoShowers]/D");
        fTree->Branch("RecoShowerMVAPion",fRecoShowerMVAPion,"RecoShowerMVAPion[NRecoShowers]/D");
        fTree->Branch("RecoShowerMVAMuon",fRecoShowerMVAMuon,"RecoShowerMVAMuon[NRecoShowers]/D");
        fTree->Branch("RecoShowerMVAProton",fRecoShowerMVAProton,"RecoShowerMVAProton[NRecoShowers]/D");
        fTree->Branch("RecoShowerMVAPhoton",fRecoShowerMVAPhoton,"RecoShowerMVAPhoton[NRecoShowers]/D");
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
        fTree->Branch("RecoShowerEnhancedPandrizzleScore",&fRecoShowerEnhancedPandrizzleScore, "RecoShowerEnhancedPandrizzleScore[NRecoShowers]/D"); 
        fTree->Branch("RecoShowerBackupPandrizzleScore",&fRecoShowerBackupPandrizzleScore, "RecoShowerBackupPandrizzleScore[NRecoShowers]/D");
        fTree->Branch("RecoShowerPandrizzleIsFilled",&fRecoShowerPandrizzleIsFilled, "RecoShowerPandrizzleIsFilled[NRecoShowers]/O");

        fTree->Branch("RecoNuVtxNShowers",&fRecoNuVtxNShowers);
        fTree->Branch("RecoNuVtxNTracks",&fRecoNuVtxNTracks);
        fTree->Branch("RecoNuVtxNChildren",&fRecoNuVtxNChildren);

        fTree->Branch("RecoEventCharge",&fRecoEventCharge);
        fTree->Branch("NumuRecoMomLep",&fNumuRecoMomLep);
        fTree->Branch("NumuRecoEHad",&fNumuRecoEHad);
        fTree->Branch("NueRecoMomLep",&fNueRecoMomLep);
        fTree->Branch("NueRecoEHad",&fNueRecoEHad);
    }

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
  art::Handle<sumdata::POTSummary> potListHandle;

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
  fSelTrackTrueEndX = kDefDoub;
  fSelTrackTrueEndY = kDefDoub;
  fSelTrackTrueEndZ = kDefDoub;
  //reco bits
  fSelTrackRecoNHits = kDefInt;
  fSelTrackRecoCompleteness = kDefDoub;
  fSelTrackRecoHitPurity = kDefDoub;
  fSelTrackRecoStartX = kDefDoub;
  fSelTrackRecoStartY = kDefDoub;
  fSelTrackRecoStartZ = kDefDoub;
  fSelTrackRecoEndX = kDefDoub;
  fSelTrackRecoEndY = kDefDoub;
  fSelTrackRecoEndZ = kDefDoub;
  fSelTrackRecoUpstreamX = kDefDoub;
  fSelTrackRecoUpstreamY = kDefDoub;
  fSelTrackRecoUpstreamZ = kDefDoub;
  fSelTrackRecoDownstreamX = kDefDoub;
  fSelTrackRecoDownstreamY = kDefDoub;
  fSelTrackRecoDownstreamZ = kDefDoub;
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
    fRecoTrackTrueEndX[i_recotrack] = kDefDoub;
    fRecoTrackTrueEndY[i_recotrack] = kDefDoub;
    fRecoTrackTrueEndZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoNHits[i_recotrack] = kDefInt;
    fRecoTrackRecoCompleteness[i_recotrack] = kDefDoub;
    fRecoTrackRecoHitPurity[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartX[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartY[i_recotrack] = kDefDoub;
    fRecoTrackRecoStartZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndX[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndY[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamX[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamY[i_recotrack] = kDefDoub;
    fRecoTrackRecoUpstreamZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamX[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamY[i_recotrack] = kDefDoub;
    fRecoTrackRecoDownstreamZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndClosestToVertexX[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndClosestToVertexY[i_recotrack] = kDefDoub;
    fRecoTrackRecoEndClosestToVertexZ[i_recotrack] = kDefDoub;
    fRecoTrackRecoLength[i_recotrack] = kDefDoub;
    fRecoTrackRecoContained[i_recotrack] = kDefInt;
    fRecoTrackRecoMomMethod[i_recotrack] = kDefInt;
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
  fSelShowerTrueEndX = kDefDoub;
  fSelShowerTrueEndY = kDefDoub;
  fSelShowerTrueEndZ = kDefDoub;
  //reco bits
  fSelShowerRecoNHits = kDefInt;
  fSelShowerRecoCompleteness = kDefDoub;
  fSelShowerRecoHitPurity = kDefDoub;
  fSelShowerRecoStartX = kDefDoub;
  fSelShowerRecoStartY = kDefDoub;
  fSelShowerRecoStartZ = kDefDoub;
  fSelShowerRecoCharge = kDefDoub;
  fSelShowerRecoMom = kDefDoub;
  fSelShowerRecoVertexX = kDefDoub;
  fSelShowerRecoVertexY = kDefDoub;
  fSelShowerRecoVertexZ = kDefDoub;
  fSelShowerDistanceToNuVertexU = kDefDoub;
  fSelShowerDistanceToNuVertexV = kDefDoub;
  fSelShowerDistanceToNuVertexW = kDefDoub;
  fSelShowerRecoNChildPFP = kDefInt;
  fSelShowerRecoNChildTrackPFP = kDefInt;
  fSelShowerRecoNChildShowerPFP = kDefInt;
  fSelShowerRecoDirX = kDefDoub;
  fSelShowerRecoDirY = kDefDoub;
  fSelShowerRecoDirZ = kDefDoub;
  fSelShowerRecoBestPlane = kDefInt;
  fSelShowerRecoLength = kDefDoub;
  fSelShowerRecoOpeningAngle = kDefDoub;

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
  fSelShowerEnhancedPandrizzleScore = kDefDoub;
  fSelShowerBackupPandrizzleScore      = kDefDoub;

  fSelShowerPandrizzleIsFilled      = 0;

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
    fRecoShowerTrueEndX[i_shower] = kDefDoub;
    fRecoShowerTrueEndY[i_shower] = kDefDoub;
    fRecoShowerTrueEndZ[i_shower] = kDefDoub;
    fRecoShowerRecoNHits[i_shower] = kDefInt;
    fRecoShowerRecoCompleteness[i_shower] = kDefDoub;
    fRecoShowerRecoHitPurity[i_shower] = kDefDoub;
    fRecoShowerRecoStartX[i_shower] = kDefDoub;
    fRecoShowerRecoStartY[i_shower] = kDefDoub;
    fRecoShowerRecoStartZ[i_shower] = kDefDoub;
    fRecoShowerRecoCharge[i_shower] = kDefDoub;
    fRecoShowerRecoMom[i_shower] = kDefDoub;
    fRecoShowerRecoVertexX[i_shower] = kDefDoub;
    fRecoShowerRecoVertexY[i_shower] = kDefDoub;
    fRecoShowerRecoVertexZ[i_shower] = kDefDoub;
    fRecoShowerDistanceToNuVertexU[i_shower] = kDefDoub;
    fRecoShowerDistanceToNuVertexV[i_shower] = kDefDoub;
    fRecoShowerDistanceToNuVertexW[i_shower] = kDefDoub;
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
    for (int i_plane = 0; i_plane < 3; i_plane++){
        fRecoShowerRecodEdx[i_shower][i_plane] = kDefDoub;
        fRecoShowerRecoEnergy[i_shower][i_plane] = kDefDoub;
    }
    fRecoShowerMVAElectron[i_shower] = kDefDoub;
    fRecoShowerMVAPion[i_shower] = kDefDoub;
    fRecoShowerMVAMuon[i_shower] = kDefDoub;
    fRecoShowerMVAProton[i_shower] = kDefDoub;
    fRecoShowerMVAPhoton[i_shower] = kDefDoub;
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
    fRecoShowerEnhancedPandrizzleScore[i_shower]   = kDefDoub;
    fRecoShowerBackupPandrizzleScore[i_shower]      = kDefDoub;
    fRecoShowerPandrizzleIsFilled[i_shower]      = 0;

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
  try
  {
      std::vector<art::Ptr<recob::Hit>> hitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
      fRecoEventCharge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, hitList);
  }
  catch(...)
  {
      return;
  }

  // Get CVN results
  art::Handle<std::vector<cvn::Result>> cvnResult;
  evt.getByLabel(fCVNModuleLabel, fCVNProductInstance, cvnResult);

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
  art::Handle<std::vector<simb::MCTruth>> mcTruthListHandle;
  std::vector<art::Ptr<simb::MCTruth>> mcList;
  if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle)){
    art::fill_ptr_vector(mcList, mcTruthListHandle);
  }

  if (mcList.size() == 0)
    return;

  if (mcList.size() > 1)
    mf::LogWarning("CCNuSelection") << "There are  " << mcList.size() << " MCTruth in this event.  Only taking the first one.";

  //Get the flux record
  art::Handle<std::vector<simb::MCFlux>> mcFluxListHandle;
  std::vector<art::Ptr<simb::MCFlux>> mcFlux;
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
  const std::vector<art::Ptr<simb::MCParticle>> associated_particles = fmpt.at(mcTruth.key());

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

///////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::FillVertexInfo(art::Event const & evt)
{
  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return;

  art::Ptr<recob::PFParticle> nu_pfp = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nu_pfp, evt, fPFParticleModuleLabel);

  fRecoNuVtxNChildren = nuChildren.size();

  for (art::Ptr<recob::PFParticle> nuChild : nuChildren)
  {
    int pdg = nuChild->PdgCode();

    if (pdg == 11) 
      fRecoNuVtxNShowers++;
    else if (pdg == 13) 
      fRecoNuVtxNTracks++;
  }

  try
  {
      art::Ptr<recob::Vertex> nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nu_pfp, evt, fPFParticleModuleLabel);

      fRecoNuVtxX = nuVertex->position().X();
      fRecoNuVtxY = nuVertex->position().Y();
      fRecoNuVtxZ = nuVertex->position().Z();
  }
  catch(...)
  {
      return;
  }

  return;
}

////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::GetRecoTrackInfo(art::Event const & evt)
{
  if(!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel); 
  std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fPFParticleModuleLabel);

  int trackCounter = 0;

  for (art::Ptr<recob::PFParticle> pfp : pfps)
  {
    if (trackCounter == kDefMaxNRecoTracks)
      break;

    if (!dune_ana::DUNEAnaPFParticleUtils::HasTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
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

    // General track variables
    recob::Track::Point_t trackStart, trackEnd;
    std::tie(trackStart, trackEnd) = current_track->Extent(); 
    fRecoTrackRecoStartX[trackCounter] = trackStart.X();
    fRecoTrackRecoStartY[trackCounter] = trackStart.Y();
    fRecoTrackRecoStartZ[trackCounter] = trackStart.Z();
    fRecoTrackRecoEndX[trackCounter] = trackEnd.X();
    fRecoTrackRecoEndY[trackCounter] = trackEnd.Y();
    fRecoTrackRecoEndZ[trackCounter] = trackEnd.Z();

    if (fRecoTrackRecoEndZ[trackCounter] > fRecoTrackRecoStartZ[trackCounter])
    {
      fRecoTrackRecoUpstreamX[trackCounter] = fRecoTrackRecoStartX[trackCounter];
      fRecoTrackRecoUpstreamY[trackCounter] = fRecoTrackRecoStartY[trackCounter];
      fRecoTrackRecoUpstreamZ[trackCounter] = fRecoTrackRecoStartZ[trackCounter];
      fRecoTrackRecoDownstreamX[trackCounter] = fRecoTrackRecoEndX[trackCounter];
      fRecoTrackRecoDownstreamY[trackCounter] = fRecoTrackRecoEndY[trackCounter];
      fRecoTrackRecoDownstreamZ[trackCounter] = fRecoTrackRecoEndZ[trackCounter];
    }
    else
    {
      fRecoTrackRecoDownstreamX[trackCounter] = fRecoTrackRecoStartX[trackCounter];
      fRecoTrackRecoDownstreamY[trackCounter] = fRecoTrackRecoStartY[trackCounter];
      fRecoTrackRecoDownstreamZ[trackCounter] = fRecoTrackRecoStartZ[trackCounter];
      fRecoTrackRecoUpstreamX[trackCounter] = fRecoTrackRecoEndX[trackCounter];
      fRecoTrackRecoUpstreamY[trackCounter] = fRecoTrackRecoEndY[trackCounter];
      fRecoTrackRecoUpstreamZ[trackCounter] = fRecoTrackRecoEndZ[trackCounter];
    }

    fRecoTrackRecoLength[trackCounter] = current_track->Length();

    try
    {
      art::Ptr<recob::Vertex> track_reco_vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, fPFParticleModuleLabel);

      fRecoTrackRecoVertexX[trackCounter] = track_reco_vertex->position().X();
      fRecoTrackRecoVertexY[trackCounter] = track_reco_vertex->position().Y();
      fRecoTrackRecoVertexZ[trackCounter] = track_reco_vertex->position().Z();
    }
    catch(...)
    {
    }

    TVector3 upstream_end(fRecoTrackRecoUpstreamX[trackCounter], fRecoTrackRecoUpstreamY[trackCounter], fRecoTrackRecoUpstreamZ[trackCounter]);
    TVector3 downstream_end(fRecoTrackRecoDownstreamX[trackCounter], fRecoTrackRecoDownstreamY[trackCounter], fRecoTrackRecoDownstreamZ[trackCounter]);
    TVector3 vertex_pos(fRecoTrackRecoVertexX[trackCounter], fRecoTrackRecoVertexY[trackCounter], fRecoTrackRecoVertexZ[trackCounter]);

    if ((vertex_pos - upstream_end).Mag() < (vertex_pos - downstream_end).Mag())
    {
      fRecoTrackRecoEndClosestToVertexX[trackCounter] = fRecoTrackRecoUpstreamX[trackCounter];
      fRecoTrackRecoEndClosestToVertexY[trackCounter] = fRecoTrackRecoUpstreamY[trackCounter];
      fRecoTrackRecoEndClosestToVertexZ[trackCounter] = fRecoTrackRecoUpstreamZ[trackCounter];
    }
    else
    {
      fRecoTrackRecoEndClosestToVertexX[trackCounter] = fRecoTrackRecoDownstreamX[trackCounter];
      fRecoTrackRecoEndClosestToVertexY[trackCounter] = fRecoTrackRecoDownstreamY[trackCounter];
      fRecoTrackRecoEndClosestToVertexZ[trackCounter] = fRecoTrackRecoDownstreamZ[trackCounter];
    }

    // Child particle info
    FillChildPFPInformation(pfp, evt, fRecoTrackRecoNChildPFP[trackCounter], fRecoTrackRecoNChildTrackPFP[trackCounter], fRecoTrackRecoNChildShowerPFP[trackCounter]);

    // Fill momentum variables
    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(current_track, evt)));
    fRecoTrackRecoContained[trackCounter] = energyRecoHandle->longestTrackContained;
    fRecoTrackRecoMomMethod[trackCounter] = energyRecoHandle->trackMomMethod;

    if (energyRecoHandle->trackMomMethod == 1)
    {   
      // momentum by range was used to calculate ENu
      fRecoTrackRecoMomContained[trackCounter] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    }
    else if (energyRecoHandle->trackMomMethod == 0)
    {
      // momentum by MCS
      fRecoTrackRecoMomMCS[trackCounter] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    }

    std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, current_track_hits, 1);
    fRecoTrackRecoCompleteness[trackCounter] = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, current_track_hits, eventHitList, g4id);
    fRecoTrackRecoHitPurity[trackCounter] = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, current_track_hits, g4id);

    if (TruthMatchUtils::Valid(g4id)){
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);

        if (matched_mcparticle)
        {
          fRecoTrackTruePDG[trackCounter] = matched_mcparticle->PdgCode();

          if (matched_mcparticle->Mother() == 0) 
            fRecoTrackTruePrimary[trackCounter] = true;
          else 
            fRecoTrackTruePrimary[trackCounter] = false;

          fRecoTrackTrueMomX[trackCounter] = matched_mcparticle->Momentum().X();
          fRecoTrackTrueMomY[trackCounter] = matched_mcparticle->Momentum().Y();
          fRecoTrackTrueMomZ[trackCounter] = matched_mcparticle->Momentum().Z();
          fRecoTrackTrueStartX[trackCounter] = matched_mcparticle->Position(0).X();
          fRecoTrackTrueStartY[trackCounter] = matched_mcparticle->Position(0).Y();
          fRecoTrackTrueStartZ[trackCounter] = matched_mcparticle->Position(0).Z();
          fRecoTrackTrueEndX[trackCounter] = matched_mcparticle->EndPosition().X();
          fRecoTrackTrueEndY[trackCounter] = matched_mcparticle->EndPosition().Y();
          fRecoTrackTrueEndZ[trackCounter] = matched_mcparticle->EndPosition().Z();
        }
    }

    // MVA PID
    art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Track>>{current_track}, evt, fPIDModuleLabel);
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

void FDSelection::CCNuSelection::FillChildPFPInformation(art::Ptr<recob::PFParticle> const pfp, art::Event const & evt, int &n_child_pfp, int &n_child_track_pfp, int &n_child_shower_pfp)
{
  n_child_pfp = 0;
  n_child_track_pfp = 0;
  n_child_shower_pfp = 0;

  std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfp, evt, fPFParticleModuleLabel);

  for (art::Ptr<recob::PFParticle> childPFP : childPFPs)
  {
    int pdg = childPFP->PdgCode();

    if (pdg == 13)
      n_child_track_pfp++;
   else if (pdg == 11)
      n_child_shower_pfp++;
   else 
      std::cout << "FillChildPFPInformation: found a child PFP with an unexpected pdg code: " << pdg << std::endl;
  }

  return;
}

///////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunTrackSelection(art::Event const & evt)
{
  // Get the selected track
  art::Ptr<recob::Track> sel_track = fRecoTrackSelector->FindSelectedTrack(evt);

  if (!sel_track.isAvailable()) 
  {
    std::cout<<"FDSelection::CCNuSelection::RunTrackSelection - no track returned from selection" << std::endl; 
    return;
  }

  std::vector<art::Ptr<recob::Hit>> sel_track_hits = dune_ana::DUNEAnaTrackUtils::GetHits(sel_track, evt, fTrackModuleLabel);
  fSelTrackRecoNHits = sel_track_hits.size();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  fSelTrackRecoCharge  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, sel_track_hits); 

  // Get general track variables
  recob::Track::Point_t trackStart, trackEnd;
  std::tie(trackStart, trackEnd) = sel_track->Extent(); 
  fSelTrackRecoStartX = trackStart.X();
  fSelTrackRecoStartY = trackStart.Y();
  fSelTrackRecoStartZ = trackStart.Z();
  fSelTrackRecoEndX = trackEnd.X();
  fSelTrackRecoEndY = trackEnd.Y();
  fSelTrackRecoEndZ = trackEnd.Z();

  if (fSelTrackRecoEndZ > fSelTrackRecoStartZ)
  {
    fSelTrackRecoUpstreamX = fSelTrackRecoStartX;
    fSelTrackRecoUpstreamY = fSelTrackRecoStartY;
    fSelTrackRecoUpstreamZ = fSelTrackRecoStartZ;
    fSelTrackRecoDownstreamX = fSelTrackRecoEndX;
    fSelTrackRecoDownstreamY = fSelTrackRecoEndY;
    fSelTrackRecoDownstreamZ = fSelTrackRecoEndZ;
  }
  else
  {
    fSelTrackRecoDownstreamX = fSelTrackRecoStartX;
    fSelTrackRecoDownstreamY = fSelTrackRecoStartY;
    fSelTrackRecoDownstreamZ = fSelTrackRecoStartZ;
    fSelTrackRecoUpstreamX = fSelTrackRecoEndX;
    fSelTrackRecoUpstreamY = fSelTrackRecoEndY;
    fSelTrackRecoUpstreamZ = fSelTrackRecoEndZ;
  }

  fSelTrackRecoLength = sel_track->Length();

  art::Ptr<recob::PFParticle> sel_track_pfp = dune_ana::DUNEAnaTrackUtils::GetPFParticle(sel_track, evt, fTrackModuleLabel);

  try
  {
      art::Ptr<recob::Vertex> track_reco_vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(sel_track_pfp, evt, fPFParticleModuleLabel);

      fSelTrackRecoVertexX = track_reco_vertex->position().X();
      fSelTrackRecoVertexY = track_reco_vertex->position().Y();
      fSelTrackRecoVertexZ = track_reco_vertex->position().Z();
  }
  catch(...)
  {
  }

  TVector3 upstream_end(fSelTrackRecoUpstreamX, fSelTrackRecoUpstreamY, fSelTrackRecoUpstreamZ);
  TVector3 downstream_end(fSelTrackRecoDownstreamX, fSelTrackRecoDownstreamY, fSelTrackRecoDownstreamZ);
  TVector3 vertex_pos(fSelTrackRecoVertexX, fSelTrackRecoVertexY, fSelTrackRecoVertexZ);

  if ((vertex_pos - upstream_end).Mag() < (vertex_pos - downstream_end).Mag())
  {
    fSelTrackRecoEndClosestToVertexX = fSelTrackRecoUpstreamX;
    fSelTrackRecoEndClosestToVertexY = fSelTrackRecoUpstreamY;
    fSelTrackRecoEndClosestToVertexZ = fSelTrackRecoUpstreamZ;
  }
  else
  {
    fSelTrackRecoEndClosestToVertexX = fSelTrackRecoDownstreamX;
    fSelTrackRecoEndClosestToVertexY = fSelTrackRecoDownstreamY;
    fSelTrackRecoEndClosestToVertexZ = fSelTrackRecoDownstreamZ;
  }

  FillChildPFPInformation(sel_track_pfp, evt, fSelTrackRecoNChildPFP, fSelTrackRecoNChildTrackPFP, fSelTrackRecoNChildShowerPFP);

  // Fill neutrino energy variables
  // Use selected track to get neutrino energy
  std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(sel_track, evt)));
  fNumuRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNumuRecoEHad = energyRecoHandle->fHadLorentzVector.E();

  fSelTrackRecoContained = energyRecoHandle->longestTrackContained; 
  fSelTrackRecoMomMethod = energyRecoHandle->trackMomMethod;

  if (energyRecoHandle->trackMomMethod == 1)
  { 
    // momentum by range was used to calculate ENu
    fSelTrackRecoMomContained = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fNumuRecoMomLep = fSelTrackRecoMomContained;
  }
  else if (energyRecoHandle->trackMomMethod == 0)
  {
    // momentum by MCS
    fSelTrackRecoMomMCS = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());
    fNumuRecoMomLep = fSelTrackRecoMomMCS;
  }

  // Get truth information
  int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sel_track_hits, 1);
  std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
  fSelTrackRecoCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, sel_track_hits, eventHitList, g4id);
  fSelTrackRecoHitPurity = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, sel_track_hits, g4id);

  if (TruthMatchUtils::Valid(g4id))
  {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
      if (matched_mcparticle)
      {
        fSelTrackTruePDG = matched_mcparticle->PdgCode();

        if (matched_mcparticle->Mother() == 0) 
          fSelTrackTruePrimary = 1;
        else 
          fSelTrackTruePrimary = 0;

        fSelTrackTrueMomX = matched_mcparticle->Momentum().X();
        fSelTrackTrueMomY = matched_mcparticle->Momentum().Y();
        fSelTrackTrueMomZ = matched_mcparticle->Momentum().Z();
        fSelTrackTrueMomT = matched_mcparticle->Momentum().T();
        fSelTrackTrueStartX = matched_mcparticle->Position(0).X();
        fSelTrackTrueStartY = matched_mcparticle->Position(0).Y();
        fSelTrackTrueStartZ = matched_mcparticle->Position(0).Z();
        fSelTrackTrueEndX = matched_mcparticle->EndPosition().X();
        fSelTrackTrueEndY = matched_mcparticle->EndPosition().Y();
        fSelTrackTrueEndZ = matched_mcparticle->EndPosition().Z();
      }
  }

  // MVA PID
  art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Track>>{sel_track}, evt, fPIDModuleLabel);
  art::Ptr<anab::MVAPIDResult> pid(findPIDResult.at(0));

  if (pid.isAvailable())
  {
    std::map<std::string,double> mvaOutMap = pid->mvaOutput;

    if (!mvaOutMap.empty())
    {
      fSelTrackMVAElectron = mvaOutMap["electron"];
      fSelTrackMVAPion = mvaOutMap["pich"];
      fSelTrackMVAMuon = mvaOutMap["muon"];
      fSelTrackMVAProton = mvaOutMap["proton"];
      fSelTrackMVAPhoton = mvaOutMap["photon"];
    }
  }

  // Pandizzle variables
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

  // DeepPan variables
  ctp::CTPResult deepPanPIDResult = fConvTrackPID.RunConvolutionalTrackPID(sel_track_pfp, evt);
  if (deepPanPIDResult.IsValid())
  {
      fSelTrackDeepPanMuVar = deepPanPIDResult.GetMuonScore();
      fSelTrackDeepPanPiVar = deepPanPIDResult.GetPionScore();
      fSelTrackDeepPanProtonVar = deepPanPIDResult.GetProtonScore();
  }
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::GetRecoShowerInfo(art::Event const & evt)
{
  if(!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel); 
  std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fPFParticleModuleLabel);

  int showerCounter = 0;

  for (art::Ptr<recob::PFParticle> pfp : pfps)
  {
    if (showerCounter == kDefMaxNRecoShowers)
      break;

    if (!dune_ana::DUNEAnaPFParticleUtils::HasShower(pfp, evt, fPFParticleModuleLabel, fShowerModuleLabel))
      continue;

    for (art::Ptr<recob::PFParticle> nuChild : nuChildren)
    {
      if (nuChild->Self() == pfp->Self())
      {
        fRecoShowerRecoIsPrimaryPFPDaughter[showerCounter] = true;
        break;
      }
    }

    const std::vector<art::Ptr<recob::Hit>> current_shower_hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPFParticleModuleLabel);
    fRecoShowerRecoNHits[showerCounter] = current_shower_hits.size();

    art::Ptr<recob::Shower> current_shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, fPFParticleModuleLabel, fShowerModuleLabel);

    // General shower variables
    fRecoShowerRecoDirX[showerCounter] = current_shower->Direction().X();
    fRecoShowerRecoDirY[showerCounter] = current_shower->Direction().Y();
    fRecoShowerRecoDirZ[showerCounter] = current_shower->Direction().Z();
    fRecoShowerRecoStartX[showerCounter] = current_shower->ShowerStart().X();
    fRecoShowerRecoStartY[showerCounter] = current_shower->ShowerStart().Y();
    fRecoShowerRecoStartZ[showerCounter] = current_shower->ShowerStart().Z();
    fRecoShowerRecoBestPlane[showerCounter] = current_shower->best_plane();
    fRecoShowerRecoLength[showerCounter] = current_shower->Length();
    fRecoShowerRecoOpeningAngle[showerCounter] = current_shower->OpenAngle();

    // Vertex info
    try
    {
        art::Ptr<recob::Vertex> shower_reco_vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, fPFParticleModuleLabel);

        fRecoShowerRecoVertexX[showerCounter] = shower_reco_vertex->position().X();
        fRecoShowerRecoVertexY[showerCounter] = shower_reco_vertex->position().Y();
        fRecoShowerRecoVertexZ[showerCounter] = shower_reco_vertex->position().Z();
    }
    catch(...)
    {
    }

    try
    {
        std::vector<art::Ptr<recob::Hit>> hitListU = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(current_shower_hits, 0);
        fRecoShowerDistanceToNuVertexU[showerCounter] = DistanceToNuVertex(evt, hitListU, TVector3(fNuX, fNuY, fNuZ));
    }
    catch(...)
    {
    } 

    try
    {
        std::vector<art::Ptr<recob::Hit>> hitListV = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(current_shower_hits, 1);
        fRecoShowerDistanceToNuVertexV[showerCounter] = DistanceToNuVertex(evt, hitListV, TVector3(fNuX, fNuY, fNuZ));
    }
    catch(...)
    {
    }

    try
    {
        std::vector<art::Ptr<recob::Hit>> hitListW = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(current_shower_hits, 2);
        fRecoShowerDistanceToNuVertexW[showerCounter] = DistanceToNuVertex(evt, hitListW, TVector3(fNuX, fNuY, fNuZ));
    }
    catch(...)
    {
    }

    // Momentum and energy
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    fRecoShowerRecoCharge[showerCounter]  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, current_shower_hits); 

    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(current_shower, evt)));
    fRecoShowerRecoMom[showerCounter] = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

    if (current_shower->dEdx().size() > 0)
    {
      for (int i_plane = 0; i_plane < 3; i_plane++)
      {
        fRecoShowerRecodEdx[showerCounter][i_plane] = current_shower->dEdx()[i_plane];
        fRecoShowerRecoEnergy[showerCounter][i_plane] = current_shower->Energy()[i_plane];
      }
    }

    FillChildPFPInformation(pfp, evt, fRecoShowerRecoNChildPFP[showerCounter], fRecoShowerRecoNChildTrackPFP[showerCounter], fRecoShowerRecoNChildShowerPFP[showerCounter]);

    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, current_shower_hits, 1);
    std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
    fRecoShowerRecoCompleteness[showerCounter] = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, current_shower_hits, eventHitList, g4id);
    fRecoShowerRecoHitPurity[showerCounter] = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, current_shower_hits, g4id);

    if (TruthMatchUtils::Valid(g4id))
    {
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        const simb::MCParticle* matched_mcparticle = pi_serv->TrackIdToParticle_P(g4id);
        if (matched_mcparticle)
        {
          fRecoShowerTruePDG[showerCounter] = matched_mcparticle->PdgCode();

          if 
            (matched_mcparticle->Mother()==0) fRecoShowerTruePrimary[showerCounter] = 1;
          else 
            fRecoShowerTruePrimary[showerCounter] = 0;

          fRecoShowerTrueMomX[showerCounter] = matched_mcparticle->Momentum().X();
          fRecoShowerTrueMomY[showerCounter] = matched_mcparticle->Momentum().Y();
          fRecoShowerTrueMomZ[showerCounter] = matched_mcparticle->Momentum().Z();
          fRecoShowerTrueMomT[showerCounter] = matched_mcparticle->Momentum().T();
          fRecoShowerTrueStartX[showerCounter] = matched_mcparticle->Position(0).X();
          fRecoShowerTrueStartY[showerCounter] = matched_mcparticle->Position(0).Y();
          fRecoShowerTrueStartZ[showerCounter] = matched_mcparticle->Position(0).Z();
          fRecoShowerTrueEndX[showerCounter] = matched_mcparticle->EndPosition().X();
          fRecoShowerTrueEndY[showerCounter] = matched_mcparticle->EndPosition().Y();
          fRecoShowerTrueEndZ[showerCounter] = matched_mcparticle->EndPosition().Z();
        }
    }

    // MVA PID
    art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower>>{current_shower}, evt, fPIDModuleLabel);
    art::Ptr<anab::MVAPIDResult> pid(findPIDResult.at(0));

    if (pid.isAvailable())
    {
      std::map<std::string,double> mvaOutMap = pid->mvaOutput;

      if (!mvaOutMap.empty())
      {
        fRecoShowerMVAElectron[showerCounter] = mvaOutMap["electron"];
        fRecoShowerMVAPion[showerCounter] = mvaOutMap["pich"];
        fRecoShowerMVAMuon[showerCounter] = mvaOutMap["muon"];
        fRecoShowerMVAProton[showerCounter] = mvaOutMap["proton"];
        fRecoShowerMVAPhoton[showerCounter] = mvaOutMap["photon"];
      }
    }

    // Pandrizzle
    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(current_shower, evt));

    fRecoShowerPandrizzleEvalRatio[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
    fRecoShowerPandrizzleConcentration[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
    fRecoShowerPandrizzleCoreHaloRatio[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
    fRecoShowerPandrizzleConicalness[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
    fRecoShowerPandrizzledEdxBestPlane[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
    fRecoShowerPandrizzleDisplacement[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
    fRecoShowerPandrizzleDCA[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
    fRecoShowerPandrizzleWideness[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
    fRecoShowerPandrizzleEnergyDensity[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);
    fRecoShowerPandrizzlePathwayLengthMin[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kPathwayLengthMin);
    fRecoShowerPandrizzleMaxShowerStartPathwayScatteringAngle2D[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxShowerStartPathwayScatteringAngle2D);
    fRecoShowerPandrizzleMaxNPostShowerStartHits[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxNPostShowerStartHits);
    fRecoShowerPandrizzleMaxPostShowerStartScatterAngle[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartScatterAngle);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyAsymmetry[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartShowerStartEnergyAsymmetry[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartShowerStartEnergyAsymmetry);
    fRecoShowerPandrizzleMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance[showerCounter] = 
      pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
    fRecoShowerPandrizzleMinPostShowerStartShowerStartMoliereRadius[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinPostShowerStartShowerStartMoliereRadius);
    fRecoShowerPandrizzleMaxPostShowerStartOpeningAngle[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartOpeningAngle);
    fRecoShowerPandrizzleMaxFoundHitRatio[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxFoundHitRatio);
    fRecoShowerPandrizzleMaxInitialGapSize[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxInitialGapSize);
    fRecoShowerPandrizzleMinLargestProjectedGapSize[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinLargestProjectedGapSize);
    fRecoShowerPandrizzleNViewsWithAmbiguousHits[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kNViewsWithAmbiguousHits);
    fRecoShowerPandrizzleAmbiguousHitMaxUnaccountedEnergy[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kAmbiguousHitMaxUnaccountedEnergy);

    fRecoShowerPandrizzleBDTMethod[showerCounter] = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);
    float pandrizzleScore(pandrizzleRecord.GetMVAScore());
    fRecoShowerBackupPandrizzleScore[showerCounter] = (std::fabs(fRecoShowerPandrizzleBDTMethod[showerCounter] - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerEnhancedPandrizzleScore[showerCounter] = (std::fabs(fRecoShowerPandrizzleBDTMethod[showerCounter] - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
    fRecoShowerPandrizzleIsFilled[showerCounter] = pandrizzleRecord.IsFilled();

    ++showerCounter;
  }

  fNRecoShowers = showerCounter;

  return;
}

///////////////////////////////////////////////////////////////////////////////////

void FDSelection::CCNuSelection::RunShowerSelection(art::Event const & evt)
{
  // Get the selected shower
  art::Ptr<recob::Shower> sel_shower = fRecoShowerSelector->FindSelectedShower(evt);

  if (!sel_shower.isAvailable()) 
  {
    std::cout << "FDSelection::CCNuSelection::RunShowerSelection - no shower selected by tool" << std::endl;
    return;
  }

  art::Ptr<recob::PFParticle> sel_pfp = dune_ana::DUNEAnaShowerUtils::GetPFParticle(sel_shower, evt, fShowerModuleLabel);

  //Get the hits for said shower
  const std::vector<art::Ptr<recob::Hit>> sel_shower_hits = dune_ana::DUNEAnaShowerUtils::GetHits(sel_shower, evt, fShowerModuleLabel);
  fSelShowerRecoNHits = sel_shower_hits.size();

  fSelShowerRecoDirX = sel_shower->Direction().X();
  fSelShowerRecoDirY = sel_shower->Direction().Y();
  fSelShowerRecoDirZ = sel_shower->Direction().Z();
  fSelShowerRecoStartX = sel_shower->ShowerStart().X();
  fSelShowerRecoStartY = sel_shower->ShowerStart().Y();
  fSelShowerRecoStartZ = sel_shower->ShowerStart().Z();
  fSelShowerRecoBestPlane = sel_shower->best_plane();
  fSelShowerRecoLength = sel_shower->Length();
  fSelShowerRecoOpeningAngle = sel_shower->OpenAngle();

  // Vertex info
  try
  {
    art::Ptr<recob::Vertex> shower_reco_vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(sel_pfp, evt, fPFParticleModuleLabel);

    fSelShowerRecoVertexX = shower_reco_vertex->position().X();
    fSelShowerRecoVertexY = shower_reco_vertex->position().Y();
    fSelShowerRecoVertexZ = shower_reco_vertex->position().Z();
  }
  catch(...)
  {
  }

  try
  {
    std::vector<art::Ptr<recob::Hit>> hitListU = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(sel_shower_hits, 0);
    fSelShowerDistanceToNuVertexU = DistanceToNuVertex(evt, hitListU, TVector3(fNuX, fNuY, fNuZ));
  }
  catch(...)
  {
  }

  try
  {
    std::vector<art::Ptr<recob::Hit>> hitListV = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(sel_shower_hits, 1);
    fSelShowerDistanceToNuVertexV = DistanceToNuVertex(evt, hitListV, TVector3(fNuX, fNuY, fNuZ));
  }
  catch(...)
  {
  }

  try
  {
    std::vector<art::Ptr<recob::Hit>> hitListW = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(sel_shower_hits, 2);
    fSelShowerDistanceToNuVertexW = DistanceToNuVertex(evt, hitListW, TVector3(fNuX, fNuY, fNuZ));
  }
  catch(...)
  {
  }

  // Momentum and energy
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  fSelShowerRecoCharge  = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, sel_shower_hits); 

  std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(sel_shower, evt)));
  fNueRecoENu = energyRecoHandle->fNuLorentzVector.E();
  fNueRecoEHad = energyRecoHandle->fHadLorentzVector.E();
  fNueRecoMomLep = sqrt(energyRecoHandle->fLepLorentzVector.Vect().Mag2());

  if (sel_shower->dEdx().size() > 0)
  {
    for (unsigned int i_plane = 0; i_plane < 3; i_plane++)
    {
      fSelShowerRecodEdx[i_plane] = sel_shower->dEdx()[i_plane];
      fSelShowerRecoEnergy[i_plane] = sel_shower->Energy()[i_plane];
    }
  }

  // Child PFP info
  FillChildPFPInformation(sel_pfp, evt, fSelShowerRecoNChildPFP, fSelShowerRecoNChildTrackPFP, fSelShowerRecoNChildShowerPFP);

  // Truth info
  int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sel_shower_hits, 1);
  std::vector<art::Ptr<recob::Hit>> eventHitList = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsModuleLabel);
  fSelShowerRecoCompleteness = FDSelectionUtils::CompletenessFromTrueParticleID(clockData, sel_shower_hits, eventHitList, g4id);
  fSelShowerRecoHitPurity = FDSelectionUtils::HitPurityFromTrueParticleID(clockData, sel_shower_hits, g4id);

  if (TruthMatchUtils::Valid(g4id))
  {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      const simb::MCParticle* matched_mcparticle = pi_serv->TrackIdToParticle_P(g4id);
      if (matched_mcparticle)
      {
        fSelShowerTruePDG = matched_mcparticle->PdgCode();

        if (matched_mcparticle->Mother() == 0) 
          fSelShowerTruePrimary = 1;
        else 
          fSelShowerTruePrimary = 0;

        fSelShowerTrueMomX = matched_mcparticle->Momentum().X();
        fSelShowerTrueMomY = matched_mcparticle->Momentum().Y();
        fSelShowerTrueMomZ = matched_mcparticle->Momentum().Z();
        fSelShowerTrueMomT = matched_mcparticle->Momentum().T();
        fSelShowerTrueStartX = matched_mcparticle->Position(0).X();
        fSelShowerTrueStartY = matched_mcparticle->Position(0).Y();
        fSelShowerTrueStartZ = matched_mcparticle->Position(0).Z();
        fSelShowerTrueEndX = matched_mcparticle->EndPosition().X();
        fSelShowerTrueEndY = matched_mcparticle->EndPosition().Y();
        fSelShowerTrueEndZ = matched_mcparticle->EndPosition().Z();
      }
  }

  // MVA PID
  art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower>>{sel_shower}, evt, fPIDModuleLabel);
  art::Ptr<anab::MVAPIDResult> pid(findPIDResult.at(0));

  if (pid.isAvailable())
  {
    std::map<std::string,double> mvaOutMap = pid->mvaOutput;

    if (!mvaOutMap.empty())
    {
      fSelShowerMVAElectron = mvaOutMap["electron"];
      fSelShowerMVAPion = mvaOutMap["pich"];
      fSelShowerMVAMuon = mvaOutMap["muon"];
      fSelShowerMVAProton = mvaOutMap["proton"];
      fSelShowerMVAPhoton = mvaOutMap["photon"];
    }
  }

  //Pandrizzle
  FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(sel_shower, evt));
  fSelShowerPandrizzleEvalRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
  fSelShowerPandrizzleConcentration = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
  fSelShowerPandrizzleCoreHaloRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
  fSelShowerPandrizzleConicalness = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
  fSelShowerPandrizzledEdxBestPlane = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
  fSelShowerPandrizzleDisplacement = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
  fSelShowerPandrizzleDCA = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
  fSelShowerPandrizzleWideness = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
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
  fSelShowerBackupPandrizzleScore = (std::fabs(fSelShowerPandrizzleBDTMethod - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
  fSelShowerEnhancedPandrizzleScore = (std::fabs(fSelShowerPandrizzleBDTMethod - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
  fSelShowerPandrizzleIsFilled = pandrizzleRecord.IsFilled();
}

//////////////////////////////////////////////////////////////////

TVector3 FDSelection::CCNuSelection::ProjectVectorOntoPlane(TVector3 vector_to_project, TVector3 plane_norm_vector){
  TVector3 projected_vector = vector_to_project - (vector_to_project.Dot(plane_norm_vector) / (plane_norm_vector.Mag() * plane_norm_vector.Mag()))*plane_norm_vector;
  return projected_vector;
}

//////////////////////////////////////////////////////////////////

double FDSelection::CCNuSelection::DistanceToNuVertex(art::Event const & evt, std::vector<art::Ptr<recob::Hit>> const artHitList, const TVector3 &nuVertex)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  double closestDistanceSquared(std::numeric_limits<double>::max());

  for (art::Ptr<recob::Hit> const& hit : artHitList)
  {
    const geo::WireID hit_WireID(hit->WireID());
    const double hit_Time(hit->PeakTime());
    const geo::View_t hit_View(hit->View());
    const geo::CryostatID cryostatID(hit_WireID.Cryostat);
    const geo::View_t pandora_View(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));

    TVector3 vertexXYZ = nuVertex;
    geo::Point_t hitXYZ = wireReadout.Wire(hit_WireID).GetCenter();
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
