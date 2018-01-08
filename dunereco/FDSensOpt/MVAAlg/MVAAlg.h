/////////////////////////////////////////////////////////////////
//  \file MVAAlg.h
//  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////
#ifndef MVAAlg_H
#define MVAAlg_H
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <memory>
#include <utility>

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "TMVA/Reader.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TFile.h"


#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RawData/RawDigit.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"



constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxShower     = 1000;  //maximum number of showers
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxVertices   = 1000;  //max number of 3D vertices
constexpr int kMaxPrimaries  = 20000; //maximum number of primary particles
constexpr int kMaxFlash      = 1000;  //maximum number of flashes
//constexpr int kMaxTrackHits  = 1000;  //maximum number of track trajectory points

namespace dunemva{



  //--------------------------------------------------------------- 
  class MVAAlg {
    public:

      MVAAlg(fhicl::ParameterSet const& p);

      virtual ~MVAAlg();

      void reconfigure(fhicl::ParameterSet const& p);

      void Run(const art::Event& evt,std::vector<double>& result, double& wgt);
      void Run(const art::Event& evt,double& result, double& wgt);
      void endSubRun(const art::SubRun& sr);

      float Norm(int ccnc, int nu0, int nu1, int subrun);
      float OscPro(int ccnc, int nu0, int nu1, float NuE);


      // from root script:

      // all tmva reader variables must be float 
      int   itype;
      float fOscPro;
      float weight;
      float evtcharge;
      float ntrack;
      float avgtrklength;
      float maxtrklength;
      float trkdedx;
      float trkrch; // (track charge)/(path charge)... path width must be set in time ticks
      float trkrt;  // (lower half of charge)/(higher half of charge)... ~1 for tracks
      float trkfr;  // (track charge)/(event charge)
      float trkpida_save;
      float nshower;
      float showerdedx;
      float eshower;
      float frshower;
      float nhitspershw;
      float shwlength;
      float shwmax;
      float fract_5_wires;
      float fract_10_wires;
      float fract_50_wires;
      float fract_100_wires;
      float shwdis;
      float shwdisx;
      float shwdisy;
      float shwdisz;
      float shwcosx;
      float shwcosy;
      float shwcosz;
      float trkcosx;
      float trkcosy;
      float trkcosz;
      float ET;

      // Additional variables
      double rawcharge;
      double wirecharge;

    private:

      void  PrepareEvent(const art::Event& event);
      void  MakeTree();
      void  CalculateInputs();

      std::ofstream fFile;

      TMVA::Reader fReader;
      std::string fMVAMethod;
      std::string fWeightFile;

      art::ServiceHandle<geo::Geometry> fGeom;
      art::ServiceHandle<art::TFileService> tfs;

      std::string fSelect;
      std::string fBeamMode;

      // ~~~~~~~~~~~~~~ from NuEAna: ~~~~~~~~~~~~~~~~

      void ResetVars();
      bool insideFidVol(const double posX, const double posY, const double posZ);

      // Declare member data here.
      TTree *fTree;
      TTree* fPOT;

      // TTree variables

      // Run information
      int run;
      int subrun;
      int event;
      float evttime;
      float taulife;
      short isdata;
      double pot;

      // Track information
      int ntracks_reco;                   //number of reconstructed tracks
      int   trkid[kMaxTrack];             //track id from recob::Track::ID(), this does not have to be the index of track
      float trkstartx[kMaxTrack];         //track start position (cm)
      float trkstarty[kMaxTrack];
      float trkstartz[kMaxTrack];
      float trkendx[kMaxTrack];           //track end position (cm)
      float trkendy[kMaxTrack];
      float trkendz[kMaxTrack];
      float trkstartdcosx[kMaxTrack];     //track start direction cosine
      float trkstartdcosy[kMaxTrack];
      float trkstartdcosz[kMaxTrack];
      float trkenddcosx[kMaxTrack];       //track end direction cosine
      float trkenddcosy[kMaxTrack];
      float trkenddcosz[kMaxTrack];
      float trklen[kMaxTrack];            //track length (cm)
      float trkke[kMaxTrack][3];          //track kinetic energy (in 3 planes)
      float trkpida[kMaxTrack][3];        //track PIDA (in 3 planes)
      int   trkbestplane[kMaxTrack];      //best plane for trkke and trkpida

      //geant information for the track
      int   trkg4id[kMaxTrack];           //geant track id for the track
      int   trkg4pdg[kMaxTrack];          //pdg of geant particle
      float trkg4startx[kMaxTrack];       //start position of geant particle
      float trkg4starty[kMaxTrack];
      float trkg4startz[kMaxTrack];
      float trkg4initdedx[kMaxTrack];     //initial dE/dx of the track using true energy (MeV/cm)

      int nhits;
      int nhits_stored;
      Short_t  hit_plane[kMaxHits];      //plane number
      Short_t  hit_wire[kMaxHits];       //wire number
      Int_t    hit_channel[kMaxHits];    //channel ID
      Short_t  hit_tpc[kMaxHits];        //tpc
      Float_t  hit_peakT[kMaxHits];      //peak time
      Float_t  hit_charge[kMaxHits];     //charge (area)
      Float_t  hit_summedADC[kMaxHits];  //summed ADC
      Float_t  hit_startT[kMaxHits];     //hit start time
      Float_t  hit_endT[kMaxHits];       //hit end time
      Int_t    hit_trkkey[kMaxHits];     //track index if hit is associated with a track
      Float_t  hit_dQds[kMaxHits];       //hit dQ/ds
      Float_t  hit_dEds[kMaxHits];       //hit dE/ds
      Float_t  hit_resrange[kMaxHits];   //hit residual range
      Int_t    hit_shwkey[kMaxHits];     //shower index if hit is associated with a shower

      // vertex information
      int infidvol;
      Short_t  nvtx;                     //number of vertices
      Float_t  vtx[kMaxVertices][3];     //vtx[3] 

      Float_t	vtxrecomc;		// distance between mc and reco vtx
      Float_t	vtxrecomcx;		
      Float_t	vtxrecomcy;		
      Float_t	vtxrecomcz;		

      // shower information
      int nshws;                         //number of showers
      int shwid[kMaxShower];             //recob::Shower::ID()
      Float_t shwdcosx[kMaxShower];      //shower direction cosine
      Float_t shwdcosy[kMaxShower];
      Float_t shwdcosz[kMaxShower];
      Float_t shwstartx[kMaxShower];     //shower start position (cm)
      Float_t shwstarty[kMaxShower];
      Float_t shwstartz[kMaxShower];
      Float_t shwenergy[kMaxShower][3];  //shower energy measured on the 3 planes (GeV)
      Float_t shwdedx[kMaxShower][3];    //shower dE/dx of the initial track measured on the 3 plane (MeV/cm)
      int shwbestplane[kMaxShower];      //recommended plane for energy and dE/dx information
      int   shwg4id[kMaxTrack];          //geant track id for the shower

      // flash information
      int    flash_total;                //total number of flashes
      Float_t flash_time[kMaxFlash];     //flash time
      Float_t flash_width[kMaxFlash];    //flash width
      Float_t flash_abstime[kMaxFlash];  //flash absolute time
      Float_t flash_YCenter[kMaxFlash];  //flash y center (cm)
      Float_t flash_YWidth[kMaxFlash];   //flash y width (cm)
      Float_t flash_ZCenter[kMaxFlash];  //flash z center (cm)
      Float_t flash_ZWidth[kMaxFlash];   //flash z width (cm)
      Float_t flash_TotalPE[kMaxFlash];  //flash total pe

      //mctruth information
      Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
      Int_t     nuPDG_truth;     //neutrino PDG code
      Int_t     ccnc_truth;      //0=CC 1=NC
      Int_t     mode_truth;      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      Float_t  enu_truth;       //true neutrino energy
      Float_t  Q2_truth;        //Momentum transfer squared
      Float_t  W_truth;         //hadronic invariant mass
      Float_t  X_truth;
      Float_t  Y_truth;
      Int_t    hitnuc_truth;    //hit nucleon
      Int_t    target_truth;    //hit nucleus 
      Float_t  nuvtxx_truth;    //neutrino vertex x
      Float_t  nuvtxy_truth;    //neutrino vertex y
      Float_t  nuvtxz_truth;    //neutrino vertex z
      Float_t  nu_dcosx_truth;  //neutrino dcos x
      Float_t  nu_dcosy_truth;  //neutrino dcos y
      Float_t  nu_dcosz_truth;  //neutrino dcos z
      Float_t  lep_mom_truth;   //lepton momentum
      Float_t  lep_dcosx_truth; //lepton dcos x
      Float_t  lep_dcosy_truth; //lepton dcos y
      Float_t  lep_dcosz_truth; //lepton dcos z
      Float_t  t0_truth;        // t0

      // === Storing Geant4 MC Truth Information ===
      int no_primaries;				//<---Number of primary Geant4 particles in the event
      int geant_list_size;				//<---Number of Geant4 particles tracked
      int pdg[kMaxPrimaries];			//<---PDG Code number of this particle
      Float_t Eng[kMaxPrimaries];			//<---Energy of the particle
      Float_t Px[kMaxPrimaries];			//<---Px momentum of the particle
      Float_t Py[kMaxPrimaries];			//<---Py momentum of the particle
      Float_t Pz[kMaxPrimaries];			//<---Pz momentum of the particle
      Float_t StartPointx[kMaxPrimaries];		//<---X position that this Geant4 particle started at
      Float_t StartPointy[kMaxPrimaries];		//<---Y position that this Geant4 particle started at
      Float_t StartPointz[kMaxPrimaries];		//<---Z position that this Geant4 particle started at
      Float_t EndPointx[kMaxPrimaries];		//<---X position that this Geant4 particle ended at
      Float_t EndPointy[kMaxPrimaries];		//<---Y position that this Geant4 particle ended at
      Float_t EndPointz[kMaxPrimaries];		//<---Z position that this Geant4 particle ended at
      Float_t Startdcosx[kMaxPrimaries];            //<---X direction cosine that Geant4 particle started at
      Float_t Startdcosy[kMaxPrimaries];            //<---Y direction cosine that Geant4 particle started at
      Float_t Startdcosz[kMaxPrimaries];            //<---Z direction cosine that Geant4 particle started at
      int NumberDaughters[kMaxPrimaries];		//<---Number of Daughters this particle has
      int TrackId[kMaxPrimaries];			//<---Geant4 TrackID number
      int Mother[kMaxPrimaries];			//<---TrackID of the mother of this particle
      int process_primary[kMaxPrimaries];		//<---Is this particle primary (primary = 1, non-primary = 1)
      std::vector<std::string> G4Process;         //<---The process which created this particle
      std::vector<std::string> G4FinalProcess;    //<---The last process which this particle went under

      //flux information
      Int_t    ptype_flux;        //Parent GEANT code particle ID
      Float_t  pdpx_flux;        //Parent X momentum at decay point (GeV)
      Float_t  pdpy_flux;        //Parent Y momentum at decay point (GeV)
      Float_t  pdpz_flux;        //Parent Z momentum at decay point (GeV)
      Int_t    pntype_flux;      //oscillated neutrino type
      Float_t  vx_flux;          //X position of hadron/muon decay (cm)
      Float_t  vy_flux;          //Y position of hadron/muon decay (cm)
      Float_t  vz_flux;          //Z position of hadron/muon decay (cm)

      //end of ttree variables

      //Module labels to get data products
      std::string fRawDigitModuleLabel;
      std::string fWireModuleLabel;
      std::string fHitsModuleLabel;
      std::string fClusterModuleLabel;
      std::string fTrackModuleLabel;
      std::string fShowerModuleLabel;
      std::string fVertexModuleLabel;
      std::string fGenieGenModuleLabel;
      std::string fPOTModuleLabel;
      std::string fFlashModuleLabel;
      std::string fCalorimetryModuleLabel;

      double fFidVolCut;

      int isinfidvol;
      int isinfidvoltruth;
      float oscpro;

      calo::CalorimetryAlg fCalorimetryAlg;

      bool fMakeAnaTree;
      bool fMakeWeightTree;

      TH1D *mva_nue_osc;
      TH1D *mva_nc;
      TH1D *mva_numu;
      TH1D *mva_nue_beam;
      TH1D *mva_nutau;

      TH1D *mva_numu_nue;
      TH1D *mva_nue_nue;
      TH1D *mva_numu_numu;
      TH1D *mva_nue_numu;
      TH1D *mva_numu_nutau;
      TH1D *mva_nue_nutau;

      TH1D *enu_nc;
      TH1D *enu_numu_nue;
      TH1D *enu_nue_nue;
      TH1D *enu_numu_numu;
      TH1D *enu_nue_numu;
      TH1D *enu_numu_nutau;
      TH1D *enu_nue_nutau;

      TTree *fWeightTree[6]; // each sample plus tree for all
      int nSamples = 5;
      //char name[5][100] = {"#nu_{e}^{osc}","NC","#nu_{#mu} CC","#nu_{e}^{beam}","#nu_{#tau} CC"};
      float events_truth[5];
      float events_reco[5];
      TH1D *enu[5];
      TH1D *enu_osc[5];
      TH1D *hvtxx[5];
      TH1D *hvtxy[5];
      TH1D *hvtxz[5];
      TH1D *htrklen[5];
      TH1D *htrkdedx[5];
      TH1D *hrch[5]; 
      TH1D *hrt[5]; 
      TH1D *hpida[5];
      TH1D *hnshw[5];
      TH1D *heshw[5];
      TH1D *hshwdedx[5];
      TH1D *hdisx[5];
      TH1D *hdisy[5];
      TH1D *hdisz[5];
      TH1D *hfrshower[5];
      TH1D *hnhitspershw[5];
      TH1D *hevtcharge[5];
      TH1D *hfr100w[5];

  }; // class MVAAlg

} // namespace dunemva

#endif // ifndef MVAAlg_H
