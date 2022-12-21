//
// Name: TrackAna_module.cc
//
// Purpose: Module TrackAnaCT.
//
// Configuration parameters.
//
//  TrackModuleLabel:   Track module label.
//  MinMCKE:            Minimum MC particle kinetic energy.
//  MatchColinearity:   Minimum colinearity for mc-track matching.
//  MatchDisp:          Maximum uv displacement for mc-track matching.
//  WMatchDisp:         maximum w displacement of mc-track matching.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/sim.h"

#include "TH2F.h"
#include "TFile.h"
#include "TMatrixD.h"

namespace {

  // Local functions.

  // Calculate distance to boundary.
  //----------------------------------------------------------------------------
  double bdist(const TVector3& pos, unsigned int tpc = 0, unsigned int /*cstat*/ = 0)
  {
    // Get geometry.

    double d3,d4;
    art::ServiceHandle<geo::Geometry> geom;
      
    auto const& tpcg = geom->TPC(geo::TPCID{0, tpc}); // FIXME (KJK): Should probably be TPCID{cstat, tpc}

    if(tpc==2 || tpc==3 || tpc==4 || tpc==5)
      {
	d3 = pos.Y() - 1.25;     // 1.25cm APA 2/3 distance to horizontal.
	//    double d3 = pos.Y() + 85.0;     // Distance to bottom.
	d4 = 113.0 - pos.Y();     // Distance to top.
	//    double d4 = 113.0 - pos.Y();     // Distance to top.
      }
    else  //tpc 0 1  6  7
      {
        d3 = pos.Y() + tpcg.HalfHeight()-15.0;     // Distance to bottom.
	//    double d3 = pos.Y() + 85.0;     // Distance to bottom.
        d4 = tpcg.HalfHeight()+15.0 - pos.Y();     // Distance to top.
	//    double d4 = 113.0 - pos.Y();     // Distance to top.
      }
    
    //    mf::LogVerbatim("output") <<"d3" << d3;
    //    mf::LogVerbatim("output") <<"d4" << d4;

    double d1 = abs(pos.X());          // Distance to right side (wires).
    double d2=2.*tpcg.HalfWidth()- abs(pos.X());
    //    mf::LogVerbatim("output") <<"d1" << d1;
    //    mf::LogVerbatim("output") <<"d2" << d2;
    //    double d2 = 226.539 - pos.X();   // Distance to left side (cathode).
    double d5 = 0.;
    double d6 = 0.;
    
    if(tpc==0 || tpc==1)
      {
	d5 = pos.Z()+1.0;                             // Distance to front.
        d6 = tpcg.Length() -1.0- pos.Z();         // Distance to back.
      }
    else if (tpc==2||tpc==3 || tpc==4 || tpc==5)
      {
	d5 = pos.Z()-51.0;                             // Distance to front.     
        d6 = tpcg.Length() +51.0- pos.Z();         // Distance to back.
      }
    else if (tpc==6 || tpc==7)
      {
	d5 = pos.Z()-103.0;                             // Distance to front.     
        d6 = tpcg.Length() +103.0- pos.Z();         // Distance to back.
	
      }
    if(d6<0){
      //      mf::LogVerbatim("output")<< "z"  <<pos.Z();
      //      mf::LogVerbatim("output")<< "Tpc" <<tpc;
      //      mf::LogVerbatim("output")<< "DetLength" <<tpcg.Length();
      
    }
    //    mf::LogVerbatim("output") <<"d5" << d5;
    //    mf::LogVerbatim("output") <<"d6" << d6;
    double result = std::min({d1, d2, d3, d4, d5, d6});
    //    mf::LogVerbatim("output")<< "bdist" << result;
    //    mf::LogVerbatim("output")<< "Height" << geom->DetHalfHeight(tpc);
    //    mf::LogVerbatim("output")<< "Width" << geom->DetHalfWidth(tpc);
    if(result<0) result=0;
    return result;
   
  }

  // Length of reconstructed track.
  //----------------------------------------------------------------------------
  double length(const recob::Track& track)
  {
    return track.Length();
  }

  // Length of MC particle.
  //----------------------------------------------------------------------------
  double length(detinfo::DetectorPropertiesData const& detProp,
                const simb::MCParticle& part, double dx,
		TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
		unsigned int /*tpc*/ = 0, unsigned int /*cstat*/ = 0)
  {
    // Get services.

    art::ServiceHandle<geo::Geometry> geom;

    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;

    for(int i = 0; i < n; ++i) {
      auto pos = part.Position(i).Vect();

      // Make fiducial cuts.  Require the particle to be within the physical volume of
      // the tpc, and also require the apparent x position to be within the expanded
      // readout frame.

      geo::TPCID tpcid=geom->FindTPCAtPosition(geo::vect::toPoint(pos));
      if (!tpcid.isValid) continue;
      
      pos.SetX(pos.X() + dx);
      
      double ticks=0;
      ticks = detProp.ConvertXToTicks(pos.X(), 0, tpcid.TPC, tpcid.Cryostat);
      if(ticks >= 0. && ticks < detProp.ReadOutWindowSize()) {
	if(first) {
	  start = pos;
	  startmom = part.Momentum(i).Vect();
	}
	else {
	  disp -= pos;
	  result += disp.Mag();
	}
	first = false;
	disp = pos;
	end = pos;
	endmom = part.Momentum(i).Vect();
      }
      
    }
    //    mf::LogVerbatim("output") << " length (MCParticle) " << result;
    return result;
  }
  
  // Fill efficiency histogram assuming binomial errors.

  void effcalc(const TH1* hnum, const TH1* hden, TH1* heff)
  {
    int nbins = hnum->GetNbinsX();
    if (nbins != hden->GetNbinsX())
      throw cet::exception("TrackAnaCT") << "effcalc[" __FILE__ "]: incompatible histograms (I)\n";
    if (nbins != heff->GetNbinsX())
      throw cet::exception("TrackAnaCT") << "effcalc[" __FILE__ "]: incompatible histograms (II)\n";

    // Loop over bins, including underflow and overflow.

    for(int ibin = 0; ibin <= nbins+1; ++ibin) {
      double num = hnum->GetBinContent(ibin);
      double den = hden->GetBinContent(ibin);
      if(den == 0.) {
	heff->SetBinContent(ibin, 0.);
	heff->SetBinError(ibin, 0.);
      }
      else {
	double eff = num / den;
	if(eff < 0.)
	  eff = 0.;
	if(eff > 1.)
	  eff = 1.;
	double err = std::sqrt(eff * (1.-eff) / den);
	heff->SetBinContent(ibin, eff);
	heff->SetBinError(ibin, err);
      }
    }

    heff->SetMinimum(0.);
    heff->SetMaximum(1.05);
    heff->SetMarkerStyle(20);
  }


class flattener : public std::vector<unsigned int> {

public:

 flattener() : std::vector<unsigned int>() {};

 flattener(const std::vector<std::vector<unsigned int> >& input)
 { Convert(input); }

 ~flattener(){}

 void Convert(const std::vector<std::vector<unsigned int> >& input) 
  {
   clear();
   size_t length=0;
   for(auto const& vec : input)
     length += vec.size();
   reserve(length);

   for(auto const& vec : input)
     for(auto const& value : vec)
	push_back(value);

  }
}; // end class flattener

} // end namespace

namespace trkf {

  class TrackAnaCT : public art::EDAnalyzer
  {
  public:

    // Embedded structs.

    // Struct for histograms that depend on reco track only.

    struct RecoHists
    {
      // Constructors.

      RecoHists();
      RecoHists(const std::string& subdir);

      // Pure reco track histograms.

      TH1F* fHstartx;      // Starting x position.
      TH1F* fHstarty;      // Starting y position.
      TH1F* fHstartz;      // Starting z position.
      TH1F* fHstartd;      // Starting distance to boundary.
      TH1F* fHendx;        // Ending x position.
      TH1F* fHendy;        // Ending y position.
      TH1F* fHendz;        // Ending z position.
      TH1F* fHendd;        // Ending distance to boundary.
      TH1F* fHtheta;       // Theta.
      TH1F* fHphi;         // Phi.
      TH1F* fHtheta_xz;    // Theta_xz.
      TH1F* fHtheta_yz;    // Theta_yz.
      TH1F* fHmom;         // Momentum.
      TH1F* fHlen;         // Length.

      // Histograms for the consituent Hits

      TH1F* fHHitChg;       // hit charge
      TH1F* fHHitWidth;     // hit width
      TH1F* fHHitPdg;       // Pdg primarily responsible.
      TH1F* fHHitTrkId;     // TrkId
      TH1F* fModeFrac;       // mode fraction
      TH1F* fNTrkIdTrks;    // # of stitched tracks in which unique TrkId appears
      TH2F* fNTrkIdTrks2;   
      TH2F* fNTrkIdTrks3;   
    };

    // Struct for mc particles and mc-matched tracks.

    struct MCHists
    {
      // Constructors.

      MCHists();
      MCHists(const std::string& subdir);

      // Reco-MC matching.

      TH2F* fHduvcosth;    // 2D mc vs. data matching, duv vs. cos(theta).
      TH1F* fHcosth;       // 1D direction matching, cos(theta).
      TH1F* fHmcu;         // 1D endpoint truth u.
      TH1F* fHmcv;         // 1D endpoint truth v.
      TH1F* fHmcw;         // 1D endpoint truth w.
      TH1F* fHupull;       // 1D endpoint u pull.
      TH1F* fHvpull;       // 1D endpoint v pull.
      TH1F* fHmcdudw;      // Truth du/dw.
      TH1F* fHmcdvdw;      // Truth dv/dw.
      TH1F* fHdudwpull;    // du/dw pull.
      TH1F* fHdvdwpull;    // dv/dw pull.

      // Histograms for matched tracks.

      TH1F* fHstartdx;     // Start dx.
      TH1F* fHstartdy;     // Start dy.
      TH1F* fHstartdz;     // Start dz.
      TH1F* fHenddx;       // End dx.
      TH1F* fHenddy;       // End dy.
      TH1F* fHenddz;       // End dz.
      TH2F* fHlvsl;        // MC vs. reco length.
      TH1F* fHdl;          // Delta(length).
      TH2F* fHpvsp;        // MC vs. reco momentum.
      TH2F* fHpvspc;       // MC vs. reco momentum (contained tracks).
      TH1F* fHdp;          // Momentum difference.
      TH1F* fHdpc;         // Momentum difference (contained tracks).
      TH1F* fHppull;       // Momentum pull.
      TH1F* fHppullc;      // Momentum pull (contained tracks).

      // Pure MC particle histograms (efficiency denominator).

      TH1F* fHmcstartx;    // Starting x position.
      TH1F* fHmcstarty;    // Starting y position.
      TH1F* fHmcstartz;    // Starting z position.
      TH1F* fHmcendx;      // Ending x position.
      TH1F* fHmcendy;      // Ending y position.
      TH1F* fHmcendz;      // Ending z position.
      TH1F* fHmctheta;     // Theta.
      TH1F* fHmcphi;       // Phi.
      TH1F* fHmctheta_xz;  // Theta_xz.
      TH1F* fHmctheta_yz;  // Theta_yz.
      TH1F* fHmcmom;       // Momentum.
      TH1F* fHmclen;       // Length.

      // Histograms for well-reconstructed matched tracks (efficiency numerator).

      TH1F* fHgstartx;     // Starting x position.
      TH1F* fHgstarty;     // Starting y position.
      TH1F* fHgstartz;     // Starting z position.
      TH1F* fHgendx;       // Ending x position.
      TH1F* fHgendy;       // Ending y position.
      TH1F* fHgendz;       // Ending z position.
      TH1F* fHgtheta;      // Theta.
      TH1F* fHgphi;        // Phi.
      TH1F* fHgtheta_xz;   // Theta_xz.
      TH1F* fHgtheta_yz;   // Theta_yz.
      TH1F* fHgmom;        // Momentum.
      TH1F* fHglen;        // Length.

      // Efficiency histograms.

      TH1F* fHestartx;     // Starting x position.
      TH1F* fHestarty;     // Starting y position.
      TH1F* fHestartz;     // Starting z position.
      TH1F* fHeendx;       // Ending x position.
      TH1F* fHeendy;       // Ending y position.
      TH1F* fHeendz;       // Ending z position.
      TH1F* fHetheta;      // Theta.
      TH1F* fHephi;        // Phi.
      TH1F* fHetheta_xz;   // Theta_xz.
      TH1F* fHetheta_yz;   // Theta_yz.
      TH1F* fHemom;        // Momentum.
      TH1F* fHelen;        // Length.


    };

    // Constructors, destructor

    explicit TrackAnaCT(fhicl::ParameterSet const& pset);
    virtual ~TrackAnaCT();

    // Overrides.

    void analyze(const art::Event& evt);
    void anaStitch(detinfo::DetectorClocksData const& clockData,
                   detinfo::DetectorPropertiesData const& detProp,
                   const art::Event& evt);
    void endJob();

  private:

    template <typename T> std::vector<size_t> fsort_indexes(const std::vector<T> &v) ;

    // Fcl Attributes.

    std::string fTrackModuleLabel;
    std::string fSpacepointModuleLabel;
    std::string fStitchModuleLabel;
    std::string fTrkSpptAssocModuleLabel;
    std::string fHitSpptAssocModuleLabel;

    int fDump;                 // Number of events to dump to debug message facility.
    double fMinMCKE;           // Minimum MC particle kinetic energy (GeV).
    double fMinMCLen;          // Minimum MC particle length in tpc (cm).
    double fMatchColinearity;  // Minimum matching colinearity.
    double fMatchDisp;         // Maximum matching displacement.
    double fWMatchDisp;        // Maximum matching displacement in the w direction.
    bool fIgnoreSign;          // Ignore sign of mc particle if true.
    bool fStitchedAnalysis;    // if true, do the whole drill-down from stitched track to assd hits

    // Histograms.

    std::map<int, MCHists> fMCHistMap;       // Indexed by pdg id.
    std::map<int, RecoHists> fRecoHistMap;   // Indexed by pdg id.

    // Statistics.

    int fNumEvent;
  };

  DEFINE_ART_MODULE(TrackAnaCT)

  // RecoHists methods.

  TrackAnaCT::RecoHists::RecoHists() :
    //
    // Purpose: Default constructor.
    //
    fHstartx(0),
    fHstarty(0),
    fHstartz(0),
    fHstartd(0),
    fHendx(0),
    fHendy(0),
    fHendz(0),
    fHendd(0),
    fHtheta(0),
    fHphi(0),
    fHtheta_xz(0),
    fHtheta_yz(0),
    fHmom(0),
    fHlen(0)
    ,fHHitChg(0)
    ,fHHitWidth(0)
    ,fHHitPdg(0)
    ,fHHitTrkId(0)
    ,fModeFrac(0)
    ,fNTrkIdTrks(0)
    ,fNTrkIdTrks2(0)
    ,fNTrkIdTrks3(0)
  {}

  TrackAnaCT::RecoHists::RecoHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.

    *this = RecoHists();

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAnaCT histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.
    auto const& cryostat = geom->Cryostat(geo::CryostatID{0});
    fHstartx = dir.make<TH1F>("xstart", "X Start Position",
                              100, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHstarty = dir.make<TH1F>("ystart", "Y Start Position",
                              100, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHstartz = dir.make<TH1F>("zstart", "Z Start Position",
                              100, 0., cryostat.Length());
    fHstartd = dir.make<TH1F>("dstart", "Start Position Distance to Boundary",
                              100, -10., cryostat.HalfWidth());
    fHendx = dir.make<TH1F>("xend", "X End Position",
                            100, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHendy = dir.make<TH1F>("yend", "Y End Position",
                            100, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHendz = dir.make<TH1F>("zend", "Z End Position",
                            100, 0., cryostat.Length());
    fHendd = dir.make<TH1F>("dend", "End Position Distance to Boundary",
                            100, -10., cryostat.HalfWidth());
    fHtheta = dir.make<TH1F>("theta", "Theta", 100, 0., 3.142);
    fHphi = dir.make<TH1F>("phi", "Phi", 100, -3.142, 3.142);
    fHtheta_xz = dir.make<TH1F>("theta_xz", "Theta_xz", 100, -3.142, 3.142);
    fHtheta_yz = dir.make<TH1F>("theta_yz", "Theta_yz", 100, -3.142, 3.142);
    fHmom = dir.make<TH1F>("mom", "Momentum", 100, 0., 10.);
    fHlen = dir.make<TH1F>("len", "Track Length", 100, 0., 3.0 * cryostat.Length());
    fHHitChg = dir.make<TH1F>("hchg", "Hit Charge (ADC counts)", 100, 0., 4000.);
    fHHitWidth = dir.make<TH1F>("hwid", "Hit Width (ticks)", 40, 0., 20.);
    fHHitPdg = dir.make<TH1F>("hpdg", "Hit Pdg code",5001, -2500.5, +2500.5);
    fHHitTrkId = dir.make<TH1F>("htrkid", "Hit Track ID", 401, -200.5, +200.5);
    fModeFrac = dir.make<TH1F>("hmodefrac", "quasi-Purity: Fraction of component tracks with the Track mode value", 20, 0.01, 1.01);
    fNTrkIdTrks = dir.make<TH1F>("hntrkids", "quasi-Efficiency: Number of stitched tracks in which TrkId appears", 20, 0., +10.0);
    fNTrkIdTrks2 = dir.make<TH2F>("hntrkids2", "Number of stitched tracks in which TrkId appears vs KE [GeV]", 20, 0., +10.0, 20, 0.0, 1.5);
    fNTrkIdTrks3 = dir.make<TH2F>("hntrkids3", "MC Track vs Reco Track, wtd by nhits on Collection Plane", 10, -0.5, 9.5, 10, -0.5, 9.5);
    fNTrkIdTrks3->GetXaxis()->SetTitle("Sorted-by-Descending-CPlane-Hits outer Track Number");
    fNTrkIdTrks3->GetYaxis()->SetTitle("Sorted-by-Descending-True-Length G4Track");
 
  }

  // MCHists methods.

  TrackAnaCT::MCHists::MCHists() :
    //
    // Purpose: Default constructor.
    //
    fHduvcosth(0),
    fHcosth(0),
    fHmcu(0),
    fHmcv(0),
    fHmcw(0),
    fHupull(0),
    fHvpull(0),
    fHmcdudw(0),
    fHmcdvdw(0),
    fHdudwpull(0),
    fHdvdwpull(0),
    fHstartdx(0),
    fHstartdy(0),
    fHstartdz(0),
    fHenddx(0),
    fHenddy(0),
    fHenddz(0),
    fHlvsl(0),
    fHdl(0),
    fHpvsp(0),
    fHpvspc(0),
    fHdp(0),
    fHdpc(0),
    fHppull(0),
    fHppullc(0),
    fHmcstartx(0),
    fHmcstarty(0),
    fHmcstartz(0),
    fHmcendx(0),
    fHmcendy(0),
    fHmcendz(0),
    fHmctheta(0),
    fHmcphi(0),
    fHmctheta_xz(0),
    fHmctheta_yz(0),
    fHmcmom(0),
    fHmclen(0),
    fHgstartx(0),
    fHgstarty(0),
    fHgstartz(0),
    fHgendx(0),
    fHgendy(0),
    fHgendz(0),
    fHgtheta(0),
    fHgphi(0),
    fHgtheta_xz(0),
    fHgtheta_yz(0),
    fHgmom(0),
    fHglen(0),
    fHestartx(0),
    fHestarty(0),
    fHestartz(0),
    fHeendx(0),
    fHeendy(0),
    fHeendz(0),
    fHetheta(0),
    fHephi(0),
    fHetheta_xz(0),
    fHetheta_yz(0),
    fHemom(0),
    fHelen(0)
  {}

  TrackAnaCT::MCHists::MCHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.

    *this = MCHists();

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAnaCT histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.
    auto const& cryostat = geom->Cryostat(geo::CryostatID{0});

    fHduvcosth = dir.make<TH2F>("duvcosth", "Delta(uv) vs. Colinearity", 
				100, 0.95, 1., 100, 0., 1.);
    fHcosth = dir.make<TH1F>("colin", "Colinearity", 100, 0.95, 1.);
    fHmcu = dir.make<TH1F>("mcu", "MC Truth U", 100, -5., 5.);
    fHmcv = dir.make<TH1F>("mcv", "MC Truth V", 100, -5., 5.);
    fHmcw = dir.make<TH1F>("mcw", "MC Truth W", 100, -20., 20.);
    fHupull = dir.make<TH1F>("dupull", "U Pull", 100, -20., 20.);
    fHvpull = dir.make<TH1F>("dvpull", "V Pull", 100, -20., 20.);
    fHmcdudw = dir.make<TH1F>("mcdudw", "MC Truth U Slope", 100, -0.2, 0.2);
    fHmcdvdw = dir.make<TH1F>("mcdvdw", "MV Truth V Slope", 100, -0.2, 0.2);
    fHdudwpull = dir.make<TH1F>("dudwpull", "U Slope Pull", 100, -10., 10.);
    fHdvdwpull = dir.make<TH1F>("dvdwpull", "V Slope Pull", 100, -10., 10.);
    fHstartdx = dir.make<TH1F>("startdx", "Start Delta x", 100, -10., 10.);
    fHstartdy = dir.make<TH1F>("startdy", "Start Delta y", 100, -10., 10.);
    fHstartdz = dir.make<TH1F>("startdz", "Start Delta z", 100, -10., 10.);
    fHenddx = dir.make<TH1F>("enddx", "End Delta x", 100, -10., 10.);
    fHenddy = dir.make<TH1F>("enddy", "End Delta y", 100, -10., 10.);
    fHenddz = dir.make<TH1F>("enddz", "End Delta z", 100, -10., 10.);
    fHlvsl = dir.make<TH2F>("lvsl", "Reco Length vs. MC Truth Length",
                            100, 0., 1.1 * cryostat.Length(), 100, 0., 1.1 * cryostat.Length());
    fHdl = dir.make<TH1F>("dl", "Track Length Minus MC Particle Length", 100, -50., 50.);
    fHpvsp = dir.make<TH2F>("pvsp", "Reco Momentum vs. MC Truth Momentum",
			    100, 0., 5., 100, 0., 5.);
    fHpvspc = dir.make<TH2F>("pvspc", "Reco Momentum vs. MC Truth Momentum (Contained Tracks)",
			     100, 0., 5., 100, 0., 5.);
    fHdp = dir.make<TH1F>("dp", "Reco-MC Momentum Difference", 100, -5., 5.);
    fHdpc = dir.make<TH1F>("dpc", "Reco-MC Momentum Difference (Contained Tracks)",
			   100, -5., 5.);
    fHppull = dir.make<TH1F>("ppull", "Momentum Pull", 100, -10., 10.);
    fHppullc = dir.make<TH1F>("ppullc", "Momentum Pull (Contained Tracks)", 100, -10., 10.);

    fHmcstartx = dir.make<TH1F>("mcxstart", "MC X Start Position",
                                10, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHmcstarty = dir.make<TH1F>("mcystart", "MC Y Start Position",
                                10, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHmcstartz = dir.make<TH1F>("mczstart", "MC Z Start Position",
                                10, 0., cryostat.Length());
    fHmcendx = dir.make<TH1F>("mcxend", "MC X End Position",
                              10, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHmcendy = dir.make<TH1F>("mcyend", "MC Y End Position",
                              10, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHmcendz = dir.make<TH1F>("mczend", "MC Z End Position",
                              10, 0., cryostat.Length());
    fHmctheta = dir.make<TH1F>("mctheta", "MC Theta", 20, 0., 3.142);
    fHmcphi = dir.make<TH1F>("mcphi", "MC Phi", 10, -3.142, 3.142);
    fHmctheta_xz = dir.make<TH1F>("mctheta_xz", "MC Theta_xz", 40, -3.142, 3.142);
    fHmctheta_yz = dir.make<TH1F>("mctheta_yz", "MC Theta_yz", 40, -3.142, 3.142);
    fHmcmom = dir.make<TH1F>("mcmom", "MC Momentum", 10, 0., 10.);
    fHmclen = dir.make<TH1F>("mclen", "MC Particle Length", 10, 0., 1.1 * cryostat.Length());

    fHgstartx = dir.make<TH1F>("gxstart", "Good X Start Position",
                               10, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHgstarty = dir.make<TH1F>("gystart", "Good Y Start Position",
                               10, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHgstartz = dir.make<TH1F>("gzstart", "Good Z Start Position",
                               10, 0., cryostat.Length());
    fHgendx = dir.make<TH1F>("gxend", "Good X End Position",
                             10, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHgendy = dir.make<TH1F>("gyend", "Good Y End Position",
                             10, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHgendz = dir.make<TH1F>("gzend", "Good Z End Position",
                             10, 0., cryostat.Length());
    fHgtheta = dir.make<TH1F>("gtheta", "Good Theta", 20, 0., 3.142);
    fHgphi = dir.make<TH1F>("gphi", "Good Phi", 10, -3.142, 3.142);
    fHgtheta_xz = dir.make<TH1F>("gtheta_xz", "Good Theta_xz", 40, -3.142, 3.142);
    fHgtheta_yz = dir.make<TH1F>("gtheta_yz", "Good Theta_yz", 40, -3.142, 3.142);
    fHgmom = dir.make<TH1F>("gmom", "Good Momentum", 10, 0., 10.);
    fHglen = dir.make<TH1F>("glen", "Good Particle Length", 10, 0., 1.1 * cryostat.Length());

    fHestartx = dir.make<TH1F>("exstart", "Efficiency vs. X Start Position",
                               10, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHestarty = dir.make<TH1F>("eystart", "Efficiency vs. Y Start Position",
                               10, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHestartz = dir.make<TH1F>("ezstart", "Efficiency vs. Z Start Position",
                               10, 0., cryostat.Length());
    fHeendx = dir.make<TH1F>("exend", "Efficiency vs. X End Position",
                             10, -2.*cryostat.HalfWidth(), 4.*cryostat.HalfWidth());
    fHeendy = dir.make<TH1F>("eyend", "Efficiency vs. Y End Position",
                             10, -cryostat.HalfHeight(), cryostat.HalfHeight());
    fHeendz = dir.make<TH1F>("ezend", "Efficiency vs. Z End Position",
                             10, 0., cryostat.Length());
    fHetheta = dir.make<TH1F>("etheta", "Efficiency vs. Theta", 20, 0., 3.142);
    fHephi = dir.make<TH1F>("ephi", "Efficiency vs. Phi", 10, -3.142, 3.142);
    fHetheta_xz = dir.make<TH1F>("etheta_xz", "Efficiency vs. Theta_xz", 40, -3.142, 3.142);
    fHetheta_yz = dir.make<TH1F>("etheta_yz", "Efficiency vs. Theta_yz", 40, -3.142, 3.142);
    fHemom = dir.make<TH1F>("emom", "Efficiency vs. Momentum", 10, 0., 10.);
    fHelen = dir.make<TH1F>("elen", "Efficiency vs. Particle Length",
                            10, 0., 1.1 * cryostat.Length());
  }

  TrackAnaCT::TrackAnaCT(const fhicl::ParameterSet& pset)
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    : EDAnalyzer(pset)
    , fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel"))
    , fSpacepointModuleLabel(pset.get<std::string>("SpacepointModuleLabel"))
    , fStitchModuleLabel(pset.get<std::string>("StitchModuleLabel"))
    , fTrkSpptAssocModuleLabel(pset.get<std::string>("TrkSpptAssocModuleLabel"))
    , fHitSpptAssocModuleLabel(pset.get<std::string>("HitSpptAssocModuleLabel"))
    , fDump(pset.get<int>("Dump"))
    , fMinMCKE(pset.get<double>("MinMCKE"))
    , fMinMCLen(pset.get<double>("MinMCLen"))
    , fMatchColinearity(pset.get<double>("MatchColinearity"))
    , fMatchDisp(pset.get<double>("MatchDisp"))
    , fWMatchDisp(pset.get<double>("WMatchDisp"))
    , fIgnoreSign(pset.get<bool>("IgnoreSign"))
    , fStitchedAnalysis(pset.get<bool>("StitchedAnalysis",false))
    , fNumEvent(0)
  {

    ///\todo Move this module to DUNE code and remove it from larreco

    art::ServiceHandle<geo::Geometry> geom;
    if(geom->DetectorName().find("dune") == std::string::npos)
      throw cet::exception("TrackAnaCT") << "TrackAnaCT should only be used with DUNE "
					 << "geometries, the name for this detector, "
					 << geom->DetectorName() << ", does not contain "
					 << "dune, bail.";

    // Report.

    mf::LogInfo("TrackAnaCT") 
      << "TrackAnaCT configured with the following parameters:\n"
      << "  TrackModuleLabel = " << fTrackModuleLabel << "\n"
      << "  StitchModuleLabel = " << fStitchModuleLabel << "\n"
      << "  TrkSpptAssocModuleLabel = " << fTrkSpptAssocModuleLabel << "\n"
      << "  HitSpptAssocModuleLabel = " << fHitSpptAssocModuleLabel << "\n"
      << "  Dump = " << fDump << "\n"
      << "  MinMCKE = " << fMinMCKE << "\n"
      << "  MinMCLen = " << fMinMCLen;
  }

  TrackAnaCT::~TrackAnaCT()
  //
  // Purpose: Destructor.
  //
  {}

  void TrackAnaCT::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    art::ServiceHandle<geo::Geometry> geom;


    ++fNumEvent;

    // Optional dump stream.

    std::unique_ptr<mf::LogInfo> pdump;
    if(fDump > 0) {
      --fDump;
      pdump = std::unique_ptr<mf::LogInfo>(new mf::LogInfo("TrackAnaCT"));
    }

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();

    // Get mc particles.

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    std::vector<const simb::MCParticle*> plist2;
    if(mc) {

//      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      sim::ParticleList const& plist = pi_serv->ParticleList();
      plist2.reserve(plist.size());
  

      if(pdump) {
	*pdump << "MC Particles\n"
	       << "       Id   pdg           x         y         z          dx        dy        dz           p\n"
	       << "-------------------------------------------------------------------------------------------\n";
      }

      // Loop over mc particles, and fill histograms that depend only
      // on mc particles.  Also, fill a secondary list of mc particles
      // that pass various selection criteria.

      for(sim::ParticleList::const_iterator ipart = plist.begin();
	  ipart != plist.end(); ++ipart) {
	const simb::MCParticle* part = (*ipart).second;
	if (!part)
	  throw cet::exception("SeedAna") << "no particle!\n";
	int pdg = part->PdgCode();
	if(fIgnoreSign)
	  pdg = std::abs(pdg);

	// Ignore everything except stable charged nonshowering particles.

	int apdg = std::abs(pdg);
	if(apdg == 13 ||     // Muon
	   apdg == 211 ||    // Charged pion
	   apdg == 321 ||    // Charged kaon
	   apdg == 2212) {   // (Anti)proton

	  // Apply minimum energy cut.

	  if(part->E() >= 0.001*part->Mass() + fMinMCKE) {

	    // Calculate the x offset due to nonzero mc particle time.

	    double mctime = part->T();                                // nsec
	    double mcdx = mctime * 1.e-3 * detProp.DriftVelocity();  // cm

	    // Calculate the length of this mc particle inside the fiducial volume.

	    TVector3 mcstart;
	    TVector3 mcend;
	    TVector3 mcstartmom;
	    TVector3 mcendmom;
            double plen = length(detProp, *part, mcdx, mcstart, mcend, mcstartmom, mcendmom);

	    // Apply minimum fiducial length cut.  Always reject particles that have
	    // zero length in the tpc regardless of the configured cut.

	    if(plen > 0. && plen > fMinMCLen) {

	      // This is a good mc particle (capable of making a track).

	      plist2.push_back(part);

	      // Dump MC particle information here.

	      if(pdump) {
		double pstart = mcstartmom.Mag();
		double pend = mcendmom.Mag();
		*pdump << "\nOffset"
		       << std::setw(3) << part->TrackId()
		       << std::setw(6) << part->PdgCode()
		       << "  " 
		       << std::fixed << std::setprecision(2) 
		       << std::setw(10) << mcdx
		       << "\nStart " 
		       << std::setw(3) << part->TrackId()
		       << std::setw(6) << part->PdgCode()
		       << "  " 
		       << std::fixed << std::setprecision(2) 
		       << std::setw(10) << mcstart[0]
		       << std::setw(10) << mcstart[1]
		       << std::setw(10) << mcstart[2];
		if(pstart > 0.) {
		  *pdump << "  "
			 << std::fixed << std::setprecision(3) 
			 << std::setw(10) << mcstartmom[0]/pstart
			 << std::setw(10) << mcstartmom[1]/pstart
			 << std::setw(10) << mcstartmom[2]/pstart;
		}
		else
		  *pdump << std::setw(32) << " ";
		*pdump << std::setw(12) << std::fixed << std::setprecision(3) << pstart;
		*pdump << "\nEnd   " 
		       << std::setw(3) << part->TrackId()
		       << std::setw(6) << part->PdgCode()
		       << "  " 
		       << std::fixed << std::setprecision(2)
		       << std::setw(10) << mcend[0]
		       << std::setw(10) << mcend[1]
		       << std::setw(10) << mcend[2];
		if(pend > 0.01) {
		  *pdump << "  " 
			 << std::fixed << std::setprecision(3) 
			 << std::setw(10) << mcendmom[0]/pend
			 << std::setw(10) << mcendmom[1]/pend
			 << std::setw(10) << mcendmom[2]/pend;
		}
		else
		  *pdump << std::setw(32)<< " ";
		*pdump << std::setw(12) << std::fixed << std::setprecision(3) << pend << "\n";
	      }

	      // Fill histograms.

	      if(fMCHistMap.count(pdg) == 0) {
		std::ostringstream ostr;
		ostr << "MC" << (fIgnoreSign ? "All" : (pdg > 0 ? "Pos" : "Neg")) << std::abs(pdg);
		fMCHistMap[pdg] = MCHists(ostr.str());
	      }
	      const MCHists& mchists = fMCHistMap[pdg];

	      double mctheta_xz = std::atan2(mcstartmom.X(), mcstartmom.Z());
	      double mctheta_yz = std::atan2(mcstartmom.Y(), mcstartmom.Z());

	      mchists.fHmcstartx->Fill(mcstart.X());
	      mchists.fHmcstarty->Fill(mcstart.Y());
	      mchists.fHmcstartz->Fill(mcstart.Z());
	      mchists.fHmcendx->Fill(mcend.X());
	      mchists.fHmcendy->Fill(mcend.Y());
	      mchists.fHmcendz->Fill(mcend.Z());
	      mchists.fHmctheta->Fill(mcstartmom.Theta());
	      mchists.fHmcphi->Fill(mcstartmom.Phi());
	      mchists.fHmctheta_xz->Fill(mctheta_xz);
	      mchists.fHmctheta_yz->Fill(mctheta_yz);
	      mchists.fHmcmom->Fill(mcstartmom.Mag());
	      mchists.fHmclen->Fill(plen);
	    }
	  }
	}
      }
    } //mc

	
    // Get tracks and spacepoints and hits

    auto trackh = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
    auto  trackvh = evt.getHandle< std::vector< art::PtrVector < recob::Track > > >(fStitchModuleLabel);

    // This new top part of TrackAnaCT between two long lines of ************s
    // is particular to analyzing Stitched Tracks.
    // *******************************************************************//

    if (trackvh.isValid() && fStitchedAnalysis) 
      {
	mf::LogDebug("TrackAnaCT") 
	  << "TrackAnaCT read "  << trackvh->size()
	  << "  vectors of Stitched PtrVectorsof tracks.";
        anaStitch(clockData, detProp, evt);
      }

    if(trackh.isValid()) {

      if(pdump) {
	*pdump << "\nReconstructed Tracks\n"
	       << "       Id  MCid           x         y         z          dx        dy        dz           p\n"
	       << "-------------------------------------------------------------------------------------------\n";
      }

      // Loop over tracks.
     

      int ntrack = trackh->size();
      for(int i = 0; i < ntrack; ++i) {
	art::Ptr<recob::Track> ptrack(trackh, i);
	const recob::Track& track = *ptrack;
	art::FindManyP<recob::Hit> fh(trackh, evt, fTrkSpptAssocModuleLabel);

	////
	///              figuring out which TPC
	///
	///
	//
	//	auto pcoll { ptrack };
	//art::FindManyP<recob::SpacePoint> fs( pcoll, evt, fTrkSpptAssocModuleLabel);
	//	auto sppt = fs.at(0);//.at(is);
	//	art::FindManyP<recob::Hit> fh( sppt, evt, fHitSpptAssocModuleLabel);
	auto hit = fh.at(0).at(0);
	geo::WireID tmpWireid=hit->WireID();
	int hit_tpc=tmpWireid.TPC;
	///
	///
	//
	//
	//
	//
	//
	//
	
	
 

	/*	
	std::vector<art::Ptr<recob::Hit> > hitlist;
        auto  hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel);
	if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);
        */

	// Calculate the x offset due to nonzero reconstructed time.

	//double recotime = track.Time() * sampling_rate(clockData);       // nsec
	double recotime = 0.;
	double trackdx = recotime * 1.e-3 * detProp.DriftVelocity();  // cm

	// Fill histograms involving reco tracks only.
	
	int ntraj = track.NumberTrajectoryPoints();
	if(ntraj > 0) {
	  TVector3 pos = track.Vertex<TVector3>();
	  TVector3 dir = track.VertexDirection<TVector3>();
	  TVector3 end = track.End<TVector3>();
	  pos[0] += trackdx;
	  end[0] += trackdx;
	  
	  double dpos = bdist(pos,hit_tpc);
	  double dend = bdist(end,hit_tpc);
	  double tlen = length(track);
	  double theta_xz = std::atan2(dir.X(), dir.Z());
	  double theta_yz = std::atan2(dir.Y(), dir.Z());
	  
	  if(fRecoHistMap.count(0) == 0)
	    fRecoHistMap[0] = RecoHists("Reco");
	  const RecoHists& rhists = fRecoHistMap[0];
	  
	  rhists.fHstartx->Fill(pos.X());
	  rhists.fHstarty->Fill(pos.Y());
	  rhists.fHstartz->Fill(pos.Z());
	  rhists.fHstartd->Fill(dpos);
	  rhists.fHendx->Fill(end.X());
	  rhists.fHendy->Fill(end.Y());
	  rhists.fHendz->Fill(end.Z());
	  rhists.fHendd->Fill(dend);
	  rhists.fHtheta->Fill(dir.Theta());
	  rhists.fHphi->Fill(dir.Phi());
	  rhists.fHtheta_xz->Fill(theta_xz);
	  rhists.fHtheta_yz->Fill(theta_yz);
	  
	  double mom = 0.;
	  if(track.NumberTrajectoryPoints() > 0)
	    mom = track.VertexMomentum();
	  rhists.fHmom->Fill(mom);
	  rhists.fHlen->Fill(tlen);

	  // Id of matching mc particle.

	  int mcid = -1;

	  // Loop over direction.  

	  for(int swap=0; swap<2; ++swap) {

	    // Analyze reversed tracks only if start momentum = end momentum.

	    if(swap != 0 && track.NumberTrajectoryPoints() > 0 &&
	       std::abs(track.VertexMomentum() - track.EndMomentum()) > 1.e-3)
	      continue;

	    // Calculate the global-to-local rotation matrix.

	    int start_point = (swap == 0 ? 0 : ntraj-1);
	    TMatrixD rot = track.GlobalToLocalRotationAtPoint<TMatrixD>(start_point);

	    // Update track data for reversed case.

	    if(swap != 0) {
	      rot(1, 0) = -rot(1, 0);
	      rot(2, 0) = -rot(2, 0);
	      rot(1, 1) = -rot(1, 1);
	      rot(2, 1) = -rot(2, 1);
	      rot(1, 2) = -rot(1, 2);
	      rot(2, 2) = -rot(2, 2);

	      pos = track.End<TVector3>();
	      dir = -track.EndDirection<TVector3>();
	      end = track.Vertex<TVector3>();
	      pos[0] += trackdx;
	      end[0] += trackdx;
	  
	      dpos = bdist(pos,hit_tpc);
	      dend = bdist(end,hit_tpc);
	      theta_xz = std::atan2(dir.X(), dir.Z());
	      theta_yz = std::atan2(dir.Y(), dir.Z());

	      if(track.NumberTrajectoryPoints() > 0)
		mom = track.EndMomentum();
	    }
	  
	    // Get covariance matrix.

	    //	    const TMatrixD& cov = (swap == 0 ? track.VertexCovariance() : track.EndCovariance());
	    
	    // Loop over track-like mc particles.

	    for(auto ipart = plist2.begin(); ipart != plist2.end(); ++ipart) {
	      const simb::MCParticle* part = *ipart;
	      if (!part)
	        throw cet::exception("SeedAna") << "no particle! [II]\n";
	      int pdg = part->PdgCode();
	      if(fIgnoreSign) pdg = std::abs(pdg);
	      auto iMCHistMap = fMCHistMap.find(pdg);
	      if (iMCHistMap == fMCHistMap.end())
	        throw cet::exception("SeedAna") << "no particle with ID=" << pdg << "\n";
	      const MCHists& mchists = iMCHistMap->second;

	      // Calculate the x offset due to nonzero mc particle time.

	      double mctime = part->T();                                 // nsec
	      double mcdx = mctime * 1.e-3 * detProp.DriftVelocity();   // cm

	      // Calculate the points where this mc particle enters and leaves the
	      // fiducial volume, and the length in the fiducial volume.

	      TVector3 mcstart;
	      TVector3 mcend;
	      TVector3 mcstartmom;
	      TVector3 mcendmom;
              double plen = length(detProp, *part, mcdx, mcstart, mcend, mcstartmom, mcendmom);

	      // Get the displacement of this mc particle in the global coordinate system.

	      TVector3 mcpos = mcstart - pos;
	    
	      // Rotate the momentum and position to the
	      // track-local coordinate system.

	      TVector3 mcmoml = rot * mcstartmom;
	      TVector3 mcposl = rot * mcpos;

	      double colinearity = mcmoml.Z() / mcmoml.Mag();
	    
	      double u = mcposl.X();
	      double v = mcposl.Y();
	      double w = mcposl.Z();
	    
	      double pu = mcmoml.X();
	      double pv = mcmoml.Y();
	      double pw = mcmoml.Z();

	      double dudw = pu / pw;
	      double dvdw = pv / pw;

	      double u0 = u - w * dudw;
	      double v0 = v - w * dvdw;
	      double uv0 = std::sqrt(u0*u0 + v0*v0);

	      mchists.fHduvcosth->Fill(colinearity, uv0);
	      if(std::abs(uv0) < fMatchDisp) {

		// Fill slope matching histograms.

		mchists.fHmcdudw->Fill(dudw);
		mchists.fHmcdvdw->Fill(dvdw);
		//		mchists.fHdudwpull->Fill(dudw / std::sqrt(cov(2,2)));
		//		mchists.fHdvdwpull->Fill(dvdw / std::sqrt(cov(3,3)));
	      }
	      mchists.fHcosth->Fill(colinearity);
	      if(colinearity > fMatchColinearity) {

		// Fill displacement matching histograms.

		mchists.fHmcu->Fill(u0);
		mchists.fHmcv->Fill(v0);
		mchists.fHmcw->Fill(w);
		//		mchists.fHupull->Fill(u0 / std::sqrt(cov(0,0)));
		//		mchists.fHvpull->Fill(v0 / std::sqrt(cov(1,1)));
	      
		if(std::abs(uv0) < fMatchDisp) {

		  // Fill matching histograms.

		  double mctheta_xz = std::atan2(mcstartmom.X(), mcstartmom.Z());
		  double mctheta_yz = std::atan2(mcstartmom.Y(), mcstartmom.Z());

		  mchists.fHstartdx->Fill(pos.X() - mcstart.X());
		  mchists.fHstartdy->Fill(pos.Y() - mcstart.Y());
		  mchists.fHstartdz->Fill(pos.Z() - mcstart.Z());
		  mchists.fHenddx->Fill(end.X() - mcend.X());
		  mchists.fHenddy->Fill(end.Y() - mcend.Y());
		  mchists.fHenddz->Fill(end.Z() - mcend.Z());
		  mchists.fHlvsl->Fill(plen, tlen);
		  mchists.fHdl->Fill(tlen - plen);
		  mchists.fHpvsp->Fill(mcstartmom.Mag(), mom);
		  double dp = mom - mcstartmom.Mag();
		  mchists.fHdp->Fill(dp);
		  //		  mchists.fHppull->Fill(dp / std::sqrt(cov(4,4)));
		  if(std::abs(dpos) >= 5. && std::abs(dend) >= 5.) {
		    mchists.fHpvspc->Fill(mcstartmom.Mag(), mom);
		    mchists.fHdpc->Fill(dp);
		    //		    mchists.fHppullc->Fill(dp / std::sqrt(cov(4,4)));
		  }

		  // Count this track as well-reconstructed if it is matched to an
		  // mc particle (which it is if get here), and if in addition the
		  // starting position (w) matches and the reconstructed track length
		  // is more than 0.5 of the mc particle trajectory length.

		  bool good = std::abs(w) <= fWMatchDisp &&
		    tlen > 0.5 * plen;
		  if(good) {
		    mcid = part->TrackId();
		    mchists.fHgstartx->Fill(mcstart.X());
		    mchists.fHgstarty->Fill(mcstart.Y());
		    mchists.fHgstartz->Fill(mcstart.Z());
		    mchists.fHgendx->Fill(mcend.X());
		    mchists.fHgendy->Fill(mcend.Y());
		    mchists.fHgendz->Fill(mcend.Z());
		    mchists.fHgtheta->Fill(mcstartmom.Theta());
		    mchists.fHgphi->Fill(mcstartmom.Phi());
		    mchists.fHgtheta_xz->Fill(mctheta_xz);
		    mchists.fHgtheta_yz->Fill(mctheta_yz);
		    mchists.fHgmom->Fill(mcstartmom.Mag());
		    mchists.fHglen->Fill(plen);
		  }
		}
	      }
	    }
	  }
	  
	  // Dump track information here.

	  if(pdump) {
	    TVector3 pos = track.Vertex<TVector3>();
	    TVector3 dir = track.VertexDirection<TVector3>();
	    TVector3 end = track.End<TVector3>();
	    pos[0] += trackdx;
	    end[0] += trackdx;
	    TVector3 enddir = track.EndDirection<TVector3>();
	    double pstart = track.VertexMomentum();
	    double pend = track.EndMomentum();
	    *pdump << "\nOffset"
		   << std::setw(3) << track.ID()
		   << std::setw(6) << mcid
		   << "  "
		   << std::fixed << std::setprecision(2) 
		   << std::setw(10) << trackdx
		   << "\nStart " 
		   << std::setw(3) << track.ID()
		   << std::setw(6) << mcid
		   << "  "
		   << std::fixed << std::setprecision(2) 
		   << std::setw(10) << pos[0]
		   << std::setw(10) << pos[1]
		   << std::setw(10) << pos[2];
	    if(pstart > 0.) {
	      *pdump << "  "
		     << std::fixed << std::setprecision(3) 
		     << std::setw(10) << dir[0]
		     << std::setw(10) << dir[1]
		     << std::setw(10) << dir[2];
	    }
	    else
	      *pdump << std::setw(32) << " ";
	    *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pstart;
	    *pdump << "\nEnd   " 
		   << std::setw(3) << track.ID()
		   << std::setw(6) << mcid
		   << "  "
		   << std::fixed << std::setprecision(2)
		   << std::setw(10) << end[0]
		   << std::setw(10) << end[1]
		   << std::setw(10) << end[2];
	    if(pend > 0.01) {
	      *pdump << "  " 
		     << std::fixed << std::setprecision(3) 
		     << std::setw(10) << enddir[0]
		     << std::setw(10) << enddir[1]
		     << std::setw(10) << enddir[2];
	    }
	    else 
	      *pdump << std::setw(32)<< " ";
	    *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pend << "\n";
	  }
	}
      }
    }   // i

  }

  void TrackAnaCT::anaStitch(detinfo::DetectorClocksData const& clockData,
                             detinfo::DetectorPropertiesData const& detProp,
                             const art::Event& evt)
  {

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<geo::Geometry> geom;

    std::map<int, std::map<int, art::PtrVector<recob::Hit>> > hitmap; // trkID, otrk, hitvec
    std::map<int, int > KEmap; // length traveled in det [cm]?, trkID want to sort by KE
    bool mc = !evt.isRealData();

    auto trackh = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
    auto sppth = evt.getHandle< std::vector< recob::SpacePoint> >(fSpacepointModuleLabel);
    auto trackvh = evt.getHandle< std::vector< art::PtrVector < recob::Track > > >(fStitchModuleLabel);

    int ntv(trackvh->size());
    
    std::vector < art::PtrVector<recob::Track> >::const_iterator cti = trackvh->begin();
    /// art::FindManyP<recob::Hit> fh(sppth, evt, fHitSpptAssocModuleLabel);
    
    if(trackh.isValid()) {
      art::FindManyP<recob::SpacePoint> fswhole(trackh, evt, fTrkSpptAssocModuleLabel);
      int nsppts_assnwhole = fswhole.size();
      std::cout << "TrackAnaCT: Number of clumps of Spacepoints from Assn for all Tracks: " << nsppts_assnwhole << std::endl;
    }
    
    if(fRecoHistMap.count(0) == 0)
      {
	fRecoHistMap[0] = RecoHists("Reco");
	std::cout << "\n" << "\t\t  TrkAna: Fresh fRecoHistMap[0] ******* \n" << std::endl;
      }
    const RecoHists& rhistsStitched = fRecoHistMap[0];
    
    std::vector < std::vector <unsigned int> >  NtrkIdsAll; 
    std::vector < double > ntvsorted;
    hitmap.clear();
    KEmap.clear();
    

    for (int o = 0; o < ntv; ++o) // o for outer
      {

	const art::PtrVector<recob::Track> pvtrack(*(cti++));
	auto it = pvtrack.begin();
	int ntrack = pvtrack.size();
	//	if (ntrack>1) 	std::cout << "\t\t  TrkAna: New Stitched Track ******* " << std::endl;
	std::vector< std::vector <unsigned int> > NtrkId_Hit; // hit IDs in inner tracks
	std::vector<unsigned int> vecMode;

	for(int i = 0; i < ntrack; ++i) {

	  const art::Ptr<recob::Track> ptrack(*(it++));
	  //	  const recob::Track& track = *ptrack;
	  auto pcoll={ ptrack };
	  art::FindManyP<recob::SpacePoint> fs( pcoll, evt, fTrkSpptAssocModuleLabel);
	  // From gdb> ptype fs, the vector of Ptr<SpacePoint>s it appears is grabbed after fs.at(0)
	  bool assns(true);
	  try {
	    // Get Spacepoints from this Track, get Hits from those Spacepoints.
	    //	    int nsppts = ptrack->NumberTrajectoryPoints();
	    
	    int nsppts_assn = fs.at(0).size();  
	    //	    if (ntrack>1) std::cout << "\t\tTrackAnaCT: Number of Spacepoints from Track.NumTrajPts(): " << nsppts << std::endl;
	    //	    if (ntrack>1)  std::cout << "\t\tTrackAnaCT: Number of Spacepoints from Assns for this Track: " << nsppts_assn << std::endl;
	    //assert (nsppts_assn == nsppts);
	    auto sppt = fs.at(0);//.at(is);
	    art::FindManyP<recob::Hit> fh( sppt, evt, fHitSpptAssocModuleLabel);
	    // Importantly, loop on all sppts, though they don't all contribute to the track.
	    // As opposed to looping on the trajectory pts, which is a lower number. 
	    // Also, important, in job in whch this runs I set TrackKal3DSPS parameter MaxPass=1, 
	    // cuz I don't want merely the sparse set of sppts as follows from the uncontained 
	    // MS-measurement in 2nd pass.
	    std::vector <unsigned int> vecNtrkIds;
	    for(int is = 0; is < nsppts_assn; ++is) {
	      int nhits = fh.at(is).size(); // should be 2 or 3: number of planes.
	      for(int ih = 0; ih < nhits; ++ih) {
		auto hit = fh.at(is).at(ih); // Our vector is after the .at(is) this time.
		if (hit->SignalType()!=geo::kCollection) continue;
		rhistsStitched.fHHitChg->Fill(hit->Integral());
		rhistsStitched.fHHitWidth->Fill(hit->RMS() * 2.);
		if (mc)
		  {
		    std::vector<sim::TrackIDE> tids = bt_serv->HitToTrackIDEs(clockData, hit);
		    // more here.
		    // Loop over track ids.
		    bool justOne(true); // Only take first trk that contributed to this hit
		    //	  std::cout  << "\t\t  TrkAna: TrkId  tids.size() ******* " << tids.size()  <<std::endl;
		    for(std::vector<sim::TrackIDE>::const_iterator itid = tids.begin();itid != tids.end(); ++itid) {
		      int trackID = std::abs(itid->trackID);
		      hitmap[trackID][o].push_back(hit);

		      if (justOne) { vecNtrkIds.push_back(trackID); justOne=false; }
		      // Add hit to PtrVector corresponding to this track id.
		      rhistsStitched.fHHitTrkId->Fill(trackID); 
		      const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(trackID);
		      rhistsStitched.fHHitPdg->Fill(part->PdgCode()); 
		      // This really needs to be indexed as KE deposited in volTPC, not just KE. EC, 24-July-2014.

		      TVector3 mcstart;
		      TVector3 mcend;
		      TVector3 mcstartmom;
		      TVector3 mcendmom;
		      double mctime = part->T();                                 // nsec
		      double mcdx = mctime * 1.e-3 * detProp.DriftVelocity();   // cm

                      double plen = length(detProp, *part, mcdx, mcstart, mcend, mcstartmom, mcendmom);

		      KEmap[(int)(1e6*plen)] = trackID; // multiple assignment but always the same, so fine.
		      //		      std::cout  << "\t\t  TrkAna: TrkId  trackID, KE [MeV] ******* " << trackID << ", " << (int)(1e3*(part->E()-part->Mass()))  <<std::endl;
		    }

		  } // mc
	      } //  hits

	    } //    spacepoints

	    if (mc)
	      {
		NtrkId_Hit.push_back(vecNtrkIds);	
		// Find the trkID mode for this i^th track
		unsigned int ii(1);
		int max(-12), n(1), ind(0);
		std::sort(vecNtrkIds.begin(),vecNtrkIds.end());
		std::vector<unsigned int> strkIds(vecNtrkIds);
		while ( ii < vecNtrkIds.size() )
		  { 
		    if (strkIds.at(ii) != strkIds.at(ii-1)) 
		      {
			n=1;
		      }
		    else
		      {
			n++; 
		      }
		    if (n>max) { max = n; ind = ii;}
		    ii++;
		  }
		// std::cout  << "\t\t  TrkAna: TrkId  ind for this track is ******* " << ind  <<std::endl;
		unsigned int mode(sim::NoParticleId);
		if (strkIds.begin()!=strkIds.end()) 
		  mode = strkIds.at(ind);
		vecMode.push_back(mode);

		//	if (ntrack>1)	std::cout  << "\t\t  TrkAna: TrkId mode for this component track is ******* " << mode <<std::endl;
		if (strkIds.size()!=0)
		  rhistsStitched.fModeFrac->Fill((double)max/(double)strkIds.size());
		else
		  rhistsStitched.fModeFrac->Fill(-1.0);
		} // mc

	  } // end try 
	  catch (cet::exception& x)  {
	    assns = false;
	  }
	  if (!assns) throw cet::exception("TrackAnaCT") << "Bad Associations. \n";

	} // i

	if (mc)
	  {
	    // one vector per o trk, for all modes of stitched i trks
	    NtrkIdsAll.push_back(vecMode); 

	    std::unique(NtrkIdsAll.back().begin(),NtrkIdsAll.back().end());
	    double sum(0.0);
	    for (auto const val :  NtrkIdsAll.back())
	      {
		//		rhistsStitched.fNTrkIdTrks3->Fill(o,val%100,hitmap[val][o].size());
		sum += hitmap[val][o].size();
	      }
	    ntvsorted.push_back(sum);

	  }

	//	
      } // o

    int vtmp(0);
	// get KEmap indices by most energetic first, least last.
	for (auto it = KEmap.rbegin(); it!=KEmap.rend(); ++it) 
	  {
	    //	    int tval = it->second; // grab trkIDs in order, since they're sorted by KE
	    //	    int ke = it->first; // grab trkIDs in order, since they're sorted by KE
	    //	    const simb::MCParticle* part = pi_serv->TrackIDToParticle(tval);
	    
	    //	    std::cout << "TrackAnaCTStitch: KEmap cntr vtmp, Length ke, tval, pdg : "  << vtmp << ", " << ke <<", " << tval <<", " << part->PdgCode() << ", " << std::endl;

	    vtmp++;
	  }

    // get o trk indices by their hits. Most populous track first, least last.
    for (auto const o : fsort_indexes(ntvsorted))
      {
	int v(0);
	// get KEmap indices by longest trajectory first, least last.
	for (auto it = KEmap.rbegin(); it!=KEmap.rend(); ++it) 
	  {
	    int val = it->second; // grab trkIDs in order, since they're sorted by KE
	    //	    const simb::MCParticle* part = pi_serv->TrackIDToParticle(val);
	    //	    std::cout << "TrackAnaCTStitch: trk o, KEmap cntr v, KE val, pdg  hitmap[val][o].size(): "  << o <<", " << v << ", " << val <<", " << part->PdgCode() << ", " << hitmap[val][o].size() << std::endl;
	    rhistsStitched.fNTrkIdTrks3->Fill(o,v,hitmap[val][o].size());
	    v++;
	  }
	
      }

    // In how many o tracks did each trkId appear? Histo it. Would like it to be precisely 1.
    // Histo it vs. particle KE.
    flattener flat(NtrkIdsAll);
    std::vector <unsigned int> &v = flat;
    //auto const it ( std::unique(v.begin(),v.end()) ); // never use this it, perhaps.
    for (auto const val :  v)
      {
	if (val != (unsigned int)sim::NoParticleId)
	  {
	    const simb::MCParticle* part = pi_serv->TrackIdToParticle_P( val ); 
	    double T(part->E() - 0.001*part->Mass());
	    rhistsStitched.fNTrkIdTrks->Fill(std::count(v.begin(),v.end(),val));
	    rhistsStitched.fNTrkIdTrks2->Fill(std::count(v.begin(),v.end(),val),T);
	  }
	else
	  {
	    rhistsStitched.fNTrkIdTrks2->Fill(-1.0,0.0);
	  }
      }
 
  }

  void TrackAnaCT::endJob()
  //
  // Purpose: End of job.
  //
  {
    // Print summary.

    mf::LogInfo("TrackAnaCT") 
      << "TrackAnaCT statistics:\n"
      << "  Number of events = " << fNumEvent;

    // Fill efficiency histograms.

    for(std::map<int, MCHists>::const_iterator i = fMCHistMap.begin();
	i != fMCHistMap.end(); ++i) {
      const MCHists& mchists = i->second;
      effcalc(mchists.fHgstartx, mchists.fHmcstartx, mchists.fHestartx);
      effcalc(mchists.fHgstarty, mchists.fHmcstarty, mchists.fHestarty);
      effcalc(mchists.fHgstartz, mchists.fHmcstartz, mchists.fHestartz);
      effcalc(mchists.fHgendx, mchists.fHmcendx, mchists.fHeendx);
      effcalc(mchists.fHgendy, mchists.fHmcendy, mchists.fHeendy);
      effcalc(mchists.fHgendz, mchists.fHmcendz, mchists.fHeendz);
      effcalc(mchists.fHgtheta, mchists.fHmctheta, mchists.fHetheta);
      effcalc(mchists.fHgphi, mchists.fHmcphi, mchists.fHephi);
      effcalc(mchists.fHgtheta_xz, mchists.fHmctheta_xz, mchists.fHetheta_xz);
      effcalc(mchists.fHgtheta_yz, mchists.fHmctheta_yz, mchists.fHetheta_yz);
      effcalc(mchists.fHgmom, mchists.fHmcmom, mchists.fHemom);
      effcalc(mchists.fHglen, mchists.fHmclen, mchists.fHelen);
    }
  }

    // Stole this from online. Returns indices sorted by corresponding vector values.
    template <typename T>
      std::vector<size_t> TrackAnaCT::fsort_indexes(const std::vector<T> &v) {
  // initialize original index locations
      std::vector<size_t> idx(v.size());
      for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
            // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2];}); // Want most occupied trks first. was <, EC, 21-July.
      return idx;
    }


}
