#include "dunereco/BlipRecoDUNE/Alg/BlipRecoAlg.h"

namespace blip {

  //###########################################################
  // Constructor
  //###########################################################
  BlipRecoAlg::BlipRecoAlg( fhicl::ParameterSet const& pset )
  : fGeom { *lar::providerFrom<geo::Geometry>() }
  {
    this->reconfigure(pset);
   
    auto const& detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    kLArDensity           = detProp.Density();
    kNominalEfield        = detProp.Efield();
    kDriftVelocity        = detProp.DriftVelocity(detProp.Efield(0),detProp.Temperature()); 
    kTickPeriod           = clockData.TPCClock().TickPeriod();
    kNominalRecombFactor  = ModBoxRecomb(fCalodEdx,kNominalEfield);
    kWion                 = 1000./util::kGeVToElectrons;
   
    // -------------------------------------------------------------------
    // Determine number cryostats, TPC, planes, wires.
    //
    // Also cache all the X tick offsets so we don't have to keep re-calculating them 
    // for every single hit. Note that 'detProp.GetXTicksOffset()' does not make intuitive 
    // sense for cases with wireplanes aren't at X~0 (i.e., SBND).  We will account 
    // for that here so that our calculated drift times for hits makes sense.
    
    int kNumChannels = 0;

    // Loop over cryostats
    std::cout<<"NCryostats: "<<fGeom.Ncryostats()<<"\n";
    for(size_t cstat=0; cstat<fGeom.Ncryostats(); cstat++){
      auto const& cryoid = geo::CryostatID(cstat);
      std::cout<<"Cryostat "<<cstat<<" has "<<fGeom.NTPC(cryoid)<<" TPCs\n";

      // Loop TPCs in cryostat 'cstat'
      for(size_t tpc=0; tpc<fGeom.NTPC(cryoid); tpc++){
        auto const& tpcid = geo::TPCID(cryoid,tpc);
        
        float halfwidth   = fGeom.DetHalfWidth(tpcid);
        float halfheight  = fGeom.DetHalfHeight(tpcid);
        float length      = fGeom.DetLength(tpcid);
        auto  tpcffcenter = fGeom.GetTPCFrontFaceCenter(tpcid);
        
        std::cout<<"CRYOSTAT "<<cstat<<" / TPC "<<tpc<<"\n";
        printf("Front face center: (%f,%f,%f)\n",tpcffcenter.X(),tpcffcenter.Y(),tpcffcenter.Z());
        std::cout<<"DetHalfWidth: "<<halfwidth<<" cm, halfheight: "<<halfheight<<" cm, length: "<<length<<"\n";
        
        // Loop planes in TPC 'tpc'
        for(size_t pl=0; pl<fGeom.Nplanes(tpcid); pl++){
          auto const& planeid = geo::PlaneID(cstat,tpc,pl);
          kNumChannels += fGeom.Nwires(planeid);
          float offset = detProp.GetXTicksOffset(pl,tpc,cstat);
          std::cout<<"  PLANE "<<pl<<":  "<<fGeom.Nwires(planeid)<<" wires, XTicksOffset: "<<offset<<"\n";
          kXTicksOffsets[cstat][tpc][pl] = 0;

          if( fApplyXTicksOffset ) {
            //kXTicksOffsets[cstat][tpc][pl] = offset;
            // subtract out the geometric time offset added to account for the 
            // distance between plane0 and X=0. This is based on code in
            // lardataalg/DetectorInfo/DetectorPropertiesStandard.cxx 
            // (as of lardataalg v9_15_01)
            auto const& cryostat  = fGeom.Cryostat(geo::CryostatID(cstat));
            auto const& tpcgeom   = cryostat.TPC(tpc);
            auto const xyz        = tpcgeom.Plane(0).GetCenter();
            const double dir((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 : -1.0);
            float x_ticks_coefficient = kDriftVelocity*kTickPeriod;
            float goofy_offset = -xyz.X() / (dir * x_ticks_coefficient);
            std::cout<<"  Offset after geometric correction: "<<offset - goofy_offset<<"\n";
            kXTicksOffsets[cstat][tpc][pl] = offset - goofy_offset;
          
          } else {
            // we still want to correct for global trigger offset though:
            kXTicksOffsets[cstat][tpc][pl] = clockData.TriggerTime()/kTickPeriod;
          }
          
          // additional ad-hoc corrections supplied by user
          kXTicksOffsets[cstat][tpc][pl] += fTimeOffset[pl];
        }
      }
    }

    // initialize custom 'bad' channel list
    fBadChanMask       .resize(8256,false);
    fBadChanMaskPerEvt = fBadChanMask;
    if( fBadChanFile != "" ) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string fullname;
      sp.find_file(fBadChanFile,fullname);
      if (fullname.empty()) {
        throw cet::exception("Bad channel list not found");
      } else {
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        while (std::getline(inFile,line)) {
          if( line.find("#") != std::string::npos ) continue;
          std::istringstream ss(line);
          int ch1, ch2;
          ss >> ch1;
          if( !(ss >> ch2) ) ch2 = ch1;
          for(int i=ch1; i<=ch2; i++) fBadChanMask[i] = true;
        }
      }
    }
    int NBadChansFromFile     = std::count(fBadChanMask.begin(),fBadChanMask.end(),true);

    EvtBadChanCount = 0;

    printf("******************************************\n");
    printf("Initializing BlipRecoAlg...\n");
    printf("  - Efield: %.4f kV/cm\n",kNominalEfield);
    printf("  - Drift velocity: %.4f cm/us\n",kDriftVelocity);
    printf("  - using dE/dx: %.2f MeV/cm\n",fCalodEdx);
    printf("  - equiv. recomb: %.4f\n",kNominalRecombFactor);
    printf("  - custom bad chans: %i\n",NBadChansFromFile);
    printf("*******************************************\n");

    // create diagnostic histograms
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory hdir = tfs->mkdir("BlipRecoAlg");
   
    if( fSaveTree ) {
      matchTree = hdir.make<TTree>("matchtree","timematch tree");
      matchTree ->Branch("event",&_event,"event/I");
      matchTree ->Branch("run",&_run,"run/I");
      matchTree ->Branch("cryo",&_cryo,"cryo/I");
      matchTree ->Branch("tpc",&_tpc,"tpc/I");
      matchTree ->Branch("mc",&_mc,"mc/O");
      matchTree ->Branch("planeA",&_planeA,"planeA/I");
      matchTree ->Branch("planeB",&_planeB,"planeB/I");
      matchTree ->Branch("overlap",&_overlap,"overlap/F");
      matchTree ->Branch("dt",&_dt,"dt/F");
      matchTree ->Branch("dtfrac",&_dtfrac,"dtfrac/F");
      matchTree ->Branch("qA",&_qA,"qA/F");
      matchTree ->Branch("qB",&_qB,"qB/F");
      matchTree ->Branch("score",&_score,"score/F");
    }

    /*
    h_chanstatus     = hdir.make<TH1D>("chanstatus","Channel status for 'channels' list",5,0,5);
    h_chanstatus     ->GetXaxis()->SetBinLabel(1, "disconnected");
    h_chanstatus     ->GetXaxis()->SetBinLabel(2, "dead");
    h_chanstatus     ->GetXaxis()->SetBinLabel(3, "lownoise");
    h_chanstatus     ->GetXaxis()->SetBinLabel(4, "noisy");
    h_chanstatus     ->GetXaxis()->SetBinLabel(5, "good");
    h_hit_chanstatus     = hdir.make<TH1D>("hit_chanstatus","Channel status of hits",5,0,5);
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(1, "disconnected");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(2, "dead");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(3, "lownoise");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(4, "noisy");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(5, "good");
    */

    //h_hit_times       = hdir.make<TH1D>("hit_peaktime","Hit peaktimes",500,-5000,5000);
    //h_chan_nhits      = hdir.make<TH1D>("chan_nhits","Untracked hits;TPC readout channel;Total hits",kNumChannels,0,kNumChannels);
    //h_chan_nclusts    = hdir.make<TH1D>("chan_nclusts","Untracked isolated hits;TPC readout channel;Total clusts",kNumChannels,0,kNumChannels);
    h_clust_nwires    = hdir.make<TH1D>("clust_nwires","Clusters (pre-cut);Wires in cluster",100,0,100);
    h_clust_timespan  = hdir.make<TH1D>("clust_timespan","Clusters (pre-cut);Time span [ticks]",300,0,300);

    int qbins = 200;
    float qmax = 100;
    
      for(int i=0; i<kNplanes; i++) {
        if( i == fCaloPlane ) continue;
        h_clust_overlap[i]    = hdir.make<TH1D>(Form("p%i_clust_overlap",i),   Form("Plane %i clusters;Overlap fraction",i),101,0,1.01);
        h_clust_dtfrac[i]     = hdir.make<TH1D>(Form("p%i_clust_dtfrac",i),    Form("Plane %i clusters;Charge-weighted mean dT/RMS",i),150,-1.5,1.5);
        h_clust_dt[i]         = hdir.make<TH1D>(Form("p%i_clust_dt",i),        Form("Plane %i clusters;dT [ticks]",i),200,-10,10);
        h_clust_q[i]     = hdir.make<TH2D>(Form("p%i_clust_charge",i),  
          Form("Pre-cut, Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3} e-]",fCaloPlane,i),
          qbins,0,qmax,qbins,0,qmax);
          h_clust_q[i]->SetOption("colz");
        h_clust_qratio[i]   = hdir.make<TH1D>(Form("p%i_clust_qratio",i),        Form("Plane %i clusters;Charge ratio",i),100,0,1.01);
        h_clust_score[i]    = hdir.make<TH1D>(Form("p%i_clust_matchscore",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);
        h_clust_truematch_overlap[i]    = hdir.make<TH1D>(Form("p%i_clust_truematch_overlap",i),   Form("Plane %i clusters;Overlap fraction",i),101,0,1.01);
        h_clust_truematch_dt[i]         = hdir.make<TH1D>(Form("p%i_clust_truematch_dt",i),        Form("Plane %i clusters;dT [ticks]",i),200,-10,10);
        h_clust_truematch_dtfrac[i]     = hdir.make<TH1D>(Form("p%i_clust_truematch_dtfrac",i),    Form("Plane %i clusters;Charge-weighted mean dT/RMS",i),120,-3,3);
        h_clust_truematch_q[i]     = hdir.make<TH2D>(Form("p%i_clust_truematch_charge",i),  
          Form("Pre-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3} e-]",fCaloPlane,i),
          qbins,0,qmax,qbins,0,qmax);
          h_clust_truematch_q[i]->SetOption("colz");
        h_clust_truematch_qratio[i]   = hdir.make<TH1D>(Form("p%i_clust_truematch_qratio",i),        Form("Plane %i clusters;Charge ratio",i),100,0,1.01);
        h_clust_truematch_score[i]    = hdir.make<TH1D>(Form("p%i_clust_truematch_matchscore",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);
        h_nmatches[i]         = hdir.make<TH1D>(Form("p%i_nmatches",i),Form("Number of plane%i matches to single collection cluster",i),20,0,20);
    
    }//endloop over planes
    
    // Efficiency as a function of energy deposited on a wire
    h_recoWireEff_denom = hdir.make<TH1D>("recoWireEff_trueCount","Collection plane;Electron energy deposited on wire [MeV];Count",150,0,1.5);
    h_recoWireEff_num   = hdir.make<TH1D>("recoWireEff","Collection plane;Electron energy deposited on wire [MeV];Hit reco efficiency",150,0,1.5);
    h_recoWireEffQ_denom = hdir.make<TH1D>("recoWireEffQ_trueCount","Collection plane;Charge deposited on wire [e-];Count",80,0,20000);
    h_recoWireEffQ_num   = hdir.make<TH1D>("recoWireEffQ","Collection plane;Charge deposited on wire [e-];Hit reco efficiency",80,0,20000);
    

  }
  
  //--------------------------------------------------------------  
  //Destructor
  BlipRecoAlg::~BlipRecoAlg()
  {
  }
  
  
  //###########################################################
  // Reconfigure fcl parameters
  //###########################################################
  void BlipRecoAlg::reconfigure( fhicl::ParameterSet const& pset ){
    
    fDebugMode          = pset.get<bool>     ("DebugMode",         false);

    fHitProducer        = pset.get<std::string>   ("HitProducer",       "gaushit");
    fTrkProducer        = pset.get<std::string>   ("TrkProducer",       "pandora");
    fGeantProducer      = pset.get<std::string>   ("GeantProducer",     "largeant");
    fSimDepProducer     = pset.get<std::string>   ("SimEDepProducer",   "ionization");
    fSimChanProducer    = pset.get<std::string>   ("SimChanProducer",   "driftWC:simpleSC");
    fSimGainFactor      = pset.get<float>         ("SimGainFactor",     -9);
    fTrueBlipMergeDist  = pset.get<float>         ("TrueBlipMergeDist", 0.3);
    fMaxHitTrkLength    = pset.get<float>               ("MaxHitTrkLength", 5);
    fDoHitFiltering     = pset.get<bool>                ("DoHitFiltering",  false);
    fMaxHitMult         = pset.get<int>                 ("MaxHitMult",      10);
    fMaxHitAmp          = pset.get<float>               ("MaxHitAmp",       200);  
    fMinHitAmp          = pset.get<std::vector<float>>  ("MinHitAmp",       {-99e9,-99e9,-99e9});
    fMaxHitRMS          = pset.get<std::vector<float>>  ("MaxHitRMS",       { 99e9, 99e9, 99e9});
    fMinHitRMS          = pset.get<std::vector<float>>  ("MinHitRMS",       {-99e9,-99e9,-99e9});
    fMaxHitRatio        = pset.get<std::vector<float>>  ("MaxHitRatio",     { 99e9, 99e9, 99e9});
    fMinHitRatio        = pset.get<std::vector<float>>  ("MinHitRatio",     {-99e9,-99e9,-99e9});
    fMaxHitGOF          = pset.get<std::vector<float>>  ("MaxHitGOF",       { 99e9, 99e9, 99e9});
    fMinHitGOF          = pset.get<std::vector<float>>  ("MinHitGOF",       {-99e9,-99e9,-99e9});
    
    fHitClustWidthFact  = pset.get<float>         ("HitClustWidthFact", 5.0);
    fHitClustWireRange  = pset.get<int>           ("HitClustWireRange", 1);
    fMaxWiresInCluster  = pset.get<int>           ("MaxWiresInCluster", 10);
    fMaxClusterSpan     = pset.get<float>         ("MaxClusterSpan",    30);
    fMinClusterCharge   = pset.get<float>         ("MinClusterCharge",  300);
    fMaxClusterCharge   = pset.get<float>         ("MaxClusterCharge",  12e6);
    
    fApplyXTicksOffset  = pset.get<bool>          ("ApplyXTicksOffset",     true);
    fTimeOffset         = pset.get<std::vector<float>>("TimeOffset", {0.,0.,0.});
    fMatchMinOverlap    = pset.get<float>         ("ClustMatchMinOverlap",  0.5 );
    fMatchSigmaFact     = pset.get<float>         ("ClustMatchSigmaFact",   1.0);
    fMatchMaxTicks      = pset.get<float>         ("ClustMatchMaxTicks",    5.0 );
    fMatchQDiffLimit    = pset.get<float>         ("ClustMatchQDiffLimit",  15e3);
    fMatchMinQRatio     = pset.get<float>         ("ClustMatchMinQRatio",   0.25);
    
    fMinMatchedPlanes   = pset.get<int>           ("MinMatchedPlanes",    2);
    fPickyBlips         = pset.get<bool>          ("PickyBlips",          false);
    fApplyTrkCylinderCut= pset.get<bool>          ("ApplyTrkCylinderCut", false);
    fCylinderRadius     = pset.get<float>         ("CylinderRadius",      15);
    
    fCaloAlg            = new calo::CalorimetryAlg( pset.get<fhicl::ParameterSet>("CaloAlg") );
    fCaloPlane          = pset.get<int>           ("CaloPlane",           2);
    fCalodEdx           = pset.get<float>         ("CalodEdx",            2.8);
    fLifetimeCorr       = pset.get<bool>          ("LifetimeCorrection",  false);
    fSCECorr            = pset.get<bool>          ("SCECorrection",       false);
    fYZUniformityCorr   = pset.get<bool>          ("YZUniformityCorrection",true);
    fModBoxA            = pset.get<float>         ("ModBoxA",             0.93);
    fModBoxB            = pset.get<float>         ("ModBoxB",             0.212);
    
    fVetoBadChannels    = pset.get<bool>          ("VetoBadChannels",     true);
    fBadChanProducer    = pset.get<std::string>   ("BadChanProducer",     "nfspl1:badchannels");
    fBadChanFile        = pset.get<std::string>   ("BadChanFile",         "");
    fMinDeadWireGap     = pset.get<int>           ("MinDeadWireGap",      1);
    
    fSaveTree           = pset.get<bool>          ("SaveTree",            true);

  }



  //###########################################################
  // Main reconstruction procedure.
  //
  // This function does EVERYTHING. The resulting collections of 
  // blip::HitClusts and blip::Blips can then be retrieved after
  // this function is run.
  //###########################################################
  void BlipRecoAlg::RunBlipReco( const art::Event& evt ) {
  
    std::cout<<"\n"
    <<"=========== BlipRecoAlg =========================\n"
    <<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"\n";
  
    //=======================================
    // Reset things
    //=======================================
    blips.clear();
    hitclust.clear();
    hitinfo.clear();
    pinfo.clear();
    trueblips.clear();
    EvtBadChanCount = 0;
  
    //=======================================
    // Get data products for this event
    //========================================
    
    // --- detector properties
    auto const& SCE_provider        = lar::providerFrom<spacecharge::SpaceChargeService>();
    auto const& chanFilt            = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    auto const& detProp             = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& clockData           = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    //auto const& lifetime_provider   = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
    //auto const& tpcCalib_provider   = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
  
    // -- geometry
    art::ServiceHandle<geo::Geometry> geom;

    // -- G4 particles
    art::Handle< std::vector<simb::MCParticle> > pHandle;
    std::vector<art::Ptr<simb::MCParticle> > plist;
    if (evt.getByLabel(fGeantProducer,pHandle))
      art::fill_ptr_vector(plist, pHandle);
 
    // -- SimEnergyDeposits
    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
    if (evt.getByLabel(fSimDepProducer,sedHandle)) 
      art::fill_ptr_vector(sedlist, sedHandle);
    
    // -- SimChannels (usually dropped in reco)
    art::Handle<std::vector<sim::SimChannel> > simchanHandle;
    std::vector<art::Ptr<sim::SimChannel> > simchanlist;
    if (evt.getByLabel(fSimChanProducer,simchanHandle)) 
      art::fill_ptr_vector(simchanlist, simchanHandle);

    // -- hits (from input module, usually track-masked subset of gaushit)
    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitProducer,hitHandle))
      art::fill_ptr_vector(hitlist, hitHandle);

    // -- hits (from gaushit), these are used in truth-matching of hits
    art::Handle< std::vector<recob::Hit> > hitHandleGH;
    std::vector<art::Ptr<recob::Hit> > hitlistGH;
    if (evt.getByLabel("gaushit",hitHandleGH))
      art::fill_ptr_vector(hitlistGH, hitHandleGH);

    // -- tracks
    art::Handle< std::vector<recob::Track> > tracklistHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrkProducer,tracklistHandle))
      art::fill_ptr_vector(tracklist, tracklistHandle);
  
    // -- associations
    art::FindManyP<recob::Track> fmtrk(hitHandle,evt,fTrkProducer);
    art::FindManyP<recob::Track> fmtrkGH(hitHandleGH,evt,fTrkProducer);
   
    // -- backtracker
    art::ServiceHandle<cheat::BackTrackerService> btService;

    std::cout
    <<"Found "<<hitlist.size()<<" hits from "<<fHitProducer<<"\n"
    <<"Found "<<tracklist.size()<<" tracks from "<<fTrkProducer<<"\n"
    <<"Found "<<plist.size()<<" MC particles from "<<fGeantProducer<<"\n"
    <<"Found "<<sedlist.size()<<" SimEnergyDeposits from "<<fSimDepProducer<<"\n"
    <<"Found "<<simchanlist.size()<<" SimChannels from "<<fSimChanProducer<<"\n";
  
    //====================================================
    // Update map of bad channels for this event
    //====================================================
    if( fVetoBadChannels ) {
      fBadChanMaskPerEvt = fBadChanMask;
      if( fBadChanProducer != "" ) { 
        std::vector<int> badChans;
        art::Handle< std::vector<int>> badChanHandle;
        if( evt.getByLabel(fBadChanProducer, badChanHandle))
          badChans = *(badChanHandle);
        for(auto& ch : badChans ) {
          EvtBadChanCount++;
          fBadChanMaskPerEvt[ch] = true;
        }
      }
    }
   
    //====================================================
    // Prep the particle inventory service for MC+overlay
    //====================================================
    if( evt.isRealData() && plist.size() ) {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      pi_serv->Rebuild(evt);
      pi_serv->provider()->PrepParticleList(evt);
    }
    

    //===============================================================
    // Map of each hit to its gaushit index (needed if the provided
    // hit collection is some filtered subset of gaushit, in order to
    // use gaushitTruthMatch later on)
    //===============================================================
    std::map< int, int > map_gh;
    // if input collection is already gaushit, this is trivial
    if( fHitProducer == "gaushit" ) {
      for(auto& h : hitlist ) map_gh[h.key()] = h.key(); 
    // ... but if not, find the matching gaushit. There's no convenient
    // hit ID, so we must loop through and compare channel/time (ugh)
    } else {
      std::map<int,std::vector<int>> map_chan_ghid;
      for(auto& gh : hitlistGH ) map_chan_ghid[gh->Channel()].push_back(gh.key());
      for(auto& h : hitlist ) {
        for(auto& igh : map_chan_ghid[h->Channel()]){
          if( hitlistGH[igh]->PeakTime() != h->PeakTime() ) continue;
          map_gh[h.key()] = igh;
          break;
        }
      }
    }
    
   
    //=====================================================
    // Record PDG for every G4 Track ID
    //=====================================================
    std::map<int,int> map_g4trkid_pdg;
    std::map<int, std::map<int,double> > map_g4trkid_chan_energy;
    std::map<int, std::map<int,double> > map_g4trkid_chan_charge;
    for(size_t i = 0; i<plist.size(); i++) 
      map_g4trkid_pdg[plist[i]->TrackId()] = plist[i]->PdgCode();
    

    //======================================================
    // Use SimChannels to make a map of the collected charge
    // for every G4 particle, instead of relying on the TDC-tick
    // matching that's done by BackTracker's other functions
    //======================================================
    std::map<int,double> map_g4trkid_charge;
    for(auto const &chan : simchanlist ) {
      if( fGeom.View(chan->Channel()) != fCaloPlane ) continue;
      //std::map<int,double> map_g4trkid_perWireEnergyDep;
      for(auto const& tdcide : chan->TDCIDEMap() ) {
        for(auto const& ide : tdcide.second) {
          // In SBND, when 'KeepEMShowerDaughters' is disabled, any
          // energy deposits from secondary particles are folded into
          // their most upstream parent, and the trackID is set to
          // negative of that parent.
          if( ide.trackID < 0 ) continue;
          double ne = ide.numElectrons;
          
          // ####################################################
          // # behavior in MicroBooNE as of Nov 2022            #
          // ####################################################
          // WireCell's detsim implements its gain "fudge factor" 
          // by scaling the SimChannel electrons (DocDB 31089)
          // instead of the electronics gain. So we need to correct 
          // for this effect to get accurate count of 'true' 
          // electrons collected on this channel.
          if( fSimGainFactor > 0 ) ne /= fSimGainFactor;
          // ####################################################
          
          map_g4trkid_charge[ide.trackID] += ne;
         
          // keep track of charge deposited per wire for efficiency plots
          map_g4trkid_chan_charge[ide.trackID][chan->Channel()] += ne; 
          if( abs(map_g4trkid_pdg[ide.trackID]) == 11 ) 
            map_g4trkid_chan_energy[ide.trackID][chan->Channel()] += ide.energy;
        
        }
      }
    }


    for(auto& m : map_g4trkid_chan_energy ) {
      for(auto& mm : m.second ) {
        if( mm.second > 0 ) h_recoWireEff_denom->Fill(mm.second);
      }
    } 
    for(auto& m : map_g4trkid_chan_charge ) {
      for(auto& mm : m.second ) {
        if( mm.second > 0 ) h_recoWireEffQ_denom->Fill(mm.second);
      }
    }
    
   

    //==================================================
    // Use G4 information to determine the "true" blips in this event.
    //==================================================
    if( plist.size() ) {
      pinfo.resize(plist.size());
      for(size_t i = 0; i<plist.size(); i++){
        BlipUtils::FillParticleInfo( *plist[i], pinfo[i], sedlist, fCaloPlane);
        if( map_g4trkid_charge[pinfo[i].trackId] ) pinfo[i].numElectrons = (int)map_g4trkid_charge[pinfo[i].trackId];
        pinfo[i].index = i;
      }
      BlipUtils::MakeTrueBlips(pinfo, trueblips);
      BlipUtils::MergeTrueBlips(trueblips, fTrueBlipMergeDist);
    }


    //=======================================
    // Map track IDs to the index in the vector
    //=======================================
    std::map<size_t,size_t> map_trkid_index;
    for(size_t i=0; i<tracklist.size(); i++) 
      map_trkid_index[tracklist.at(i)->ID()] = i;

    //=======================================
    // Fill vector of hit info
    //========================================
    hitinfo.resize(hitlist.size());
    std::map<int, std::map<int, std::map<int,std::vector<int> >>> cryo_tpc_plane_hitsMap;
    int nhits_untracked = 0;

    for(size_t i=0; i<hitlist.size(); i++){
      auto const& thisHit = hitlist[i];
      auto const& wireid  = thisHit->WireID();
      int   chan    = thisHit->Channel();
      int   cstat   = wireid.Cryostat;
      int   tpc     = wireid.TPC;
      int   plane   = wireid.Plane;
      int   wire    = wireid.Wire;
        
      hitinfo[i].hitid        = i;
      hitinfo[i].cryo         = cstat;
      hitinfo[i].tpc          = tpc;
      hitinfo[i].plane        = plane;
      hitinfo[i].chan         = chan;
      hitinfo[i].wire         = wire;
      hitinfo[i].amp          = thisHit->PeakAmplitude();
      hitinfo[i].rms          = thisHit->RMS();
      hitinfo[i].integralADC  = thisHit->Integral();
      hitinfo[i].sigmaintegral = thisHit->SigmaIntegral();
      hitinfo[i].sumADC       = thisHit->SummedADC();
      hitinfo[i].charge       = fCaloAlg->ElectronsFromADCArea(thisHit->Integral(),plane);
      hitinfo[i].gof          = thisHit->GoodnessOfFit() / thisHit->DegreesOfFreedom();
      hitinfo[i].peakTime     = thisHit->PeakTime();
      hitinfo[i].driftTime    = thisHit->PeakTime()-kXTicksOffsets[cstat][tpc][plane]; //detProp.GetXTicksOffset(wireid);

      if( plist.size() ) {
        int truthid       = -9;
        float truthidfrac = -9;
        std::vector<sim::TrackIDE> trackIDEs = art::ServiceHandle<cheat::BackTrackerService>()->HitToTrackIDEs(clockData,thisHit);
        float maxe = 0;
        float ne = 0;
        for(size_t i = 0; i < trackIDEs.size(); ++i){
          ne += (float)trackIDEs[i].numElectrons;
          if( trackIDEs[i].energy > maxe ) {
            maxe = trackIDEs[i].energy;
            truthidfrac = trackIDEs[i].energyFrac;
            truthid = trackIDEs[i].trackID;
          }
        }
        
        // Save the results
        hitinfo[i].g4trkid  = truthid;
        hitinfo[i].g4pdg    = map_g4trkid_pdg[truthid];
        hitinfo[i].g4frac   = truthidfrac; 
        
        if( map_g4trkid_chan_energy[hitinfo[i].g4trkid][chan] > 0 ) 
          h_recoWireEff_num->Fill(map_g4trkid_chan_energy[hitinfo[i].g4trkid][chan]);
        if( map_g4trkid_chan_charge[hitinfo[i].g4trkid][chan] > 0 ) 
          h_recoWireEffQ_num->Fill(map_g4trkid_chan_charge[hitinfo[i].g4trkid][chan]);

      }//endif MC
      

      // find associated track
      if( fHitProducer == "gaushit" && fmtrk.isValid() ) {
        if(fmtrk.at(i).size()) hitinfo[i].trkid = fmtrk.at(i)[0]->ID();
      
      // if the hit collection didn't have associations made
      // to the tracks, try gaushit instead
      } else if ( fmtrkGH.isValid() && map_gh.size() ) {
        int gi = map_gh[i];
        if (fmtrkGH.at(gi).size()) hitinfo[i].trkid= fmtrkGH.at(gi)[0]->ID(); 
      }

      // add to the map
      cryo_tpc_plane_hitsMap[cstat][tpc][plane].push_back(i);
      if( hitinfo[i].trkid < 0 ) nhits_untracked++;

    }//endloop over hits
    

    //=================================================================
    // Blip Reconstruction
    //================================================================
    //  
    //  Procedure
    //  [x] Look for hits that were not included in a track 
    //  [x] Filter hits based on hit width, etc
    //  [x] Merge together closely-spaced hits on same wires and adjacent wires
    //  [x] Plane-to-plane time matching
    //  [x] Wire intersection check to get XYZ
    //  [x] Create "blip" object

    // Create a series of masks that we'll update as we go along
    std::vector<bool> hitIsTracked(hitlist.size(),  false);
    std::vector<bool> hitIsGood(hitlist.size(),     true);
    std::vector<bool> hitIsClustered(hitlist.size(),false);
    
    
    // Basic track inclusion cut: exclude hits that were tracked
    for(size_t i=0; i<hitlist.size(); i++){
      if( hitinfo[i].trkid < 0 ) continue;
      auto it = map_trkid_index.find(hitinfo[i].trkid);
      if( it == map_trkid_index.end() ) continue;
      int trkindex = it->second;
      if( tracklist[trkindex]->Length() > fMaxHitTrkLength ) {
        hitIsTracked[i] = true;
        hitIsGood[i] = false;
      }
    }
        

    // Filter based on hit properties. For hits that are a part of
    // multi-gaussian fits (multiplicity > 1), need to re-think this.
    if( fDoHitFiltering ) {
      for(size_t i=0; i<hitlist.size(); i++){
        if( !hitIsGood[i] ) continue;
        hitIsGood[i] = false;
        auto& hit = hitlist[i];
        int plane = hit->WireID().Plane;
        if( hitinfo[i].gof        <= fMinHitGOF[plane] ) continue;
        if( hitinfo[i].gof        >= fMaxHitGOF[plane] ) continue;
        if( hit->RMS()            <= fMinHitRMS[plane] ) continue;
        if( hit->RMS()            >= fMaxHitRMS[plane] ) continue;
        if( hit->PeakAmplitude()  <= fMinHitAmp[plane] ) continue;
        if( hit->PeakAmplitude()  >= fMaxHitAmp )        continue;
        if( hit->Multiplicity()   >= fMaxHitMult )       continue;
        //float hit_ratio = hit->RMS() / hit->PeakAmplitude();
        //if( hit_ratio             < fMinHitRatio[plane] ) continue;
        //if( hit_ratio             > fMaxHitRatio[plane] ) continue;
        
        // we survived the gauntlet of cuts -- hit is good!
        hitIsGood[i] = true;
      }
    }

    // ---------------------------------------------------
    // Hit clustering
    // ---------------------------------------------------
    std::map<int,std::map<int,std::map<int,std::vector<int> >>> cryo_tpc_planeclustsmap;
    //std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
   
    for(auto const& tpc_plane_hitsMap : cryo_tpc_plane_hitsMap ) {
    
      for(auto const& plane_hitsMap : tpc_plane_hitsMap.second ) {
        //std::cout<<"Looking at TPC "<<plane_hitsMap.first<<", which has hits appearing in "<<plane_hitsMap.second.size()<<" planes\n";

        for(auto const& planehits : plane_hitsMap.second){
          //std::cout<<"Looking at TPC "<<plane_hitsMap.first<<", plane "<<planehits.first<<", which has "<<planehits.second.size()<<" hits\n"; 

          for(auto const& hi : planehits.second ){
            
            // skip hits flagged as bad, or already clustered
            if( !hitIsGood[hi] || hitIsClustered[hi] ) continue; 
            
            // initialize a new cluster with this hit as seed
            std::vector<blip::HitInfo> hitinfoVec;
            std::set<int> hitIDs;
            
            hitinfoVec    .push_back(hitinfo[hi]);
            hitIDs        .insert(hi);
            int startWire = hitinfo[hi].wire;
            int endWire   = hitinfo[hi].wire;
            hitIsClustered[hi] = true;

            // see if we can add other hits to it; continue until 
            // no new hits can be lumped in with this clust
            int hitsAdded;
            do{
              hitsAdded = 0;  
              for(auto const& hj : planehits.second ) {
                
                if( !hitIsGood[hj] || hitIsClustered[hj] ) continue; 

                // skip hits outside overall cluster wire range
                int w1 = hitinfo[hj].wire - fHitClustWireRange;
                int w2 = hitinfo[hj].wire + fHitClustWireRange;
                if( w2 < startWire    || w1 > endWire ) continue;
                
                // check for proximity with every other hit added
                // to this cluster so far
                for(auto const& hii : hitIDs ) {

                  if( hitinfo[hii].wire > w2 ) continue;
                  if( hitinfo[hii].wire < w1 ) continue;
                  
                  float t1 = hitinfo[hj].peakTime;
                  float t2 = hitinfo[hii].peakTime;
                  float rms_sum = (hitinfo[hii].rms + hitinfo[hj].rms);
                  if( fabs(t1-t2) > fHitClustWidthFact * rms_sum ) continue;

                  hitinfoVec.push_back(hitinfo[hj]);
                  startWire = std::min( hitinfo[hj].wire, startWire );
                  endWire   = std::max( hitinfo[hj].wire, endWire );
                  hitIDs.insert(hj);
                  hitIsClustered[hj] = true;
                  hitsAdded++;
                  break;
                }
              

              }
            } while ( hitsAdded!=0 );
            
            blip::HitClust hc = BlipUtils::MakeHitClust(hitinfoVec);
            float span = hc.EndTime - hc.StartTime;
            h_clust_nwires->Fill(hc.NWires);
            h_clust_timespan->Fill(span);
            
            // basic cluster checks
            if( span      <= 0                )   continue;
            if( span      > fMaxClusterSpan   )   continue;
            if( hc.NWires > fMaxWiresInCluster )  continue;
            if( hc.Charge < fMinClusterCharge )   continue;
            if( hc.Charge > fMaxClusterCharge )   continue;
            
            //std::cout<<"Making a new cluster on plane "<<planehits.first<<"\n";
            //std::cout<<"span "<<span<<" ticks, "<<hc.NWires<<" wires, "<<hc.Charge<<" electrons\n";
           
            // Exclude cluster if it is *entirely* on bad channels
            if( fVetoBadChannels ) {
              int nbadchanhits = 0;
              for(auto const& hitID : hc.HitIDs ) {
                int chan = hitinfo[hitID].chan;
                if( chanFilt.Status(chan) < 4 ||
                  fBadChanMaskPerEvt[chan] ) nbadchanhits++;
              }
              if( nbadchanhits == hc.NHits ) continue;
            }
            
            // measure wire separation to nearest dead region
            // (0 = directly adjacent)
            for(size_t dw=1; dw<=5; dw++){
              int  w1   = hc.StartWire-dw;
              int  w2   = hc.EndWire+dw;
              bool flag = false;
              // treat edges of wireplane as "dead"
              //if( w1 < 0 || w2 >= (int)fGeom.Nwires(hc.Plane) )
              if( w1 < 0 || w2 >= (int)fGeom.Nwires(geo::PlaneID(0,hc.TPC,hc.Plane)))
                flag=true;
              //otherwise, use channel filter service
              else {
                int ch1 = fGeom.PlaneWireToChannel(geo::WireID(0,hc.TPC,hc.Plane,w1));
                int ch2 = fGeom.PlaneWireToChannel(geo::WireID(0,hc.TPC,hc.Plane,w2));
                if( chanFilt.Status(ch1)<2 ) flag=true;
                if( chanFilt.Status(ch2)<2 ) flag=true;
              }
              if( flag ) { hc.DeadWireSep = dw-1; break; }
            }
            //std::cout<<"DeadWireSep "<<hc.DeadWireSep<<"\n";
           
            // veto this cluster if the gap between it and the
            // nearest dead wire (calculated above) isn't big enough
            if( fMinDeadWireGap > 0 && hc.DeadWireSep < fMinDeadWireGap ) continue;
            
            // **************************************
            // assign the ID, then go back and encode this 
            // cluster ID into the hit information
            // **************************************
            int idx = (int)hitclust.size();
            hc.ID = idx;
            //tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(idx);
            cryo_tpc_planeclustsmap[hc.Cryostat][hc.TPC][hc.Plane].push_back(idx);
            for(auto const& hitID : hc.HitIDs) hitinfo[hitID].clustid = hc.ID;
            // ... and find the associated truth-blip
            if( hc.G4IDs.size() ) {
              for(size_t j=0; j< trueblips.size(); j++){
                //std::cout<<"Looking at true blip "<<j<<" which has LeadG4ID of "<<trueblips[j].LeadG4ID<<"\n";
                if( hc.G4IDs.count(trueblips[j].LeadG4ID)) {
                  hc.isTruthMatched = true;
                  hc.EdepID         = trueblips[j].ID;
                  break;
                }
              }
            }
            
            // finally, add the finished cluster to the stack
            hitclust.push_back(hc);

          }
        }//loop over planes
      }//loop over TPCs
    }//loop over cryostats
    //std::cout<<"All done with clustering\n";
     


    // =============================================================================
    // Plane matching and 3D blip formation
    // =============================================================================

    // --------------------------------------
    // Method 1A: Require match between calo plane ( typically collection) and
    //            1 or 2 induction planes. For every hitclust on the calo plane,
    //            do the following:
    //              1. Loop over hitclusts in one of the other planes (same TPC)
    //              3. Find closest-matched clust and add it to the histclust group
    //              4. Repeat for remaining plane(s)
    
    
    float _matchQDiffLimit= (fMatchQDiffLimit <= 0 ) ? std::numeric_limits<float>::max() : fMatchQDiffLimit;
    //float _matchMinQRatio = (fMatchMinQRatio  <= 0 ) ? 0 : fMatchMinQRatio;
    
    for(auto& cryoMap : cryo_tpc_planeclustsmap ) { // cryostat loop 
    auto cryo = cryoMap.first;
    auto& tpcMaps = cryoMap.second;
    //std::cout<<"Looking in cryostat "<<cryo<<", which was clusters in "<<tpcMaps.size()<<" TPCs...\n";
  
    for(auto& tpcMap : tpcMaps ) { // loop on TPCs
      auto tpc = tpcMap.first;
      auto& planeMap = tpcMap.second;
      //std::cout<<"Performing cluster matching in TPC "<<tpc<<", which has clusters in "<<planeMap.size()<<" planes\n";

      if( planeMap.find(fCaloPlane) != planeMap.end() ){
        int   planeA              = fCaloPlane;
        auto&  hitclusts_planeA   = planeMap[planeA];
        //std::cout<<"using plane "<<fCaloPlane<<" as reference/calo plane ("<<planeMap[planeA].size()<<" clusts)\n";
        for(auto& i : hitclusts_planeA ) {
          auto& hcA = hitclust[i];

          // initiate hit-cluster group
          std::vector<blip::HitClust> hcGroup;
          hcGroup.push_back(hcA);

          // for each of the other planes, make a map of potential matches
          std::map<int, std::set<int>> cands;
         
          // map of cluster ID <--> match metrics
          std::map<int, float> map_clust_dtfrac;
          std::map<int, float> map_clust_dt;
          std::map<int, float> map_clust_overlap;
          std::map<int, float> map_clust_score;

          // ---------------------------------------------------
          // loop over other planes
          for(auto&  hitclusts_planeB : planeMap ) {
            int planeB = hitclusts_planeB.first;
            if( planeB == planeA ) continue;
            
            //std::cout<<"  checking clusters on plane "<<planeB<<"\n";

            // Loop over all non-matched clusts on this plane
            for(auto const& j : hitclusts_planeB.second ) {
              auto& hcB = hitclust[j];
              if( hcB.isMatched ) continue;
              
              // *******************************************
              // Check that the two central wires intersect
              // *******************************************
              const geo::WireID wireA(cryo,tpc,planeA,hcA.CenterWire);
              const geo::WireID wireB(cryo,tpc,planeB,hcB.CenterWire);
             
              geo::Point_t xyz;
              if( !art::ServiceHandle<geo::Geometry>()
                ->WireIDsIntersect(wireA,wireB,xyz)) continue;
              // Save intersect location, so we don't have to
              // make another call to the Geometry service later
              hcA.IntersectLocations[hcB.ID] = xyz; //xloc;
              hcB.IntersectLocations[hcA.ID] = xyz; //xloc;

              // ***********************************
              // Calculate the cluster overlap
              // ***********************************
              float overlapFrac = BlipUtils::CalcHitClustsOverlap(hcA,hcB);

              // *******************************************
              // Calculate time difference for start/end, and
              // check that Q-weighted means are comparable
              // *******************************************
              float dt_start  = (hcB.StartTime - hcA.StartTime);
              float dt_end    = (hcB.EndTime   - hcA.EndTime);
              float dt        = ( fabs(dt_start) < fabs(dt_end) ) ? dt_start : dt_end;
              float sigmaT    = std::sqrt(pow(hcA.RMS,2)+pow(hcB.RMS,2));
              float dtfrac    = (hcB.Time - hcA.Time) / sigmaT;

              // *******************************************
              // Check relative charge between clusters
              // *******************************************
              float qdiff     = fabs(hcB.Charge-hcA.Charge);
              float ratio     = std::min(hcA.Charge,hcB.Charge)/std::max(hcA.Charge,hcB.Charge);

              // *******************************************
              // Combine metrics into a consolidated score
              // *******************************************
              float score = 0;
              if( overlapFrac ) score=overlapFrac*ratio*exp(-fabs(dt)/float(fMatchMaxTicks));

              // If both clusters are matched to the same MC truth particle,
              // set flag to fill special diagnostic histograms...
              bool trueFlag = (hcA.isTruthMatched && (hcA.EdepID == hcB.EdepID)) ? true : false;
             
              // Diagnostic histograms
              h_clust_overlap[planeB] ->Fill(overlapFrac);
              h_clust_dt[planeB]      ->Fill(dt);
              h_clust_dtfrac[planeB]->Fill(dtfrac);
              double q1=hcA.Charge/1000.;
              double q2=hcB.Charge/1000.;
              h_clust_q[planeB]       ->Fill( q1, q2 );
              if( score > 0 ) h_clust_score[planeB]   ->Fill(score);
              if( trueFlag ) {
                h_clust_truematch_overlap[planeB] ->Fill(overlapFrac);
                h_clust_truematch_dt[planeB]      ->Fill(dt);
                h_clust_truematch_dtfrac[planeB]  ->Fill(dtfrac);
                h_clust_truematch_q[planeB]       ->Fill(q1,q2);
                if( score > 0 ) h_clust_truematch_score[planeB]  ->Fill(score);
              }
              
              if( fSaveTree ) {
                _run      = evt.id().run();
                _event    = evt.id().event();
                _mc       = trueFlag;
                _cryo     = cryo;
                _tpc      = tpc;
                _planeA   = planeA;
                _planeB   = planeB;
                _overlap  = overlapFrac;
                _dt       = dt;
                _dtfrac   = dtfrac;
                _qA       = q1;
                _qB       = q2;
                _score    = score;
                matchTree->Fill(); 
              }
              
              if( overlapFrac   < fMatchMinOverlap  ) continue;
              if( fabs(dt)      > fMatchMaxTicks    ) continue;
              if( fabs(dtfrac)  > fMatchSigmaFact   ) continue;
              if( qdiff         > _matchQDiffLimit  ) {
                if(trueFlag) h_clust_truematch_qratio[planeB]  ->Fill(ratio);
                h_clust_qratio[planeB]  ->Fill(ratio);
                if( ratio < fMatchMinQRatio ) continue;
              }
                
              // **************************************************
              // We made it through the cuts -- the match is good!
              // **************************************************
              map_clust_dt[j]       = dt;
              map_clust_dtfrac[j]   = dtfrac;
              map_clust_overlap[j]  = overlapFrac;
              map_clust_score[j]    = score;
              cands[planeB]         .insert(j);
            
            }
              
          }//endloop over other planes

          // ---------------------------------------------------
          // loop over the candidates found on each plane
          // and select the one with the largest score
          if( cands.size() ) {
            for(auto& c : cands ) {
              int plane = c.first;
              h_nmatches[plane]->Fill((int)c.second.size());
              float bestScore   = -9;
              int   bestID      = -9;
              for(auto cid : c.second) {
                if( map_clust_score[cid] > bestScore ) {
                  bestScore = map_clust_score[cid];
                  bestID = cid;
                }
              }
              if( bestID >= 0 ) hcGroup.push_back(hitclust[bestID]);

            }

            // ----------------------------------------
            // make our new blip, but if it isn't valid, forget it and move on
            blip::Blip newBlip = BlipUtils::MakeBlip(hcGroup,detProp,clockData);
            if( !newBlip.isValid ) continue;
            if( newBlip.NPlanes < fMinMatchedPlanes ) continue;
            
            // ---------------------------------------
            // does this qualify as a "picky" blip?
            bool picky = ( newBlip.NPlanes >= 3 && newBlip.SigmaYZ < 1. );

            // ----------------------------------------
            // save matching information
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].isMatched = true;
              for(auto hit : hitclust[hc.ID].HitIDs) hitinfo[hit].ismatch = true;
            }//end loop over clusters

            if( fPickyBlips && !picky ) continue;
            
            
            // ----------------------------------------
            // apply cylinder cut 
            for(auto& trk : tracklist ){
              if( trk->Length() < fMaxHitTrkLength ) continue;
              auto& a = trk->Vertex();
              auto& b = trk->End();
              TVector3 p1(a.X(), a.Y(), a.Z() );
              TVector3 p2(b.X(), b.Y(), b.Z() );
              // TO-DO: if this track starts or ends at a TPC boundary, 
              // we should extend p1 or p2 to outside the AV to avoid blind spots
              TVector3 bp = newBlip.Position;
              float d = BlipUtils::DistToLine(p1,p2,bp);
              if( d > 0 ) {
                // update closest trkdist
                if( newBlip.ProxTrkDist < 0 || d < newBlip.ProxTrkDist ) {
                  newBlip.ProxTrkDist = d;
                  newBlip.ProxTrkID = trk->ID();
                }
                // need to do some math to figure out if this is in
                // the 45 degreee "cone" relative to the start/end 
                if( !newBlip.inCylinder && d < fCylinderRadius ) {
                  float angle1 = asin( d / (p1-bp).Mag() ) * 180./3.14159;
                  float angle2 = asin( d / (p2-bp).Mag() ) * 180./3.14159;
                  if( angle1 < 45. && angle2 < 45. ) newBlip.inCylinder = true;
                }
              }
            }//endloop over trks
           
            if( fApplyTrkCylinderCut && newBlip.inCylinder ) continue;
            
            // ----------------------------------------
            // if we made it this far, the blip is good!
            // associate this blip with the hits and clusters within it
            newBlip.ID = blips.size();
            blips.push_back(newBlip);
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].BlipID = newBlip.ID;
              for( auto& h : hc.HitIDs ) hitinfo[h].blipid = newBlip.ID;
            }

  
          }//endif ncands > 0
        }//endloop over caloplane ("Plane A") clusters
      }//endif calo plane has clusters
    }//endloop over TPCs
    }//endloop over Cryostats

    //std::cout<<"Found "<<hitclust.size()<<" clusters and "<<blips.size()<<" blips\n";
    

    //for(size_t i=0; i<hitlist.size(); i++){
      //if (hitinfo[i].trkid >= 0 ) continue;
      //h_chan_nhits->Fill(fGeom.PlaneWireToChannel(geo::WireID(0,hitinfo[i].tpc,hitinfo[i].plane,hitinfo[i].wire)));
      //int clustid = hitinfo[i].clustid;
      //if( clustid >= 0 ) {
        //if( hitclust[clustid].NWires > 1 ) continue;
        //h_chan_nclusts->Fill(fGeom.PlaneWireToChannel(geo::WireID(0,hitinfo[i].tpc,hitinfo[i].plane,hitinfo[i].wire)));
      //}
      //if( hitinfo[i].ismatch    ) continue;
      //if( hitclust[clustid].NWires > 1 ) continue;
      //h_chan_nclusts->Fill(fGeom.PlaneWireToChannel(hitinfo[i].plane,hitinfo[i].wire));
    //}

    
    //*************************************************************************
    // Loop over the vector of blips and perform calorimetry calculations
    //*************************************************************************
    for(size_t i=0; i<blips.size(); i++){
      auto& blip = blips[i];
      blip.Charge = blip.clusters[fCaloPlane].Charge;
      
      // ***** MICROBOONE ****************************
      // --- YZ uniformity correction ---
      // Correct for charge-collection non-uniformity based on Y/Z position
      // (taken from CalibrationdEdx_module)
      //if( fYZUniformityCorr ) blip.Charge *= tpcCalib_provider.YZdqdxCorrection(fCaloPlane,blip.Position.Y(),blip.Position.Z());
      // *********************************************

      // ================================================================================
      // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
      // like associating blip with some nearby track/shower and using its tagged T0)
      //    Method 1: Assume a dE/dx for electrons, use that + local E-field to get recomb.
      //    Method 2: ESTAR lookup table method ala ArgoNeuT
      // ================================================================================
      float depEl   = std::max(0.0,(double)blip.Charge);
      float Efield  = kNominalEfield;

      // --- Lifetime correction ---
      // Ddisabled by default. Without knowing real T0 of a blip, attempting to 
      // apply this correction can do more harm than good! Note lifetime is in
      // units of 'ms', not microseconds, hence the 1E-3 conversion factor.
      if( fLifetimeCorr && blip.Time>0 ) depEl *= exp( 1e-3*blip.Time/detProp.ElectronLifetime());
      
      // --- SCE corrections ---
      geo::Point_t point( blip.Position.X(),blip.Position.Y(),blip.Position.Z() );
      if( fSCECorr ) {

        // 1) Spatial correction
        if( SCE_provider->EnableCalSpatialSCE() ) {
          // TODO: Deal with cases where X falls outside AV (diffuse out-of-time signal)
          //       For example, maybe re-assign to center of drift volume?
          geo::Vector_t loc_offset = SCE_provider->GetCalPosOffsets(point,blip.TPC);
          point.SetXYZ(point.X()-loc_offset.X(),point.Y()+loc_offset.Y(),point.Z()+loc_offset.Z());
        }
      
        // 2) E-field correction
        //
        // notes:
        //   - GetEfieldOffsets(xyz) and GetCalEfieldOffsets(xyz) return the exact
        //     same underlying E-field offset map; the only difference is the former
        //     is used in the simulation, and the latter in reconstruction (??).
        //   - The SpaceCharge service must have 'EnableCorSCE' and 'EnableCalEfieldSCE'
        //     enabled in order to use GetCalEfieldOffsets
        //   - Blips may appear to be outside the active volume if T0-corrections aren't 
        //     applied to the reconstructed 'X'. SCE map should return (0,0,0) in this case.
        if( SCE_provider->EnableCalEfieldSCE() ) {
          auto const field_offset = SCE_provider->GetCalEfieldOffsets(point,blip.TPC); 
          Efield = detProp.Efield()*std::hypot(1+field_offset.X(),field_offset.Y(),field_offset.Z());;
        }

      }
      
      // METHOD 1
      float recomb  = ModBoxRecomb(fCalodEdx,Efield);
      blip.Energy   = depEl * (1./recomb) * kWion;
      
      // ================================================
      // Save the true blip into the object;
      // at least one cluster needs to match.
      // ================================================
      std::set<int> set_edepids;
      for(auto& hc : blip.clusters ) {
        if( !hc.isValid )   continue; 
        if( hc.EdepID < 0 ) continue;
        set_edepids.insert( hc.EdepID );
      }
      if( set_edepids.size() == 1 ) 
        blip.truth = trueblips[*set_edepids.begin()];
    
    }//endloop over blip vector

  }//End main blip reco function
 
  
  
  //###########################################################
  float BlipRecoAlg::ModBoxRecomb(float dEdx, float Efield) {
    float rho = kLArDensity;
    float Xi = fModBoxB * dEdx / ( Efield * rho );
    return log(fModBoxA+Xi)/Xi;
  }

  float BlipRecoAlg::dQdx_to_dEdx(float dQdx_e, float Efield){
    float rho = kLArDensity;
    float beta  = fModBoxB / (rho * Efield);
    float alpha = fModBoxA;
    return ( exp( beta * kWion * dQdx_e ) - alpha ) / beta;
  }
  
  float BlipRecoAlg::Q_to_E(float Q, float Efield){
    if( Efield != kNominalEfield )  return kWion * (Q/ModBoxRecomb(fCalodEdx,Efield));
    else                            return kWion * (Q/kNominalRecombFactor);
  }
 
  //###########################################################
  void BlipRecoAlg::PrintConfig() {
  
    printf("BlipRecoAlg Configurations\n\n");
    printf("  Input hit collection      : %s\n",          fHitProducer.c_str());
    printf("  Input trk collection      : %s\n",          fTrkProducer.c_str());
    printf("  Max wires per cluster     : %i\n",          fMaxWiresInCluster);
    printf("  Max cluster timespan      : %.1f ticks\n",    fMaxClusterSpan);
    printf("  Min cluster overlap       : %4.1f\n",       fMatchMinOverlap);
    printf("  Clust match sigma-factor  : %4.1f\n",       fMatchSigmaFact);
    printf("  Clust match max dT        : %4.1f ticks\n", fMatchMaxTicks);
    printf("  Charge diff limit         : %.1fe3\n",fMatchQDiffLimit/1000);
    printf("  Charge ratio min          : %.1f\n",fMatchMinQRatio);      
    for(int i=0; i<kNplanes; i++){
    if( i == fCaloPlane ) continue;
    printf("  pl%i matches per cand      : %4.2f\n",       i,h_nmatches[i]->GetMean());
    }
    printf("\n");
    
  }
  
}
