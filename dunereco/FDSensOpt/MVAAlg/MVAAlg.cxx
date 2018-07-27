/////////////////////////////////////////////////////////////////
//  \file MVAAlg.cxx
//  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////

#include "dune/FDSensOpt/MVAAlg/MVAAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::utils::TrackPitchInView() 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TF1.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

//--------------------------------------------------------------------------------
dunemva::MVAAlg::MVAAlg( fhicl::ParameterSet const& p )
  : fReader("")
    , fCalorimetryAlg (p.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fMakeAnaTree    (p.get<bool>("MakeAnaTree"))
    , fMakeWeightTree (p.get<bool>("MakeWeightTree"))
{

  this->reconfigure(p);
  if(fMakeAnaTree) this->MakeTree();

  fFile.open("events.txt"); // for displaying selected events

  if(!fMakeWeightTree){

    // branch, name, unit, type
    fReader.AddVariable("evtcharge", &evtcharge);
    fReader.AddVariable("ntrack", &ntrack);
    fReader.AddVariable("maxtrklength", &maxtrklength);
    fReader.AddVariable("avgtrklength", &avgtrklength);
    fReader.AddVariable("trkdedx", &trkdedx);
    fReader.AddVariable("trkrch", &trkrch);
    fReader.AddVariable("trkrt", &trkrt);
    fReader.AddVariable("trkfr", &trkfr);
    fReader.AddVariable("trkpida", &trkpida_save);
    fReader.AddVariable("fract_5_wires", &fract_5_wires);
    fReader.AddVariable("fract_10_wires", &fract_10_wires);
    fReader.AddVariable("fract_50_wires", &fract_50_wires);
    fReader.AddVariable("fract_100_wires", &fract_100_wires);
    fReader.AddVariable("trkcosx", &trkcosx);
    fReader.AddVariable("trkcosy", &trkcosy);
    fReader.AddVariable("trkcosz", &trkcosz);
    fReader.AddVariable("ET", &ET);

    // Nice to plot and verify sig/back sample composition
    //fReader.AddVariable("ccnc", &ccnc);
    //fReader.AddVariable("NuPdg", &NuPdg); // must be floats

    if(fSelect=="nue"){
      fReader.AddVariable("nshower", &nshower);
      fReader.AddVariable("showerdedx", &showerdedx);
      fReader.AddVariable("eshower", &eshower);
      fReader.AddVariable("frshower", &frshower);
      fReader.AddVariable("nhitspershw", &nhitspershw);
      fReader.AddVariable("shwlength", &shwlength);
      fReader.AddVariable("shwmax", &shwmax);
      fReader.AddVariable("shwdisx", &shwdisx);
      fReader.AddVariable("shwdisy", &shwdisy);
      fReader.AddVariable("shwdisz", &shwdisz);
      fReader.AddVariable("shwcosx", &shwcosx);
      fReader.AddVariable("shwcosy", &shwcosy);
      fReader.AddVariable("shwcosz", &shwcosz);
    }


    /*
       fReader.AddVariable("evtcharge", &evtcharge);
       fReader.AddVariable("ntrack", &ntrack);
       fReader.AddVariable("maxtrklength", &maxtrklength);
       fReader.AddVariable("trkdedx", &trkdedx);
       fReader.AddVariable("trkrch", &trkrch);
       fReader.AddVariable("trkrt", &trkrt);
       fReader.AddVariable("trkfr", &trkfr);
       fReader.AddVariable("trkpida", &trkpida_save);
       fReader.AddVariable("fract_5_wires", &fract_5_wires);
       fReader.AddVariable("fract_10_wires", &fract_10_wires);
       fReader.AddVariable("fract_50_wires", &fract_50_wires);
       fReader.AddVariable("fract_100_wires", &fract_100_wires);
       fReader.AddVariable("trkcosx", &trkcosx);
       fReader.AddVariable("trkcosy", &trkcosy);
       fReader.AddVariable("trkcosz", &trkcosz);
       fReader.AddVariable("ET", &ET);
       fReader.AddVariable("nshower", &nshower);
       fReader.AddVariable("showerdedx", &showerdedx);
       fReader.AddVariable("eshower", &eshower);
       fReader.AddVariable("frshower", &frshower);
       fReader.AddVariable("nhitspershw", &nhitspershw);
       fReader.AddVariable("shwlength", &shwlength);
       fReader.AddVariable("shwmax", &shwmax);
       fReader.AddVariable("shwdisx", &shwdisx);
       fReader.AddVariable("shwdisy", &shwdisy);
       fReader.AddVariable("shwdisz", &shwdisz);
       fReader.AddVariable("shwcosx", &shwcosx);
       fReader.AddVariable("shwcosy", &shwcosy);
       fReader.AddVariable("shwcosz", &shwcosz);
       */


    fMVAMethod =p.get< std::string >("MVAMethod");

    std::string weightFileFull;
    fWeightFile=p.get< std::string >("WeightFile");
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fWeightFile, weightFileFull);
    fReader.BookMVA(fMVAMethod, weightFileFull);

  }

  fWeightTree[0] = tfs->make<TTree>("all","a Tree for MVA");
  fWeightTree[1] = tfs->make<TTree>("sig","a Tree for MVA");
  fWeightTree[2] = tfs->make<TTree>("nc","a Tree for MVA");
  fWeightTree[3] = tfs->make<TTree>("cc","a Tree for MVA");
  for (int i = 0; i<4; ++i){
    fWeightTree[i]->Branch("run",&run,"run/I");
    fWeightTree[i]->Branch("subrun",&subrun,"subrun/I");
    fWeightTree[i]->Branch("event",&event,"event/I");
    fWeightTree[i]->Branch("itype",&itype,"itype/I");

    // all tmva reader variables must be floats
    fWeightTree[i]->Branch("oscpro",&fOscPro,"oscpro/F");
    fWeightTree[i]->Branch("weight",&weight,"weight/F");
    fWeightTree[i]->Branch("evtcharge",&evtcharge,"evtcharge/F");
    fWeightTree[i]->Branch("ntrack",&ntrack,"ntrack/F");
    fWeightTree[i]->Branch("maxtrklength",&maxtrklength,"maxtrklength/F");
    fWeightTree[i]->Branch("avgtrklength",&avgtrklength,"avgtrklength/F");
    fWeightTree[i]->Branch("trkdedx",&trkdedx,"trkdedx/F");
    fWeightTree[i]->Branch("trkrch",&trkrch,"trkrch/F");
    fWeightTree[i]->Branch("trkrt",&trkrt,"trkrt/F");
    fWeightTree[i]->Branch("trkfr",&trkfr,"trkfr/F");
    fWeightTree[i]->Branch("trkpida",&trkpida_save,"trkpida/F");
    fWeightTree[i]->Branch("nshower",&nshower,"nshower/F");
    fWeightTree[i]->Branch("showerdedx",&showerdedx,"showerdedx/F");
    fWeightTree[i]->Branch("eshower",&eshower,"eshower/F");
    fWeightTree[i]->Branch("frshower",&frshower,"frshower/F");
    fWeightTree[i]->Branch("nhitspershw",&nhitspershw,"nhitspershw/F");
    fWeightTree[i]->Branch("shwlength",&shwlength,"shwlength/F");
    fWeightTree[i]->Branch("shwmax",&shwmax,"shwmax/F");
    fWeightTree[i]->Branch("fract_5_wires",&fract_5_wires,"fract_5_wires/F");
    fWeightTree[i]->Branch("fract_10_wires",&fract_10_wires,"fract_10_wires/F");
    fWeightTree[i]->Branch("fract_50_wires",&fract_50_wires,"fract_50_wires/F");
    fWeightTree[i]->Branch("fract_100_wires",&fract_100_wires,"fract_100_wires/F");
    fWeightTree[i]->Branch("shwdis",&shwdis,"shwdis/F");
    fWeightTree[i]->Branch("shwdisx",&shwdisx,"shwdisx/F");
    fWeightTree[i]->Branch("shwdisy",&shwdisy,"shwdisy/F");
    fWeightTree[i]->Branch("shwdisz",&shwdisz,"shwdisz/F");
    fWeightTree[i]->Branch("shwcosx",&shwcosx,"shwcosx/F");
    fWeightTree[i]->Branch("shwcosy",&shwcosy,"shwcosy/F");
    fWeightTree[i]->Branch("shwcosz",&shwcosz,"shwcosz/F");
    fWeightTree[i]->Branch("trkcosx",&trkcosx,"trkcosx/F");
    fWeightTree[i]->Branch("trkcosy",&trkcosy,"trkcosy/F");
    fWeightTree[i]->Branch("trkcosz",&trkcosz,"trkcosz/F");
    fWeightTree[i]->Branch("ET",&ET,"ET/F");

    // Nice to double check samples through to the end
    //fWeightTree[i]->Branch("ccnc",&ccnc,"ccnc/F");
    //fWeightTree[i]->Branch("NuPdg",&NuPdg,"NuPdg/F");
  }

  int mva_bins = 80;
  mva_nue_osc  = tfs->make<TH1D>("mva_nue_osc","mva_nue_osc",mva_bins,-1,1);
  mva_nc       = tfs->make<TH1D>("mva_nc","mva_nc",mva_bins,-1,1);
  mva_numu     = tfs->make<TH1D>("mva_numu","mva_numu",mva_bins,-1,1);
  mva_nue_beam = tfs->make<TH1D>("mva_nue_beam","mva_nue_beam",mva_bins,-1,1);
  mva_nutau    = tfs->make<TH1D>("mva_nutau","mva_nutau",mva_bins,-1,1);

  mva_numu_nue    = tfs->make<TH1D>("mva_numu_nue","mva_numu_nue",mva_bins,-1,1);
  mva_nue_nue    = tfs->make<TH1D>("mva_nue_nue","mva_nue_nue",mva_bins,-1,1);
  mva_numu_numu    = tfs->make<TH1D>("mva_numu_numu","mva_numu_numu",mva_bins,-1,1);
  mva_nue_numu    = tfs->make<TH1D>("mva_nue_numu","mva_nue_numu",mva_bins,-1,1);
  mva_numu_nutau    = tfs->make<TH1D>("mva_numu_nutau","mva_numu_nutau",mva_bins,-1,1);
  mva_nue_nutau    = tfs->make<TH1D>("mva_nue_nutau","mva_nue_nutau",mva_bins,-1,1);

  enu_nc       = tfs->make<TH1D>("enu_nc","enu_nc",100,0,100);
  enu_numu_nue    = tfs->make<TH1D>("enu_numu_nue","enu_numu_nue",100,0,100);
  enu_nue_nue    = tfs->make<TH1D>("enu_nue_nue","enu_nue_nue",100,0,100);
  enu_numu_numu    = tfs->make<TH1D>("enu_numu_numu","enu_numu_numu",100,0,100);
  enu_nue_numu    = tfs->make<TH1D>("enu_nue_numu","enu_nue_numu",100,0,100);
  enu_numu_nutau    = tfs->make<TH1D>("enu_numu_nutau","enu_numu_nutau",100,0,100);
  enu_nue_nutau    = tfs->make<TH1D>("enu_nue_nutau","enu_nue_nutau",100,0,100);

  mva_nue_osc->Sumw2();
  mva_nc->Sumw2();
  mva_numu->Sumw2();
  mva_nue_beam->Sumw2();
  mva_nutau->Sumw2();

  mva_numu_nue->Sumw2();
  mva_nue_nue->Sumw2();
  mva_numu_numu->Sumw2();
  mva_nue_numu->Sumw2();
  mva_numu_nutau->Sumw2();
  mva_nue_nutau->Sumw2();

  enu_nc->Sumw2();
  enu_numu_nue->Sumw2();
  enu_nue_nue->Sumw2();
  enu_numu_numu->Sumw2();
  enu_nue_numu->Sumw2();
  enu_numu_nutau->Sumw2();
  enu_nue_nutau->Sumw2();


  char name[5][100] = {"#nu_{e}^{osc}","NC","#nu_{#mu} CC","#nu_{e}^{beam}","#nu_{#tau} CC"};

  for (int i = 0; i<nSamples; ++i){
    events_truth[i] = 0;
    events_reco[i] = 0;
    enu[i] = tfs->make<TH1D>(Form("enu_%d",i),name[i],80,0,20);
    enu[i]->Sumw2();
    enu_osc[i] = tfs->make<TH1D>(Form("enu_osc_%d",i),name[i],160,0,20);
    enu_osc[i]->Sumw2();
    hvtxx[i] = tfs->make<TH1D>(Form("hvtxx_%d",i), name[i],100,-10,10);
    hvtxy[i] = tfs->make<TH1D>(Form("hvtxy_%d",i), name[i],100,-10,10);
    hvtxz[i] = tfs->make<TH1D>(Form("hvtxz_%d",i), name[i],100,-10,10);
    hvtxx[i]->Sumw2();
    hvtxy[i]->Sumw2();
    hvtxz[i]->Sumw2();
    htrklen[i] = tfs->make<TH1D>(Form("htrklen_%d",i), name[i],100,0,600);
    htrklen[i]->Sumw2();
    htrkdedx[i] = tfs->make<TH1D>(Form("htrkdedx_%d",i), name[i], 100,0,20);
    htrkdedx[i]->Sumw2();
    hrch[i] = tfs->make<TH1D>(Form("hrch_%d",i),name[i],100,0,1);
    hrch[i]->Sumw2();
    hrt[i] = tfs->make<TH1D>(Form("hrt_%d",i),name[i],100,0,1);
    hrt[i]->Sumw2();
    hpida[i] = tfs->make<TH1D>(Form("hpida_%d",i),name[i],100,0,50);
    hpida[i]->Sumw2();
    hnshw[i] = tfs->make<TH1D>(Form("hnshw_%d",i),name[i],20,0,20);
    hnshw[i]->Sumw2();
    heshw[i] = tfs->make<TH1D>(Form("heshw_%d",i),name[i],100,0,10);
    heshw[i]->Sumw2();
    hshwdedx[i] = tfs->make<TH1D>(Form("hshwdedx_%d",i),name[i],100,0,20);
    hshwdedx[i]->Sumw2();
    hdisx[i] = tfs->make<TH1D>(Form("hdisx_%d",i),name[i],100,-20,20);
    hdisx[i]->Sumw2();
    hdisy[i] = tfs->make<TH1D>(Form("hdisy_%d",i),name[i],100,-20,20);
    hdisy[i]->Sumw2();
    hdisz[i] = tfs->make<TH1D>(Form("hdisz_%d",i),name[i],100,-20,20);
    hdisz[i]->Sumw2();
    hfrshower[i] = tfs->make<TH1D>(Form("hfrshower_%d",i),name[i],110,0,1.1);
    hfrshower[i]->Sumw2();
    hnhitspershw[i] = tfs->make<TH1D>(Form("hnhitspershw_%d",i),name[i],100,0,20);
    hnhitspershw[i]->Sumw2();
    hevtcharge[i] = tfs->make<TH1D>(Form("hevtcharge_%d",i),name[i],100,0,3e6);
    hevtcharge[i]->Sumw2();
    hfr100w[i] = tfs->make<TH1D>(Form("hfr100w_%d",i),name[i],110,0,1.1);
    hfr100w[i]->Sumw2();
  }

  fPOT = tfs->make<TTree>("pottree","pot tree");
  fPOT->Branch("pot",&pot,"pot/D");
  fPOT->Branch("run",&run,"run/I");
  fPOT->Branch("subrun",&subrun,"subrun/I");

}

dunemva::MVAAlg::~MVAAlg(){
  fFile.close();
}

//--------------------------------------------------------------------------------
void dunemva::MVAAlg::reconfigure(fhicl::ParameterSet const& p){

  fRawDigitModuleLabel    =   p.get< std::string >("RawDigitModuleLabel");
  fWireModuleLabel        =   p.get< std::string >("WireModuleLabel");
  fHitsModuleLabel        =   p.get< std::string >("HitsModuleLabel");
  fTrackModuleLabel       =   p.get< std::string >("TrackModuleLabel");
  fShowerModuleLabel      =   p.get< std::string >("ShowerModuleLabel");
  fClusterModuleLabel     =   p.get< std::string >("ClusterModuleLabel");
  fVertexModuleLabel      =   p.get< std::string >("VertexModuleLabel");
  fGenieGenModuleLabel    =   p.get< std::string >("GenieGenModuleLabel");
  fPOTModuleLabel         =   p.get< std::string >("POTModuleLabel"); 
  fFlashModuleLabel       =   p.get< std::string >("FlashModuleLabel");
  fCalorimetryModuleLabel =   p.get< std::string >("CalorimetryModuleLabel");
  fSelect                 =   p.get< std::string >("Select");
  fBeamMode               =   p.get< std::string >("BeamMode","FHC");
  fFidVolCut              =   p.get< double      >("FidVolCut");
}

//--------------------------------------------------------------------------------
void dunemva::MVAAlg::Run( const art::Event & evt, double& result, double& wgt ){

  this->ResetVars();
  this->PrepareEvent(evt); // does not reset vars, make sure to reset after evaluating
  this->CalculateInputs();
  if(!fMakeWeightTree){

    if (isinfidvol){
      result = fReader.EvaluateMVA(fMVAMethod);
      mf::LogVerbatim("MVASelect") << fMVAMethod
        << " returned " << result;

      // Fill histogram of MVA value for each type

      // itype is... (see instantiation of "name")
      // 0 for oscillated NuE CC
      // 1 for NC background
      // 2 for NuMu CC (oscillated component negligible)
      // 3 for beam NuE
      // 4 for NuTau CC

      if (itype==0) mva_nue_osc->Fill(result,1.);
      if (itype==1) mva_nc->Fill(result,1.);
      if (itype==2) mva_numu->Fill(result,1.);
      if (itype==3) mva_nue_beam->Fill(result,1.);
      if (itype==4) mva_nutau->Fill(result,1.);

      if (ccnc_truth==1){
      }
      else if (std::abs(pntype_flux)==14&&std::abs(nuPDG_truth)==12){
        mva_numu_nue->Fill(result, oscpro);
      }
      else if (std::abs(pntype_flux)==12&&std::abs(nuPDG_truth)==12){
        mva_nue_nue->Fill(result, oscpro);
      }
      else if (std::abs(pntype_flux)==14&&std::abs(nuPDG_truth)==14){
        mva_numu_numu->Fill(result, oscpro);
      }
      else if (std::abs(pntype_flux)==12&&std::abs(nuPDG_truth)==14){
        mva_nue_numu->Fill(result, oscpro);
      }
      else if (std::abs(pntype_flux)==14&&std::abs(nuPDG_truth)==16){
        mva_numu_nutau->Fill(result, oscpro);
      }
      else if (std::abs(pntype_flux)==12&&std::abs(nuPDG_truth)==16){
        mva_nue_nutau->Fill(result, oscpro);
      }
    }
    else{//not in fiducial volume
      result = -2.;
    }
    if (isinfidvoltruth){
      if (ccnc_truth==1){
        enu_nc->Fill(enu_truth, oscpro);
      }
      else if (std::abs(pntype_flux)==14&&std::abs(nuPDG_truth)==12){
        enu_numu_nue->Fill(enu_truth, oscpro);
      }
      else if (std::abs(pntype_flux)==12&&std::abs(nuPDG_truth)==12){
        enu_nue_nue->Fill(enu_truth, oscpro);
      }
      else if (std::abs(pntype_flux)==14&&std::abs(nuPDG_truth)==14){
        enu_numu_numu->Fill(enu_truth, oscpro);
      }
      else if (std::abs(pntype_flux)==12&&std::abs(nuPDG_truth)==14){
        enu_nue_numu->Fill(enu_truth, oscpro);
      }
      else if (std::abs(pntype_flux)==14&&std::abs(nuPDG_truth)==16){
        enu_numu_nutau->Fill(enu_truth, oscpro);
      }
      else if (std::abs(pntype_flux)==12&&std::abs(nuPDG_truth)==16){
        enu_nue_nutau->Fill(enu_truth, oscpro);
      }
    }
  }

  wgt = weight;
}


//--------------------------------------------------------------------------------
void dunemva::MVAAlg::endSubRun(const art::SubRun& sr){
  art::Handle< sumdata::POTSummary > potListHandle;

  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    pot = potListHandle->totpot;
  else
    pot = 0.;
  if (fPOT) fPOT->Fill();
}


//--------------------------------------------------------------------------------
void dunemva::MVAAlg::CalculateInputs( ){

  itype = -99999;
  fOscPro = -99999;
  weight = -99999;
  evtcharge = -99999;
  ntrack = -99999;
  maxtrklength = -99999;
  avgtrklength = -99999;
  trkdedx = -99999;
  trkrch = -99999;
  trkrt = -99999;
  trkfr = -99999;
  trkpida_save = -99999;
  nshower = -99999;
  showerdedx = -99999;
  eshower = -99999;
  frshower = -99999;
  nhitspershw = -99999;
  shwlength = -99999;
  shwmax = -99999;
  fract_5_wires = -99999;
  fract_10_wires = -99999;
  fract_50_wires = -99999;
  fract_100_wires = -99999;
  shwdis = -99999;
  ET = 0;

  // itype is... (see instantiation of "name")
  // 0 for NuE signal
  // 1 for NC background
  // 2 for NuMu CC
  // 3 for beam NuE
  // 4 for NuTau CC

  // ccnc  = ccnc_truth;
  // pntype_flux is the original type
  // nuPDG_truth is the type at the FD
  //NuPdg = std::abs(nuPDG_truth); 

  itype = -1;
  if (ccnc_truth==1){
    itype = 1; // NC background
  }
  else if (std::abs(pntype_flux)!=12&&std::abs(nuPDG_truth)==12){
    // oscillated NuE (from NuMu)
    itype = 0;
  }
  else if (std::abs(nuPDG_truth)==12){
    // intrinsic NuE
    itype = 3;
  }
  else if (std::abs(nuPDG_truth)==14){
    // intrinsic NuMu (oscillated from NuE is negligible)
    itype = 2;
  }
  else if (std::abs(nuPDG_truth)==16){
    // NuTau background
    itype = 4; 
  }
  if (itype == -1){
    //std::cout << "Unknown type: "
    //	  << pntype_flux << " to " << nuPDG_truth
    //	  << " (ccnc=" << ccnc_truth << ")"
    //	  << std::endl;
  }      

  oscpro = this->OscPro(ccnc_truth,pntype_flux,nuPDG_truth,enu_truth);
  float norm = this->Norm(ccnc_truth,pntype_flux,nuPDG_truth, subrun);
  weight = oscpro*norm;
  fOscPro = oscpro;


  // Choose the most upstream vertex, even if it is just 
  //   the beginning of a single track
  bool findvtx = false;
  float vtxx = 0;
  float vtxy = 0;
  float vtxz = 100000;
  if (nvtx>0){
    for (int i = 0; i<nvtx&&i<kMaxVertices; ++i){
      if (vtx[i][2]<vtxz){
        vtxx = vtx[i][0];
        vtxy = vtx[i][1];
        vtxz = vtx[i][2];
        findvtx = true;
      }
    }
  }
  else if (ntracks_reco>0){
    for (int i = 0; i<ntracks_reco; ++i){
      if (trkstartz[i]<vtxz){
        vtxx = trkstartx[i];
        vtxy = trkstarty[i];
        vtxz = trkstartz[i];
        findvtx = true;
      }
    }
  }

  if (std::abs(nuvtxx_truth)<360-50&&
      std::abs(nuvtxy_truth)<600-50&&
      nuvtxz_truth>50&&nuvtxz_truth<1394-150){
    isinfidvoltruth = 1;
    //if (enu_truth>20) continue;
    //if (enu_truth<0.5) continue;
    //if (ccnc_truth==1&&Y_truth*enu_truth<0.5) continue;
    events_truth[itype]+=oscpro*norm;
    enu[itype]->Fill(enu_truth,norm);
    enu_osc[itype]->Fill(enu_truth,oscpro*norm);
    if (findvtx){
      hvtxx[itype]->Fill(vtxx-nuvtxx_truth, oscpro*norm);
      hvtxy[itype]->Fill(vtxy-nuvtxy_truth, oscpro*norm);
      hvtxz[itype]->Fill(vtxz-nuvtxz_truth, oscpro*norm);
    }
  }

  // get avg trk length
  avgtrklength = 0.;
  for (int i = 0; i<ntracks_reco; ++i){
    avgtrklength += trklen[i];
  }
  if(ntracks_reco)
    avgtrklength /= ntracks_reco;
  ntrack = 1.*ntracks_reco;

  // Try to approximate transverse energy
  unsigned int ntracks_notvtx(0);
  if(findvtx && ntracks_reco>0){	

    double Edir[3] = {0.,0.,0.};
    double Etot(0.);
    for (int i = 0; i<nshws; ++i){

      // don't count particles not out of event vertex
      // VERY rough approximation here..
      if( std::sqrt(   std::pow(vtxx-shwstartx[i],2)
            + std::pow(vtxy-shwstarty[i],2)
            + std::pow(vtxz-shwstartz[i],2) > 2. ) ) continue;

      Edir[0] += std::abs(shwenergy[i][shwbestplane[i]])*shwdcosx[i];
      Edir[1] += std::abs(shwenergy[i][shwbestplane[i]])*shwdcosy[i];
      Edir[2] += std::abs(shwenergy[i][shwbestplane[i]])*shwdcosz[i];
      Etot    += std::abs(shwenergy[i][shwbestplane[i]]);
    }
    for (int i = 0; i<ntracks_reco; ++i){

      // don't count particles not out of event vertex
      // VERY rough approximation here..
      if( std::sqrt(   std::pow(vtxx-trkstartx[i],2)
            + std::pow(vtxy-trkstarty[i],2)
            + std::pow(vtxz-trkstartz[i],2) > 2. ) ){
        //if(trklen[i]<140) ntracks_notvtx++; this made it worse
        ntracks_notvtx++;
        continue;
      }

      Edir[0] += std::abs(trkke[i][trkbestplane[i]])*trkstartdcosx[i];
      Edir[1] += std::abs(trkke[i][trkbestplane[i]])*trkstartdcosy[i];
      Edir[2] += std::abs(trkke[i][trkbestplane[i]])*trkstartdcosz[i];
      Etot    += std::abs(trkke[i][trkbestplane[i]]);
    }

    if(Etot!=0.){
      Edir[0] /= Etot;
      Edir[1] /= Etot;
      Edir[2] /= Etot;
    }	

    ET = std::sqrt( Edir[0]*Edir[0] + Edir[1]*Edir[1] );
    //ET /= Etot;
    if(ET>1||ET<0){
      std::cout << ET << ", " << Etot << std::endl;
    }

    if( ET > 0.9 ){
      fFile << run << " " << subrun << " " << event << "\n";
    }

  } else {
    //overflow
    ET = 0;
  }

  // making loose assumptions about origin and dimensions...
  if( findvtx&&
      std::abs(vtxx) < 360-50 &&
      std::abs(vtxy) < 600-50 &&
      vtxz > 50 &&
      vtxz < 1394-150 ){//in the fiducial volume

    isinfidvol = 1;
    events_reco[itype] += oscpro*norm;
    //ntrack = ntracks_reco;
    //ntrack = ntracks_notvtx;
    maxtrklength = 0;
    int itrk = -1;
    for (int i = 0; i<ntracks_reco; ++i){
      if (trklen[i]>maxtrklength){
        maxtrklength = trklen[i];
        itrk = i;
      }
    }
    //if (itype==2&&maxtrklength<5) std::cout<<run<<" "<<subrun<<" "<<event<<std::endl;
    htrklen[itype]->Fill(maxtrklength,oscpro*norm);
    //if (maxtrklength>300) continue;

    //std::vector<std::map<int,float> > vtrkdedx(3); don't use until Reco saves dEds per hit
    std::vector<std::vector<std::map<int,float> > > vtrktick(fGeom->NTPC());
    for (size_t i = 0; i<fGeom->NTPC(); ++i){
      vtrktick[i].resize(3);
    }

    trkdedx = 0;
    trkrch = 0;
    trkrt = 0;
    trkfr = 0;
    trkpida_save = 0;
    evtcharge = 0;

    nshower = 1.*nshws;
    hnshw[itype]->Fill(nshws,oscpro*norm);

    eshower = 0;
    showerdedx = 0;
    int ishw = -1;
    for (int i = 0; i<nshws; ++i){
      if (shwbestplane[i]>=0&&shwbestplane[i]<3){
        if (shwenergy[i][shwbestplane[i]]>eshower){
          eshower = shwenergy[i][shwbestplane[i]];
          showerdedx = shwdedx[i][shwbestplane[i]];
          ishw = i;
        }
      }
    }
    shwdis = 1000;
    shwdisx = 1000;
    shwdisy = 1000;
    shwdisz = 1000;
    shwcosx = -2;
    shwcosy = -2;
    shwcosz = -2;
    if (ishw!=-1){
      shwdis = sqrt(pow(shwstartx[ishw]-vtxx,2)+
          pow(shwstarty[ishw]-vtxy,2)+
          pow(shwstartz[ishw]-vtxz,2));
      shwdisx = shwstartx[ishw]-vtxx;
      shwdisy = shwstarty[ishw]-vtxy;
      shwdisz = shwstartz[ishw]-vtxz;
      shwcosx = shwdcosx[ishw];
      shwcosy = shwdcosy[ishw];
      shwcosz = shwdcosz[ishw];
    }
    trkcosx = -2;
    trkcosy = -2;
    trkcosz = -2;
    if (itrk!=-1){
      trkcosx = trkstartdcosx[itrk];
      trkcosy = trkstartdcosy[itrk];
      trkcosz = trkstartdcosz[itrk];
    }
    int shwwire0 = 100000;
    int shwwire1 = -1;
    int offset = 0;
    //int lastwire = -1;
    for (int i = 0; i<nhits && i<40000; ++i){
      if (itrk!=-1){
        if (hit_trkkey[i]==itrk){
          // dont use hit_dEds until is is saved in Reco
          // for now it is alwasy -9999
          //vtrkdedx[hit_plane[i]][hit_wire[i]] = hit_dEds[i];
          vtrktick[hit_tpc[i]][hit_plane[i]][hit_wire[i]] = hit_peakT[i];
        }
      }
      if (hit_plane[i]==2){
        evtcharge += hit_charge[i]*exp(hit_peakT[i]*0.5/taulife);
        //evtcharge += hit_charge[i];
        //	    if (run==20000000&&event==4137){
        //	      std::cout<<hit_wire[i]<<" "<<shwwire0<<" "<<shwwire1<<std::endl;
        //	    }
        if (hit_shwkey[i]==ishw){
          //if (lastwire == -1) lastwire = hit_wire[i];
          //if (hit_wire[i]<lastwire) offset = 479;
          offset = (hit_tpc[i]/4)*fGeom->Nwires(2);
          if (hit_wire[i]+offset<shwwire0){
            shwwire0 = hit_wire[i]+offset;
          }
          if (hit_wire[i]+offset>shwwire1){
            shwwire1 = hit_wire[i]+offset;
          }
          //lastwire = hit_wire[i];
          //	      if (hit_trkkey[i]>=0){
          //		if (sqrt(pow(trkstartx[hit_trkkey[i]]-vtxx,2)+
          //			 pow(trkstarty[hit_trkkey[i]]-vtxy,2)+
          //			 pow(trkstartz[hit_trkkey[i]]-vtxz,2))<shwdis)
          //		  shwdis = sqrt(pow(trkstartx[hit_trkkey[i]]-vtxx,2)+
          //				pow(trkstarty[hit_trkkey[i]]-vtxy,2)+
          //				pow(trkstartz[hit_trkkey[i]]-vtxz,2));
          //		if (sqrt(pow(trkendx[hit_trkkey[i]]-vtxx,2)+
          //			 pow(trkendy[hit_trkkey[i]]-vtxy,2)+
          //			 pow(trkendz[hit_trkkey[i]]-vtxz,2))<shwdis)
          //		  shwdis = sqrt(pow(trkendx[hit_trkkey[i]]-vtxx,2)+
          //				pow(trkendy[hit_trkkey[i]]-vtxy,2)+
          //				pow(trkendz[hit_trkkey[i]]-vtxz,2));
          //	      }
        }
      }
    }
    hevtcharge[itype]->Fill(evtcharge,oscpro*norm);
    std::vector<float> shwph;
    if (shwwire1>=shwwire0){
      shwph.resize(shwwire1-shwwire0+1,0);
      shwlength = shwwire1-shwwire0+1;
    }
    else{
      shwlength = 0;
    }
    //	if (run==20000001&&event==4394){
    //	  std::cout<<ishw<<" "<<shwwire0<<" "<<shwwire1<<std::endl;
    //	}

    if (itrk!=-1){ // itrk is track index of longest track
      int ipl[2] = {-1,-1}; // top 2 with most hits
      int maxhits[2] = {0,0};

      // find the two planes with the most hits
      for (int i = 0; i<3; ++i){
        // set maxhits[0]
        int tothits = 0;
        for (size_t j = 0; j<fGeom->NTPC(); ++j){
          tothits += vtrktick[j][i].size();
        }
        if (tothits>maxhits[0]){ // vtrktick[i] is a map from wire number to charge 
          if (ipl[1]!=-1){ // dont do first time (i=0)
            maxhits[1] = maxhits[0];
            ipl[1] = ipl[0];
          }
          maxhits[0] = tothits;
          ipl[0] = i;
        }
        // set maxhits[1], isn't called first time (i=0)
        else if (tothits>maxhits[1]){
          maxhits[1] = tothits;
          ipl[1] = i;
        }
      }

      //int totwires = 0;
      //	  int startw[2] = {-1,-1};
      //	  for (int i = 0; i<2; ++i){
      //	    if (ipl[i]==-1) continue;
      //	    int startw[i] = int(vtrkdedx[ipl[i]].begin()->first+0.3*(vtrkdedx[ipl[i]].rbegin()->first-vtrkdedx[ipl[i]].begin()->first));
      //	  }

      /* can't calculate using vtrkdedx yet since PANDORA didn't save hit_dEds for each hit
         for (int i = 0; i<2; ++i){
         if (ipl[i]==-1) continue;
         for (auto const &j : vtrkdedx[ipl[i]]){
         trkdedx += j.second;
         ++totwires;
         }
         }
         if (totwires) trkdedx/=(totwires);
         */

      if(itrk==-1) trkdedx = 0;
      else trkdedx = trkke[itrk][trkbestplane[itrk]]/trklen[itrk];

      float qall = -1;
      float qtrk = -1;
      float pida = -1;
      float ipl_evtcharge = 0;
      int npida = 0;
      std::vector<float> vhitq;
      for (int i = 0; i<nhits && i<40000; ++i){ // only 40000 hits saved after Reco
        for (int j = 0; j<2; ++j){  // planes with top 2 # of hits (ipl[j])
          if (hit_plane[i] == ipl[j]){
            if (hit_trkkey[i] == itrk){
              if (hit_resrange[i]<30&&hit_resrange[i]>0.6){
                pida += hit_dEds[i]*pow(hit_resrange[i],0.42);
                ++npida;
              }
            } 
            if (vtrktick[hit_tpc[i]][ipl[j]].find(hit_wire[i])!=vtrktick[hit_tpc[i]][ipl[j]].end()){
              if (std::abs(vtrktick[hit_tpc[i]][ipl[j]][hit_wire[i]]-hit_peakT[i])<200){
                vhitq.push_back(hit_charge[i]);
                qall+=hit_charge[i];
                if (hit_trkkey[i] == itrk){
                  qtrk+=hit_charge[i];
                } // if hit is actually in the longest track
              } // if event-hit is within 200 ticks of track-hit on this wire 
            } // if event-hit wire is in track-hit map (wire-->peak time), track has a hit on this wire

            // sum the denominator to get trk charge fraction of event charge
            // (normal evtcharge sums only on collection plane)
            ipl_evtcharge += hit_charge[i];

          } // if hit is in top 2 planes
        } // planes with top 2 # of hits
      } // all nhits in event

      // sort all of the hit charges for hits within a 200tick width of the track
      std::sort(vhitq.begin(), vhitq.end());
      float q1 = 0;
      float q2 = 0;
      for (size_t i = 0; i<vhitq.size(); ++i){
        if (i<vhitq.size()/2){
          q1+=vhitq[i];
        }
        else{
          q2+=vhitq[i];
        }
      }
      hrch[itype]->Fill(q1/q2,oscpro*norm);
      hrt[itype]->Fill(TMath::Min(float(0.999999),qtrk/qall),oscpro*norm);
      hpida[itype]->Fill(pida/npida,oscpro*norm);
      if(q2!=0) trkrch = q1/q2;
      trkrt = qtrk/qall;
      if(ipl_evtcharge!=0) trkfr = qtrk/ipl_evtcharge;

      if(trkfr>1.){
        std::cout << "Erroneous track charge fraction: " << trkfr
          << " = " << qtrk << " / " << ipl_evtcharge << std::endl;
      }

      //trkpida = pida/npida;
      trkpida_save = trkpida[itrk][trkbestplane[itrk]];
      if (trkpida_save>100.) trkpida_save = 100.;
      if (trkpida_save<0.) trkpida_save = 0.;
    }//itrk!=-1
    else{
      hrch[itype]->Fill(0.,oscpro*norm);
      hrt[itype]->Fill(0.,oscpro*norm);
      hpida[itype]->Fill(0.,oscpro*norm);
    }
    if (trkdedx<0.)   trkdedx = 0.;
    if (trkdedx>100.) trkdedx = 100.;
    htrkdedx[itype]->Fill(trkdedx, oscpro*norm);

    frshower = 0;
    int totalshwhits = 0;
    std::map<int,int> shwwires;
    offset = 0;
    //lastwire = -1;

    for (int i = 0; i<nhits && i<40000; ++i){
      if (hit_plane[i]==2&&hit_shwkey[i]==ishw){
        frshower+=hit_charge[i]*exp(hit_peakT[i]*0.5/taulife);
        shwwires[hit_wire[i]] = 1;
        ++totalshwhits;
        offset = (hit_tpc[i]/4)*fGeom->Nwires(2);
        //if (lastwire ==-1) lastwire = hit_wire[i];
        //if (hit_wire[i]<lastwire) offset = 479;
        shwph[hit_wire[i]+offset-shwwire0] += hit_charge[i]*exp(hit_peakT[i]*0.5/taulife);
        //lastwire = hit_wire[i];
      }
    }
    frshower/=evtcharge;
    if (shwwires.size()){
      nhitspershw = float(totalshwhits)/shwwires.size();
    }
    else 
      nhitspershw = 0;

    shwmax = 0;
    float maxshw = 0;
    fract_5_wires = 0;
    fract_10_wires = 0;
    fract_50_wires = 0;
    fract_100_wires = 0;
    float totshwph = 0;
    for (size_t i = 0; i<shwph.size(); ++i){
      if (shwph[i]>maxshw){
        shwmax = i;
        maxshw = shwph[i];
      }
      totshwph += shwph[i];
      float ph_5_wires = 0;
      float ph_10_wires = 0;
      float ph_50_wires = 0;
      float ph_100_wires = 0;
      for (size_t j = i; j<i+100&&j<shwph.size(); ++j){
        if (j<i+5){
          ph_5_wires += shwph[j];
        }
        if (j<i+10){
          ph_10_wires += shwph[j];
        }
        if (j<i+50){
          ph_50_wires += shwph[j];
        }
        if (j<i+100){
          ph_100_wires += shwph[j];
        }
      }
      if (ph_5_wires>fract_5_wires){
        fract_5_wires = ph_5_wires;
      }
      if (ph_10_wires>fract_10_wires){
        fract_10_wires = ph_10_wires;
      }
      if (ph_50_wires>fract_50_wires){
        fract_50_wires = ph_50_wires;
      }
      if (ph_100_wires>fract_100_wires){
        fract_100_wires = ph_100_wires;
      }
    }
    if (shwlength) shwmax/=shwlength;
    if (totshwph) {
      fract_5_wires/=totshwph;
      fract_10_wires/=totshwph;
      fract_50_wires/=totshwph;
      fract_100_wires/=totshwph;
    }
    //if (frshower<0.1&&itype==0) std::cout<<run<<" "<<subrun<<" "<<event<<std::endl;
    if (ishw!=-1){
      heshw[itype]->Fill(eshower,oscpro*norm);
      if (showerdedx<0.)   showerdedx = 0.;
      if (showerdedx>100.) showerdedx = 100.;
      if (eshower>0.5){
        hshwdedx[itype]->Fill(showerdedx,oscpro*norm);
        hfrshower[itype]->Fill(frshower,oscpro*norm);
        hnhitspershw[itype]->Fill(nhitspershw,oscpro*norm);
        hfr100w[itype]->Fill(fract_100_wires,oscpro*norm);
        hdisx[itype]->Fill(shwstartx[ishw]-trkstartx[itrk],oscpro*norm);
        hdisy[itype]->Fill(shwstarty[ishw]-trkstarty[itrk],oscpro*norm);
        hdisz[itype]->Fill(shwstartz[ishw]-trkstartz[itrk],oscpro*norm);
      }
    }


    //	if (maxtrklength<300&&
    //	    evtcharge>1e5&&
    //	    evtcharge<5e6&&
    //	    nshower>0&&
    //	    eshower>0.5&&
    //	    subrun<666){
    //	  if (itype==0) fWeightTree[1]->Fill();
    //	  if (itype==1) fWeightTree[2]->Fill();
    //	  if (itype==2) fWeightTree[3]->Fill();
    //	  fWeightTree[0]->Fill();
    //	}

    if(fMakeWeightTree){

      // Fill Weight Trees
      //
      // fWeightTree[0] is all
      // fWeightTree[1] is signal
      // fWeightTree[2] is NC background
      // fWeightTree[3] is CC background

      // itype is... (see instantiation of "name")
      // 0 for oscillated NuE CC
      // 1 for NC background
      // 2 for NuMu CC (oscillated component negligible)
      // 3 for beam NuE
      // 4 for NuTau CC

      if(fSelect=="numu"){
        if (itype==2) fWeightTree[1]->Fill(); // signal (numucc)
        if (itype==1) fWeightTree[2]->Fill(); // NC background
        if (itype==0) fWeightTree[3]->Fill(); // CC background (nuecc)
      } 
      else if(fSelect=="nue"){
        if (itype==0) fWeightTree[1]->Fill(); // signal (nuecc)
        if (itype==1) fWeightTree[2]->Fill(); // NC background
        if (itype==2) fWeightTree[3]->Fill(); // CC background (numucc)
      }
      fWeightTree[0]->Fill();

    } // if making weight tree

  } //in fiducial volume
  else {
    mf::LogVerbatim("MVASelect") << "  Not found in fiducial volume. nvtx=" << nvtx << " True vtx = ("
      << nuvtxx_truth <<", "<< nuvtxy_truth <<", "<< nuvtxz_truth << "), Reco vtx = ("<<vtxx<<", "<<vtxy<<", "<<vtxz<<")";
  }

  if(itype == -99999) mf::LogVerbatim("MVASelect") << "  itype not set";
  if(weight == -9999) mf::LogVerbatim("MVASelect") << "  weight not set";
  if(evtcharge == -9999) mf::LogVerbatim("MVASelect") << "  evtcharge not set";
  if(ntrack == -9999) mf::LogVerbatim("MVASelect") << "  ntrack not set";
  if(maxtrklength == -9999) mf::LogVerbatim("MVASelect") << "  maxtrklength not set";
  if(avgtrklength == -9999) mf::LogVerbatim("MVASelect") << "  avgtrklength not set";
  if(trkdedx == -9999) mf::LogVerbatim("MVASelect") << "  trkdedx not set";
  if(trkrch == -9999) mf::LogVerbatim("MVASelect") << "  trkrch not set";
  if(trkrt == -9999) mf::LogVerbatim("MVASelect") << "  trkrt not set";
  if(trkfr == -9999) mf::LogVerbatim("MVASelect") << "  trkfr not set";
  if(trkpida_save == -9999) mf::LogVerbatim("MVASelect") << "  trkpida_save not set";
  if(nshower == -9999) mf::LogVerbatim("MVASelect") << "  nshower not set";
  if(showerdedx == -9999) mf::LogVerbatim("MVASelect") << "  showerdedx not set";
  if(eshower == -9999) mf::LogVerbatim("MVASelect") << "  eshower not set";
  if(frshower == -9999) mf::LogVerbatim("MVASelect") << "  frshower not set";
  if(nhitspershw == -9999) mf::LogVerbatim("MVASelect") << "  nhitspershw not set";
  if(shwlength == -9999) mf::LogVerbatim("MVASelect") << "  shwlength not set";
  if(shwmax == -9999) mf::LogVerbatim("MVASelect") << "  shwmax not set";
  if(fract_5_wires == -9999) mf::LogVerbatim("MVASelect") << "  fract_5_wiresnot set";
  if(fract_10_wires == -9999) mf::LogVerbatim("MVASelect") << "  fract_10_wiresnot set";
  if(fract_50_wires == -9999) mf::LogVerbatim("MVASelect") << "  fract_50_wiresnot set";
  if(fract_100_wires == -9999) mf::LogVerbatim("MVASelect") << "  fract_100_wiresnot set";
  if(shwdis == -9999) mf::LogVerbatim("MVASelect") << "  shwdisnot set";
  //if(ET == 0) mf::LogVerbatim("MVASelect") << "  ET not set";


} // CalculateInputs()


//--------------------------------------------------------------------------------
float dunemva::MVAAlg::Norm(int ccnc, int nu0, int nu1, int subrun){

  float norm = 1;
  float mass = 814.308*1.4e3/1e6; //kt

  float pot_nu = 1.09125e+23;
  float pot_nue = 1.05555e+23;
  float pot_nutau = 2.83966e+23;
  if(fBeamMode == "RHC"){
    pot_nu = 1.89753e+23;
    pot_nue = 1.80693e+23;
    pot_nutau = 4.19518e+23;
  }

  float potpermwyr = 1.47e21/1.07;
  float pot = 150*potpermwyr/mass;//kt*mw*yr

  if (ccnc==1){
    norm = pot/(pot_nu+pot_nue+pot_nutau);
  }
  else if (std::abs(nu0)==14&&std::abs(nu1)==12){
    norm = pot/pot_nue;
  }
  else if (std::abs(nu0)==12&&std::abs(nu1)==12){
    norm = pot/pot_nu;
  }
  else if (std::abs(nu0)==14&&std::abs(nu1)==14){
    norm = pot/pot_nu;
  }
  else if (std::abs(nu0)==12&&std::abs(nu1)==14){
    norm = pot/pot_nutau;
  }
  else if (std::abs(nu0)==14&&std::abs(nu1)==16){
    norm = pot/pot_nutau;
  }
  else if (std::abs(nu0)==12&&std::abs(nu1)==16){
    norm = pot/pot_nue;
  }
  else{
    std::cout << "Unknown oscillation: "
      << nu0 << " to " << nu1 << std::endl;
  }
  return norm;
} // Norm()


//--------------------------------------------------------------------------------
float dunemva::MVAAlg::OscPro(int ccnc, int nu0, int nu1, float NuE){

  float osc_dm2=0.002457;
  float osc_L=1300.; 
  //float sinth23=pow(sin(0.67),2);//sin(3.1415926535/4.0);
  //float sin2th13=0.094;

  //FIXME osceq is created on the heap which was memory leak bound.  For now, a delete has been called at the end of the function.  This needs a more long term solution
  TF1 *osceq = 0;
  if(!osceq){
    osceq = new TF1("f2","sin(1.267*[0]*[1]/x)*sin(1.267*[0]*[1]/x)",0.,120.);
    osceq->SetParameters(osc_dm2,osc_L);
  }

  float OscProb = 0 ;
  float NumuToNutau ;
  float NumuToNue ;
  float NueSurvival ;
  float NumuSurvival ;
  float NueToNutau  ;
  float NueToNumu;

  //float th13 = 0.156;
  //float th23 = 0.670;
  float th13 = 0.148;
  float th23 = 0.738;

  //float th13 = 0.;
  //  float th23 = 0.;

  //NumuToNutau = 4.*sinth23*(1.-sinth23)*pow(1-sin2th13/4,2) ;
  //NumuToNutau *= osceq->Eval(TMath::Abs(NuE)) ;

  //NumuToNue = sinth23*sin2th13*osceq->Eval(TMath::Abs(NuE)) ;
  NumuToNue = pow(sin(th23),2)*pow(sin(2*th13),2)*osceq->Eval(TMath::Abs(NuE)) ;
  //NueSurvival = 1.- sin2th13*osceq->Eval(TMath::Abs(NuE)) ;
  NueSurvival = 1-pow(sin(2*th13),2)*osceq->Eval(TMath::Abs(NuE)) ;
  //NueSurvival = 1.;

  //NumuSurvival = 1. - NumuToNutau - NumuToNue ;
  //NumuSurvival = 1.;
  NumuSurvival = 1-(pow(cos(th13),2)*pow(sin(2*th23),2)+pow(sin(th23),4)*pow(sin(2*th13),2))*osceq->Eval(TMath::Abs(NuE)) ;

  //NumuSurvival = 1.;
  //NueToNutau = (1.-sinth23)*sin2th13*osceq->Eval(TMath::Abs(NuE)) ;

  NumuToNutau = 1 - NumuSurvival - NumuToNue;
  NueToNutau = 1 - NueSurvival - NumuToNue;

  NueToNumu = NumuToNue;  

  if (ccnc==1){
    OscProb = 1;
  }
  else if (std::abs(nu0)==14&&std::abs(nu1)==12){
    OscProb = NumuToNue;
  }
  else if (std::abs(nu0)==12&&std::abs(nu1)==12){
    OscProb = NueSurvival;
  }
  else if (std::abs(nu0)==14&&std::abs(nu1)==14){
    OscProb = NumuSurvival;
  }
  else if (std::abs(nu0)==12&&std::abs(nu1)==14){
    OscProb = NueToNumu;
  }
  else if (std::abs(nu0)==14&&std::abs(nu1)==16){
    OscProb = NumuToNutau;
  }
  else if (std::abs(nu0)==12&&std::abs(nu1)==16){
    OscProb = NueToNutau;
  }
  else{
    std::cout << "Unknown oscillation: "
      << nu0 << " to " << nu1 << std::endl;
  }
  delete osceq;

  return OscProb;
} // OscPro()



//--------------------------------------------------------------------------------
void dunemva::MVAAlg::PrepareEvent(const art::Event& evt){

  //std::cout << " ~~~~~~~~~~~~~~~ MVA: Getting Event Reco ~~~~~~~~~~~~~~ " << std::endl;

  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();
  taulife = detprop->ElectronLifetime();
  isdata = evt.isRealData();

  // * Raw Digits
  art::Handle<std::vector<raw::RawDigit> > rawListHandle;
  std::vector<art::Ptr<raw::RawDigit> > rawlist;
  if (evt.getByLabel(fRawDigitModuleLabel,rawListHandle))
    art::fill_ptr_vector(rawlist, rawListHandle);

  // * wires
  art::Handle< std::vector<recob::Wire>> wireListHandle;
  std::vector<art::Ptr<recob::Wire>> wirelist;
  if (evt.getByLabel(fWireModuleLabel, wireListHandle))
    art::fill_ptr_vector(wirelist, wireListHandle);

  // * hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  // * tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

  // * vertices
  art::Handle< std::vector<recob::Vertex> > vtxListHandle;
  std::vector<art::Ptr<recob::Vertex> > vtxlist;
  if (evt.getByLabel(fVertexModuleLabel,vtxListHandle))
    art::fill_ptr_vector(vtxlist, vtxListHandle);

  // * showers
  art::Handle<std::vector<recob::Shower>> shwListHandle;
  std::vector<art::Ptr<recob::Shower>> shwlist;
  if (evt.getByLabel(fShowerModuleLabel,shwListHandle))
    art::fill_ptr_vector(shwlist, shwListHandle);

  // * flashes
  art::Handle< std::vector<recob::OpFlash> > flashListHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashlist;
  if (evt.getByLabel(fFlashModuleLabel, flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle);

  // * associations
  art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmhs(hitListHandle, evt, fTrackModuleLabel);
  art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

  // charge from raw digits
  rawcharge = 0;
  /* Comment for now as it is too slow
     for (size_t i = 0; i<rawlist.size(); ++i){
     if (fGeom->SignalType(rawlist[i]->Channel()) == geo::kCollection){
     double pedestal = rawlist[i]->GetPedestal();
     for (size_t j = 0; j<rawlist[i]->NADC(); ++j){
     rawcharge += rawlist[i]->ADC(j)-pedestal;
     }
     }
     }
     */
  //charge from wires
  wirecharge = 0;
  for (size_t i = 0; i<wirelist.size(); ++i){
    if (fGeom->SignalType(wirelist[i]->Channel()) == geo::kCollection){
      const recob::Wire::RegionsOfInterest_t& signalROI = wirelist[i]->SignalROI();
      for(const auto& range : signalROI.get_ranges()){
        const std::vector<float>& signal = range.data();
        raw::TDCtick_t roiFirstBinTick = range.begin_index();
        for (size_t j = 0; j<signal.size(); ++j){
          wirecharge += signal[j]*exp((j+roiFirstBinTick)*0.5/taulife);
        }
      }
    }
  }

  //hit information
  nhits = hitlist.size();
  nhits_stored = std::min(nhits, kMaxHits);
  for (int i = 0; i < nhits && i < kMaxHits ; ++i){//loop over hits
    hit_channel[i] = hitlist[i]->Channel();
    hit_plane[i]   = hitlist[i]->WireID().Plane;
    hit_wire[i]    = hitlist[i]->WireID().Wire;
    hit_tpc[i]     = hitlist[i]->WireID().TPC;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Integral();
    hit_summedADC[i] = hitlist[i]->SummedADC();
    hit_startT[i] = hitlist[i]->PeakTimeMinusRMS();
    hit_endT[i] = hitlist[i]->PeakTimePlusRMS();

  }

  //track information
  ntracks_reco=tracklist.size();

  TVector3 larStart;
  TVector3 larEnd;
  for(int i=0; i<std::min(int(tracklist.size()),kMaxTrack);++i){
    recob::Track::Point_t trackStart, trackEnd;
    std::tie(trackStart, trackEnd) = tracklist[i]->Extent(); 
    larStart = tracklist[i]->VertexDirection();
    larEnd = tracklist[i]->EndDirection();

    trkid[i]       = tracklist[i]->ID();
    trkstartx[i]      = trackStart.X();
    trkstarty[i]      = trackStart.Y();
    trkstartz[i]      = trackStart.Z();
    trkendx[i]        = trackEnd.X();
    trkendy[i]        = trackEnd.Y();
    trkendz[i]        = trackEnd.Z();
    trkstartdcosx[i]  = larStart[0];
    trkstartdcosy[i]  = larStart[1];
    trkstartdcosz[i]  = larStart[2];
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    trklen[i]         = tracklist[i]->Length();
    if (fmthm.isValid()){
      auto vhit = fmthm.at(i);
      auto vmeta = fmthm.data(i);
      for (size_t h = 0; h < vhit.size(); ++h){
        if (vhit[h].key()<kMaxHits){
          hit_trkkey[vhit[h].key()] = tracklist[i].key();
          if (vmeta[h]->Dx()){
            hit_dQds[vhit[h].key()] = vhit[h]->Integral()*fCalorimetryAlg.LifetimeCorrection(vhit[h]->PeakTime())/vmeta[h]->Dx();
            hit_dEds[vhit[h].key()] = fCalorimetryAlg.dEdx_AREA(vhit[h], vmeta[h]->Dx());
          }
          hit_resrange[vhit[h].key()] = tracklist[i]->Length(vmeta[h]->Index());
        }
      }//loop over all hits
    }//fmthm is valid
    else if (fmth.isValid()){
      std::vector< art::Ptr<recob::Hit> > vhit = fmth.at(i);
      for (size_t h = 0; h < vhit.size(); ++h){
        if (vhit[h].key()<kMaxHits){
          hit_trkkey[vhit[h].key()] = tracklist[i].key();
        }
      }
    }
    if (fmcal.isValid()){
      unsigned maxnumhits = 0;
      std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
      for (auto const& calo : calos){
        if (calo->PlaneID().isValid){
          trkke[i][calo->PlaneID().Plane] = calo->KineticEnergy();
          if (calo->dEdx().size()>maxnumhits){
            maxnumhits = calo->dEdx().size();
            trkbestplane[i] = calo->PlaneID().Plane;
          }
          double pida = 0;
          int used_trkres = 0;
          for (size_t ip = 0; ip<calo->dEdx().size(); ++ip){
            if (calo->ResidualRange()[ip]<30){
              pida += calo->dEdx()[ip]*pow(calo->ResidualRange()[ip],0.42);
              ++used_trkres;
            }
          }
          if (used_trkres) pida/=used_trkres;
          trkpida[i][calo->PlaneID().Plane] = pida;
        }
      }
    }
    if (!isdata&&fmth.isValid()){
      // Find true track for each reconstructed track
      int TrackID = 0;
      std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(i);

      std::map<int,double> trkide;
      for(size_t h = 0; h < allHits.size(); ++h){
        art::Ptr<recob::Hit> hit = allHits[h];
        std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
        for(size_t e = 0; e < TrackIDs.size(); ++e){
          trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
        }	    
      }
      // Work out which IDE despoited the most charge in the hit if there was more than one.
      double maxe = -1;
      double tote = 0;
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
        tote += ii->second;
        if ((ii->second)>maxe){
          maxe = ii->second;
          TrackID = ii->first;
        }
      }
      // Now have trackID, so get PdG code and T0 etc.
      const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
      if (particle){
        trkg4id[i] = TrackID;
        trkg4pdg[i] = particle->PdgCode();
        trkg4startx[i] = particle->Vx();
        trkg4starty[i] = particle->Vy();
        trkg4startz[i] = particle->Vz();
        float sum_energy = 0;
        int numhits = 0;
        //std::map<float,float> hite;
        double x = 0;
        double y = 0;
        double z = 0;
        double mindis = 1e10;
        //find the closest point to the neutrino vertex
        for(size_t h = 0; h < allHits.size(); ++h){
          art::Ptr<recob::Hit> hit = allHits[h];
          if (hit->WireID().Plane==2){
            std::vector<art::Ptr<recob::SpacePoint> > spts = fmhs.at(hit.key());
            if (spts.size()){
              double dis = sqrt(pow(spts[0]->XYZ()[0]-trkg4startx[i],2)+
                  pow(spts[0]->XYZ()[1]-trkg4starty[i],2)+
                  pow(spts[0]->XYZ()[2]-trkg4startz[i],2));
              if (dis<mindis){
                mindis = dis;
                x = spts[0]->XYZ()[0];
                y = spts[0]->XYZ()[1];
                z = spts[0]->XYZ()[2];
              }
            }
          }
        }
        for(size_t h = 0; h < allHits.size(); ++h){
          art::Ptr<recob::Hit> hit = allHits[h];
          if (hit->WireID().Plane==2){
            std::vector<art::Ptr<recob::SpacePoint> > spts = fmhs.at(hit.key());
            if (spts.size()){
              if (sqrt(pow(spts[0]->XYZ()[0]-x,2)+
                    pow(spts[0]->XYZ()[1]-y,2)+
                    pow(spts[0]->XYZ()[2]-z,2))<3){
                std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
                float toten = 0;
                for(size_t e = 0; e < TrackIDs.size(); ++e){
                  //sum_energy += TrackIDs[e].energy;
                  toten+=TrackIDs[e].energy;
                }
                if (toten){
                  sum_energy += toten;
                  ++numhits;
                }
              }
            }
          }
        }

        float pitch = 0;
        float dis1 = sqrt(pow(trkstartx[i]-trkg4startx[i],2)+
            pow(trkstarty[i]-trkg4starty[i],2)+
            pow(trkstartz[i]-trkg4startz[i],2));
        float dis2 = sqrt(pow(trkendx[i]-trkg4startx[i],2)+
            pow(trkendy[i]-trkg4starty[i],2)+
            pow(trkendz[i]-trkg4startz[i],2));
        if (dis1<dis2){
          try{
            pitch = lar::util::TrackPitchInView(*(tracklist[i]),geo::kZ,0);
          }
          catch(...){
            pitch = 0;
          }
        }
        else{
          try{
            pitch = lar::util::TrackPitchInView(*(tracklist[i]),geo::kZ,tracklist[i]->NumberTrajectoryPoints()-1);
          }
          catch(...){
            pitch = 0;
          }
        }
        if ( pitch && numhits ) {
          trkg4initdedx[i] = sum_energy/(numhits*pitch);
        }
        else{
          trkg4initdedx[i] = 0;
        }
      }//if (particle)
    }//MC
  }

  //vertex information
  nvtx = vtxlist.size();
  for (int i = 0; i < nvtx && i < kMaxVertices ; ++i){//loop over hits
    Double_t xyz[3] = {};
    vtxlist[i]->XYZ(xyz);
    for (size_t j = 0; j<3; ++j) vtx[i][j] = xyz[j];
  }

  //shower information
  if (shwListHandle.isValid()){
    art::FindManyP<recob::Hit> fmsh(shwListHandle, evt, fShowerModuleLabel);

    nshws = shwlist.size();

    for (int i = 0; i<std::min(int(shwlist.size()),kMaxShower); ++i){
      shwid[i] = shwlist[i]->ID();
      shwdcosx[i] = shwlist[i]->Direction().X(); 
      shwdcosy[i] = shwlist[i]->Direction().Y(); 
      shwdcosz[i] = shwlist[i]->Direction().Z(); 
      shwstartx[i] = shwlist[i]->ShowerStart().X();
      shwstarty[i] = shwlist[i]->ShowerStart().Y();
      shwstartz[i] = shwlist[i]->ShowerStart().Z();
      for (size_t j = 0; j<(shwlist[i]->Energy()).size(); ++j){
        shwenergy[i][j] = shwlist[i]->Energy()[j];
      }      
      for (size_t j = 0; j<(shwlist[i]->dEdx()).size(); ++j){
        shwdedx[i][j] = shwlist[i]->dEdx()[j];
      }
      shwbestplane[i] = shwlist[i]->best_plane();
      if (fmsh.isValid()){
        auto vhit = fmsh.at(i);
        for (size_t h = 0; h < vhit.size(); ++h){
          if (vhit[h].key()<kMaxHits){
            hit_shwkey[vhit[h].key()] = shwlist[i].key();
          }
        }
      }
      if (!isdata&&fmsh.isValid()){
        // Find true track for each reconstructed track
        int TrackID = 0;
        std::vector< art::Ptr<recob::Hit> > allHits = fmsh.at(i);
        std::map<int,double> trkide;
        for(size_t h = 0; h < allHits.size(); ++h){
          art::Ptr<recob::Hit> hit = allHits[h];
          std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
          for(size_t e = 0; e < TrackIDs.size(); ++e){
            trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
          }	    
        }
        // Work out which IDE despoited the most charge in the hit if there was more than one.
        double maxe = -1;
        double tote = 0;
        for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
          tote += ii->second;
          if ((ii->second)>maxe){
            maxe = ii->second;
            TrackID = ii->first;
          }
        }
        // Now have trackID, so get PdG code and T0 etc.
        const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
        if (particle){
          shwg4id[i] = TrackID;
        }
      }
    }
  }

  // flash information
  flash_total = flashlist.size();
  for ( int f = 0; f < std::min(flash_total,kMaxHits); ++f ) {
    flash_time[f]      = flashlist[f]->Time();
    flash_width[f]     = flashlist[f]->TimeWidth();
    flash_abstime[f]   = flashlist[f]->AbsTime();
    flash_YCenter[f]   = flashlist[f]->YCenter();
    flash_YWidth[f]    = flashlist[f]->YWidth();
    flash_ZCenter[f]   = flashlist[f]->ZCenter();
    flash_ZWidth[f]    = flashlist[f]->ZWidth();
    flash_TotalPE[f]   = flashlist[f]->TotalPE();
  }

  if (!isdata){

    // * MC truth information
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
    std::vector<art::Ptr<simb::MCFlux> > fluxlist;
    if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
      art::fill_ptr_vector(fluxlist, mcfluxListHandle);


    mcevts_truth=mclist.size();
    if (mcevts_truth){
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      if (mctruth->Origin() == simb::kBeamNeutrino){
        nuPDG_truth  = mctruth->GetNeutrino().Nu().PdgCode();
        ccnc_truth   = mctruth->GetNeutrino().CCNC();
        mode_truth   = mctruth->GetNeutrino().Mode();
        Q2_truth     = mctruth->GetNeutrino().QSqr();
        W_truth      = mctruth->GetNeutrino().W();
        X_truth      = mctruth->GetNeutrino().X();
        Y_truth      = mctruth->GetNeutrino().Y();
        hitnuc_truth = mctruth->GetNeutrino().HitNuc();
        target_truth = mctruth->GetNeutrino().Target();
        enu_truth    = mctruth->GetNeutrino().Nu().E();
        nuvtxx_truth = mctruth->GetNeutrino().Nu().Vx();
        nuvtxy_truth = mctruth->GetNeutrino().Nu().Vy();
        nuvtxz_truth = mctruth->GetNeutrino().Nu().Vz();
        if (mctruth->GetNeutrino().Nu().P()){
          nu_dcosx_truth = mctruth->GetNeutrino().Nu().Px()/mctruth->GetNeutrino().Nu().P();
          nu_dcosy_truth = mctruth->GetNeutrino().Nu().Py()/mctruth->GetNeutrino().Nu().P();
          nu_dcosz_truth = mctruth->GetNeutrino().Nu().Pz()/mctruth->GetNeutrino().Nu().P();
        }
        lep_mom_truth = mctruth->GetNeutrino().Lepton().P();
        if (mctruth->GetNeutrino().Lepton().P()){
          lep_dcosx_truth = mctruth->GetNeutrino().Lepton().Px()/mctruth->GetNeutrino().Lepton().P();
          lep_dcosy_truth = mctruth->GetNeutrino().Lepton().Py()/mctruth->GetNeutrino().Lepton().P();
          lep_dcosz_truth = mctruth->GetNeutrino().Lepton().Pz()/mctruth->GetNeutrino().Lepton().P();
        }

        if (mctruth->NParticles()){
          simb::MCParticle particle = mctruth->GetParticle(0);
          t0_truth = particle.T();
        }


        float mindist2 = 9999; // cm;
        TVector3 nuvtx(nuvtxx_truth, nuvtxy_truth, nuvtxz_truth);
        infidvol = insideFidVol(nuvtxx_truth, nuvtxy_truth, nuvtxz_truth); 
        //find the closest reco vertex to the neutrino mc truth
        if (infidvol)
        {
          // vertex is when at least two tracks meet
          for(size_t i = 0; i < vtxlist.size(); ++i){ // loop over vertices
            Double_t xyz[3] = {};
            vtxlist[i]->XYZ(xyz);
            TVector3 vtxreco(xyz);
            float dist2 = pma::Dist2(vtxreco, nuvtx);
            if (dist2 < mindist2)
            {
              mindist2 = dist2;
              vtxrecomc = std::sqrt(dist2);
              vtxrecomcx = vtxreco.X() - nuvtxx_truth;
              vtxrecomcy = vtxreco.Y() - nuvtxy_truth;
              vtxrecomcz = vtxreco.Z() - nuvtxz_truth;
            }
          }

          // two endpoints of tracks are somehow also vertices...
          for (size_t i = 0; i < tracklist.size(); ++i){ // loop over tracks
            float dist2 = pma::Dist2(tracklist[i]->Vertex(), nuvtx);
            if (dist2 < mindist2)
            {
              mindist2 = dist2;
              vtxrecomc = std::sqrt(dist2);
              vtxrecomcx = tracklist[i]->Vertex().X() - nuvtxx_truth;
              vtxrecomcy = tracklist[i]->Vertex().Y() - nuvtxy_truth;
              vtxrecomcz = tracklist[i]->Vertex().Z() - nuvtxz_truth;

            }
            dist2 = pma::Dist2(tracklist[i]->End(), nuvtx);
            if (dist2 < mindist2)
            {
              mindist2 = dist2;
              vtxrecomc = std::sqrt(dist2);
              vtxrecomcx = tracklist[i]->End().X() - nuvtxx_truth;
              vtxrecomcy = tracklist[i]->End().Y() - nuvtxy_truth;
              vtxrecomcz = tracklist[i]->End().Z() - nuvtxz_truth;

            }
          }
        }
      }//is neutrino
    } // if numver of mctruth si not zero

    if (fluxlist.size()){
      ptype_flux  = fluxlist[0]->fptype;
      pdpx_flux   = fluxlist[0]->fpdpx;
      pdpy_flux   = fluxlist[0]->fpdpy;
      pdpz_flux   = fluxlist[0]->fpdpz;
      pntype_flux = fluxlist[0]->fntype;
      vx_flux     = fluxlist[0]->fvx;
      vy_flux     = fluxlist[0]->fvy;
      vz_flux     = fluxlist[0]->fvz;
    }

    //save g4 particle information
    std::vector<const simb::MCParticle* > geant_part;

    // ### Looping over all the Geant4 particles from the BackTrackerService ###
    for(size_t p = 0; p < plist.size(); ++p) 
    {
      // ### Filling the vector with MC Particles ###
      geant_part.push_back(plist.Particle(p)); 
    }

    //std::cout<<"No of geant part= "<<geant_part.size()<<std::endl;

    // ### Setting a string for primary ###
    std::string pri("primary");

    int primary=0;
    int geant_particle=0;

    // ############################################################
    // ### Determine the number of primary particles from geant ###
    // ############################################################
    for( unsigned int i = 0; i < geant_part.size(); ++i ){
      geant_particle++;
      // ### Counting the number of primary particles ###
      if(geant_part[i]->Process()==pri)
      { primary++;}
    }//<---End i loop


    // ### Saving the number of primary particles ###
    no_primaries=primary;
    // ### Saving the number of Geant4 particles ###
    geant_list_size=geant_particle;

    // ### Looping over all the Geant4 particles ###
    for( unsigned int i = 0; i < geant_part.size(); ++i ){

      // ### If this particle is primary, set = 1 ###
      if(geant_part[i]->Process()==pri)
      {process_primary[i]=1;}
      // ### If this particle is not-primary, set = 0 ###
      else
      {process_primary[i]=0;}

      // ### Saving the particles mother TrackID ###
      Mother[i]=geant_part[i]->Mother();
      // ### Saving the particles TrackID ###
      TrackId[i]=geant_part[i]->TrackId();
      // ### Saving the PDG Code ###
      pdg[i]=geant_part[i]->PdgCode();
      // ### Saving the particles Energy ###
      Eng[i]=geant_part[i]->E();

      // ### Saving the Px, Py, Pz info ###
      Px[i]=geant_part[i]->Px();
      Py[i]=geant_part[i]->Py();
      Pz[i]=geant_part[i]->Pz();

      // ### Saving the Start and End Point for this particle ###
      StartPointx[i]=geant_part[i]->Vx();
      StartPointy[i]=geant_part[i]->Vy();
      StartPointz[i]=geant_part[i]->Vz();
      EndPointx[i]=geant_part[i]->EndPosition()[0];
      EndPointy[i]=geant_part[i]->EndPosition()[1];
      EndPointz[i]=geant_part[i]->EndPosition()[2];

      // ### Saving the processes for this particle ###
      //std::cout<<"finding proc"<<std::endl;
      G4Process.push_back( geant_part[i]->Process() );
      G4FinalProcess.push_back( geant_part[i]->EndProcess() );
      //std::cout<<"found proc"<<std::endl;
      //      std::cout << "ID " << TrackId[i] << ", pdg " << pdg[i] << ", Start X,Y,Z " << StartPointx[i] << ", " << StartPointy[i] << ", " << StartPointz[i]
      //		<< ", End XYZ " << EndPointx[i] << ", " << EndPointy[i] << ", " << EndPointz[i] << ", Start Proc " << G4Process[i] << ", End Proc " << G4FinalProcess[i]
      //		<< std::endl;

      // ### Saving the Start direction cosines for this particle ###
      Startdcosx[i] = geant_part[i]->Momentum(0).Px() / geant_part[i]->Momentum(0).P();
      Startdcosy[i] = geant_part[i]->Momentum(0).Py() / geant_part[i]->Momentum(0).P();
      Startdcosz[i] = geant_part[i]->Momentum(0).Pz() / geant_part[i]->Momentum(0).P();
      // ### Saving the number of Daughters for this particle ###
      NumberDaughters[i]=geant_part[i]->NumberDaughters();

    } //geant particles


  }//is neutrino

  if(fMakeAnaTree) 
    fTree->Fill();

}




bool dunemva::MVAAlg::insideFidVol(const double posX, const double posY, const double posZ) 
{

  double vtx[3] = {posX, posY, posZ};
  bool inside = false;

  geo::TPCID idtpc = fGeom->FindTPCAtPosition(vtx);

  if (fGeom->HasTPC(idtpc))
  {		
    const geo::TPCGeo& tpcgeo = fGeom->GetElement(idtpc);
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

    for (size_t c = 0; c < fGeom->Ncryostats(); c++)
    {
      const geo::CryostatGeo& cryostat = fGeom->Cryostat(c);
      for (size_t t = 0; t < cryostat.NTPC(); t++)
      {	
        const geo::TPCGeo& tpcg = cryostat.TPC(t);
        if (tpcg.MinX() < minx) minx = tpcg.MinX();
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX(); 
        if (tpcg.MinY() < miny) miny = tpcg.MinY();
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
      }
    }	


    //x
    double dista = fabs(minx - posX);
    double distb = fabs(posX - maxx); 
    if ((posX > minx) && (posX < maxx) &&
        (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    //y
    dista = fabs(maxy - posY);
    distb = fabs(posY - miny);
    if (inside && (posY > miny) && (posY < maxy) &&
        (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    else inside = false;

    //z
    dista = fabs(maxz - posZ);
    distb = fabs(posZ - minz);
    if (inside && (posZ > minz) && (posZ < maxz) &&
        (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    else inside = false;
  }

  return inside;
}

void dunemva::MVAAlg::ResetVars(){

  G4Process.clear();
  G4FinalProcess.clear();

  run = -9999;
  subrun = -9999;
  event = -9999;
  evttime = -9999;
  taulife = 0;
  isdata = -9999;

  ntracks_reco = 0;
  for (int i = 0; i < kMaxTrack; ++i){
    trkid[i] = -9999;
    trkstartx[i] = -9999;
    trkstarty[i] = -9999;
    trkstartz[i] = -9999;
    trkendx[i] = -9999;
    trkendy[i] = -9999;
    trkendz[i] = -9999;
    trkstartdcosx[i] = -9999;
    trkstartdcosy[i] = -9999;
    trkstartdcosz[i] = -9999;
    trkenddcosx[i] = -9999;
    trkenddcosy[i] = -9999;
    trkenddcosz[i] = -9999;
    trklen[i] = -9999;
    trkbestplane[i] = -9999;
    for (int j = 0; j<3; ++j){
      trkke[i][j] = -9999;
      trkpida[i][j] = -9999;
    }
    trkg4id[i] = -9999;
    trkg4pdg[i] = -9999;
    trkg4startx[i] = -9999;
    trkg4starty[i] = -9999;
    trkg4startz[i] = -9999;
    trkg4initdedx[i] = -9999;
  }

  nshws = 0;
  for (int i = 0; i<kMaxShower; ++i){
    shwid[i] = -9999;
    shwdcosx[i] = -9999;
    shwdcosy[i] = -9999;
    shwdcosz[i] = -9999;
    shwstartx[i] = -9999;
    shwstarty[i] = -9999;
    shwstartz[i] = -9999;
    for (int j = 0; j<3; ++j){
      shwenergy[i][j] = -9999;
      shwdedx[i][j] = -9999;
    }
    shwbestplane[i] = -9999;
    trkg4id[i] = -9999;
  }

  flash_total = 0;
  for (int f = 0; f < kMaxFlash; ++f) {
    flash_time[f]    = -9999;
    flash_width[f]   = -9999;
    flash_abstime[f] = -9999;
    flash_YCenter[f] = -9999;
    flash_YWidth[f]  = -9999;
    flash_ZCenter[f] = -9999;
    flash_ZWidth[f]  = -9999;
    flash_TotalPE[f] = -9999;
  }

  nhits = 0;
  nhits_stored = 0;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -9999;
    hit_wire[i] = -9999;
    hit_tpc[i] = -9999;
    hit_channel[i] = -9999;
    hit_peakT[i] = -9999;
    hit_charge[i] = -9999;
    hit_summedADC[i] = -9999;
    hit_startT[i] = -9999;
    hit_endT[i] = -9999;
    hit_trkkey[i] = -9999;
    hit_dQds[i] = -9999;
    hit_dEds[i] = -9999;
    hit_resrange[i] = -9999;
    hit_shwkey[i] = -9999;
  }

  infidvol = 0;
  nvtx = 0;
  for (int i = 0; i<kMaxVertices; ++i){
    vtx[i][0] = -9999;
    vtx[i][1] = -9999;
    vtx[i][2] = -9999;
  }
  vtxrecomc = 9999;
  vtxrecomcx = 9999;
  vtxrecomcy = 9999;
  vtxrecomcz = 9999;

  mcevts_truth = -9999; 
  nuPDG_truth = -9999;  
  ccnc_truth = -9999;    
  mode_truth = -9999;   
  enu_truth = -9999;     
  Q2_truth = -9999;       
  W_truth = -9999;        
  X_truth = -9999;
  Y_truth = -9999;
  hitnuc_truth = -9999;
  target_truth = -9999;
  nuvtxx_truth = -9999;   
  nuvtxy_truth = -9999;   
  nuvtxz_truth = -9999;   
  nu_dcosx_truth = -9999; 
  nu_dcosy_truth = -9999; 
  nu_dcosz_truth = -9999; 
  lep_mom_truth = -9999;  
  lep_dcosx_truth = -9999;
  lep_dcosy_truth = -9999;
  lep_dcosz_truth = -9999;
  t0_truth = -9999;

  no_primaries = -99999;
  geant_list_size=-9999;
  for (int i = 0; i<kMaxPrimaries; ++i){
    pdg[i] = -99999;
    Eng[i] = -99999;
    Px[i] = -99999;
    Py[i] = -99999;
    Pz[i] = -99999;
    StartPointx[i] = -99999;
    StartPointy[i] = -99999;
    StartPointz[i] = -99999;
    EndPointx[i] = -99999;
    EndPointy[i] = -99999;
    EndPointz[i] = -99999;
    Startdcosx[i]= -99999;
    Startdcosy[i]= -99999;
    Startdcosz[i]= -99999;
    NumberDaughters[i] = -99999;
    Mother[i] = -99999;
    TrackId[i] = -99999;
    process_primary[i] = -99999;
  }

  ptype_flux = -99999;
  pdpx_flux = -99999;
  pdpy_flux = -99999;
  pdpz_flux = -99999;
  pntype_flux = -99999;
  vx_flux = -99999;
  vy_flux = -99999;
  vz_flux = -99999;

  isinfidvol = 0;
  isinfidvoltruth = 0;
  oscpro = 0;
}


void dunemva::MVAAlg::MakeTree(){

  fTree = tfs->make<TTree>("nueana","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/F");
  fTree->Branch("taulife",&taulife,"taulife/F");
  fTree->Branch("isdata",&isdata,"isdata/S");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("trkid",trkid,"trkid[ntracks_reco]/I");
  fTree->Branch("trkstartx",trkstartx,"trkstartx[ntracks_reco]/F");
  fTree->Branch("trkstarty",trkstarty,"trkstarty[ntracks_reco]/F");
  fTree->Branch("trkstartz",trkstartz,"trkstartz[ntracks_reco]/F");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/F");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/F");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/F");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/F");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/F");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/F");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/F");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/F");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/F");
  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/F");
  fTree->Branch("trkbestplane",trkbestplane,"trkbestplane[ntracks_reco]/I");
  fTree->Branch("trkke",trkke,"trkke[ntracks_reco][3]/F");
  fTree->Branch("trkpida",trkpida,"trkpida[ntracks_reco][3]/F");
  fTree->Branch("trkg4id",trkg4id,"trkg4id[ntracks_reco]/I");
  fTree->Branch("trkg4pdg",trkg4pdg,"trkg4pdg[ntracks_reco]/I");
  fTree->Branch("trkg4startx",trkg4startx,"trkg4startx[ntracks_reco]/F");
  fTree->Branch("trkg4starty",trkg4starty,"trkg4starty[ntracks_reco]/F");
  fTree->Branch("trkg4startz",trkg4startz,"trkg4startz[ntracks_reco]/F");
  fTree->Branch("trkg4initdedx",trkg4initdedx,"trkg4initdedx[ntracks_reco]/F");
  fTree->Branch("nshws",&nshws,"nshws/I");
  fTree->Branch("shwid",shwid,"shwid[nshws]/I");
  fTree->Branch("shwdcosx",shwdcosx,"shwdcosx[nshws]/F");
  fTree->Branch("shwdcosy",shwdcosy,"shwdcosy[nshws]/F");
  fTree->Branch("shwdcosz",shwdcosz,"shwdcosz[nshws]/F");
  fTree->Branch("shwstartx",shwstartx,"shwstartx[nshws]/F");
  fTree->Branch("shwstarty",shwstarty,"shwstarty[nshws]/F");
  fTree->Branch("shwstartz",shwstartz,"shwstartz[nshws]/F");
  fTree->Branch("shwenergy",shwenergy,"shwenergy[nshws][3]/F");
  fTree->Branch("shwdedx",shwdedx,"shwdedx[nshws][3]/F");
  fTree->Branch("shwbestplane",shwbestplane,"shwbestplane[nshws]/I");
  fTree->Branch("shwg4id",shwg4id,"shwg4id[nshws]/I");
  fTree->Branch("flash_total"  ,&flash_total ,"flash_total/I");
  fTree->Branch("flash_time"   ,flash_time   ,"flash_time[flash_total]/F");
  fTree->Branch("flash_width"  ,flash_width  ,"flash_width[flash_total]/F");
  fTree->Branch("flash_abstime",flash_abstime,"flash_abstime[flash_total]/F");
  fTree->Branch("flash_YCenter",flash_YCenter,"flash_YCenter[flash_total]/F");
  fTree->Branch("flash_YWidth" ,flash_YWidth ,"flash_YWidth[flash_total]/F");
  fTree->Branch("flash_ZCenter",flash_ZCenter,"flash_ZCenter[flash_total]/F");
  fTree->Branch("flash_ZWidth" ,flash_ZWidth ,"flash_ZWidth[flash_total]/F");
  fTree->Branch("flash_TotalPE",flash_TotalPE,"flash_TotalPE[flash_total]/F");
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("nhits_stored",&nhits_stored,"nhits_stored/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits_stored]/S");
  fTree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits_stored]/S");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits_stored]/S");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits_stored]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits_stored]/F");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits_stored]/F");
  fTree->Branch("hit_summedADC",hit_summedADC,"hit_summedADC[nhits_stored]/F");
  fTree->Branch("hit_startT",hit_startT,"hit_startT[nhits_stored]/F");
  fTree->Branch("hit_endT",hit_endT,"hit_endT[nhits_stored]/F");
  fTree->Branch("hit_trkkey",hit_trkkey,"hit_trkkey[nhits_stored]/I");
  fTree->Branch("hit_dQds",hit_dQds,"hit_dQds[nhits_stored]/F");
  fTree->Branch("hit_dEds",hit_dEds,"hit_dEds[nhits_stored]/F");
  fTree->Branch("hit_resrange",hit_resrange,"hit_resrange[nhits_stored]/F");
  fTree->Branch("hit_shwkey",hit_shwkey,"hit_shwkey[nhits_stored]/I");
  fTree->Branch("infidvol",&infidvol,"infidvol/I");
  fTree->Branch("nvtx",&nvtx,"nvtx/S");
  fTree->Branch("vtx",vtx,"vtx[nvtx][3]/F");
  fTree->Branch("vtxrecomc",&vtxrecomc,"vtxrecomc/F");
  fTree->Branch("vtxrecomcx",&vtxrecomcx,"vtxrecomcx/F");
  fTree->Branch("vtxrecomcy",&vtxrecomcy,"vtxrecomcy/F");
  fTree->Branch("vtxrecomcz",&vtxrecomcz,"vtxrecomcz/F");
  fTree->Branch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
  fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  fTree->Branch("enu_truth",&enu_truth,"enu_truth/F");
  fTree->Branch("Q2_truth",&Q2_truth,"Q2_truth/F");
  fTree->Branch("W_truth",&W_truth,"W_truth/F");
  fTree->Branch("X_truth",&X_truth,"X_truth/F");
  fTree->Branch("Y_truth",&Y_truth,"Y_truth/F");
  fTree->Branch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
  fTree->Branch("target_truth",&target_truth,"target_truth/I");
  fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/F");
  fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/F");
  fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/F");
  fTree->Branch("nu_dcosx_truth",&nu_dcosx_truth,"nu_dcosx_truth/F");
  fTree->Branch("nu_dcosy_truth",&nu_dcosy_truth,"nu_dcosy_truth/F");
  fTree->Branch("nu_dcosz_truth",&nu_dcosz_truth,"nu_dcosz_truth/F");
  fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/F");
  fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/F");
  fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/F");
  fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/F");
  fTree->Branch("t0_truth",&t0_truth,"t0_truth/F");
  fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",&geant_list_size,"geant_list_size/I");
  fTree->Branch("pdg",pdg,"pdg[geant_list_size]/I");
  fTree->Branch("Eng",Eng,"Eng[geant_list_size]/F");
  fTree->Branch("Px",Px,"Px[geant_list_size]/F");
  fTree->Branch("Py",Py,"Py[geant_list_size]/F");
  fTree->Branch("Pz",Pz,"Pz[geant_list_size]/F");
  fTree->Branch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
  fTree->Branch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
  fTree->Branch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
  fTree->Branch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
  fTree->Branch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
  fTree->Branch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
  fTree->Branch("Startdcosx",Startdcosx,"Startdcosx[geant_list_size]/F");
  fTree->Branch("Startdcosy",Startdcosy,"Startdcosy[geant_list_size]/F");
  fTree->Branch("Startdcosz",Startdcosz,"Startdcosz[geant_list_size]/F");
  fTree->Branch("N(((((((((umberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
  fTree->Branch("Mother",Mother,"Mother[geant_list_size]/I");
  fTree->Branch("TrackId",TrackId,"TrackId[geant_list_size]/I");
  fTree->Branch("process_primary",process_primary,"process_primary[geant_list_size]/I");
  fTree->Branch("G4Process",&G4Process);//,"G4Process[geant_list_size]");
  fTree->Branch("G4FinalProcess",&G4FinalProcess);//,"G4FinalProcess[geant_list_size]");
  fTree->Branch("ptype_flux",&ptype_flux,"ptype_flux/I");
  fTree->Branch("pdpx_flux",&pdpx_flux,"pdpx_flux/F");
  fTree->Branch("pdpy_flux",&pdpy_flux,"pdpy_flux/F");
  fTree->Branch("pdpz_flux",&pdpz_flux,"pdpz_flux/F");
  fTree->Branch("pntype_flux",&pntype_flux,"pntype_flux/I");
  fTree->Branch("vx_flux",&vx_flux,"vx_flux/F");
  fTree->Branch("vy_flux",&vy_flux,"vy_flux/F");
  fTree->Branch("vz_flux",&vz_flux,"vz_flux/F");

}
