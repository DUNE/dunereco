////////////////////////////////////////////////////////////////////////
//
// \file MVASelect_module.cc
//
// tylerdalion@gmail.com
//
///////////////////////////////////////////////////////////////////////

#ifndef MVASelect_H
#define MVASelect_H

// Generic C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

#include "NuReweight/art/NuReweight.h"

#include "SimulationBase/GTruth.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCFlux.h"
#include "larcore/SummaryData/POTSummary.h"
#include "dune/FDSensOpt/MVAAlg/MVAAlg.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"

constexpr int knRwgts = 41; // -2 to 2 sigma in 0.1 steps

namespace dunemva {

class MVASelect : public art::EDAnalyzer {
  
public:
  
  explicit MVASelect(fhicl::ParameterSet const& pset);
  virtual ~MVASelect();
  void beginJob() override;
  void beginSubRun(const art::SubRun& sr) override;
  void endSubRun(const art::SubRun& sr) override;
  void reconfigure(fhicl::ParameterSet const& pset) override;
  void analyze(art::Event const & evt) override;

  
private:
  
  MVAAlg fMVAAlg;
  std::vector<double> fMVAResult;
  std::vector<std::string> fMVAMethods;
  int fnMethods;
  double fmvaResults[10];

  bool fReweight;
  unsigned int fRun,fSubrun,fEvent;
  double fWeight;
  TTree* fTree;  

  void FillINuke( rwgt::ReweightLabel_t label, double sig, double wgt);
  std::vector< rwgt::ReweightLabel_t > fINukeLabel;

  // INuke reweights
  double fSigs[knRwgts];
  double finuke_MFP_pi[knRwgts];
  double finuke_MFP_N[knRwgts];
  double finuke_FrCEx_pi[knRwgts];
  double finuke_FrElas_pi[knRwgts];
  double finuke_FrInel_pi[knRwgts];
  double finuke_FrAbs_pi[knRwgts];
  double finuke_FrPiProd_pi[knRwgts];
  double finuke_FrCEx_N[knRwgts];
  double finuke_FrElas_N[knRwgts];
  double finuke_FrInel_N[knRwgts];
  double finuke_FrAbs_N[knRwgts];
  double finuke_FrPiProd_N[knRwgts];

  // variables for counting particles
  int nP;
  int nN;
  int nPip;
  int nPim;
  int nPi0;
  int nKp;
  int nKm;
  int nK0;
  int nEM;
  int nOtherHad;
  int nNucleus;
  int nUNKNOWN;  


  double fQ2; 
  double fEtrue; 

  int fIsCoh; // 1=is coherent, 0=otherwise
  int fCC;    // 1=is CC, 0=otherwise
  int fNC;    // 1=is NC, 0=otherwise
  int fEvClass_reco; // 0=NuMuCC, 1=NuECC, 2=NC
  int fNuPdg;
  int fBeamPdg;

  // trut nu vertex for fiducial cuts
  double nuvtxx_truth;
  double nuvtxy_truth;
  double nuvtxz_truth;

  int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  int fCCNC; // 0=CC, 1=NC

  double ftotPOT;


}; // class MVASelect

 
//------------------------------------------------------------------------------
MVASelect::MVASelect(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
    , fMVAAlg(pset)
{
  this->reconfigure(pset);

  fINukeLabel.emplace_back(rwgt::fReweightMFP_pi);
  fINukeLabel.emplace_back(rwgt::fReweightMFP_N);
  fINukeLabel.emplace_back(rwgt::fReweightFrCEx_pi);
  fINukeLabel.emplace_back(rwgt::fReweightFrElas_pi);
  fINukeLabel.emplace_back(rwgt::fReweightFrInel_pi);
  fINukeLabel.emplace_back(rwgt::fReweightFrAbs_pi);
  fINukeLabel.emplace_back(rwgt::fReweightFrPiProd_pi);
  fINukeLabel.emplace_back(rwgt::fReweightFrCEx_N);
  fINukeLabel.emplace_back(rwgt::fReweightFrElas_N);
  fINukeLabel.emplace_back(rwgt::fReweightFrInel_N);
  fINukeLabel.emplace_back(rwgt::fReweightFrAbs_N);
  fINukeLabel.emplace_back(rwgt::fReweightFrPiProd_N);
}

dunemva::MVASelect::~MVASelect(){}

//------------------------------------------------------------------------------
void MVASelect::reconfigure(fhicl::ParameterSet const& pset) 
{
  fMVAMethods=pset.get<std::vector<std::string> >("MVAMethods");
  fReweight=pset.get<bool>("Reweight");
  fMVAAlg.reconfigure(pset);
  return;
}


//------------------------------------------------------------------------------
void MVASelect::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree =tfs->make<TTree>("MVASelection","Results");

  fTree->Branch("run",         &fRun,        "run/I");
  fTree->Branch("subrun",      &fSubrun,     "subrun/I");
  fTree->Branch("event",       &fEvent,      "event/I");
  fTree->Branch("nmvamethods",  fnMethods,   "nmvamethods/I");
  fTree->Branch("mvaresults",   fmvaResults, "mvaresults[10]/D");
  fTree->Branch("weight",      &fWeight,     "weight/D");

  // particle counts
  fTree->Branch("nP",        &nP,         "nP/I");
  fTree->Branch("nN",        &nN,         "nN/I");
  fTree->Branch("nipip",     &nPip,       "nipip/I");
  fTree->Branch("nipim",     &nPim,       "nipim/I");
  fTree->Branch("nipi0",     &nPi0,       "nipi0/I");
  fTree->Branch("nikp",      &nKp,        "nikp/I");
  fTree->Branch("nikm",      &nKm,        "nikm/I");
  fTree->Branch("nik0",      &nK0,        "nik0/I");
  fTree->Branch("niem",      &nEM,        "niem/I");
  fTree->Branch("niother",   &nOtherHad,  "niother/I");
  fTree->Branch("nNucleus",  &nNucleus,   "nNucleus/I");
  fTree->Branch("nUNKNOWN",  &nUNKNOWN,   "nUNKNOWN/I");

  fTree->Branch("Q2",           &fQ2,           "Q2/D");
  fTree->Branch("Ev",           &fEtrue,        "Ev/D");
  fTree->Branch("Ev_reco",      &fEtrue,        "Ev_reco/D");
  fTree->Branch("EvClass_reco", &fEvClass_reco, "EvClass_reco/I");
  fTree->Branch("coh",          &fIsCoh,        "coh/I");
  fTree->Branch("cc",           &fCC,           "cc/I");
  fTree->Branch("nc",           &fNC,           "nc/I");
  fTree->Branch("neu",          &fNuPdg,        "neu/I");
  fTree->Branch("beamPdg",      &fBeamPdg,      "beamPdg/I");

  fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/D");
  fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/D");
  fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/D");

  fTree->Branch("totpot",       &ftotPOT,       "totpot/D");

  fTree->Branch("sigs",             fSigs,             "sigs[41]/D");
  fTree->Branch("f_int_pi_mfp",     finuke_MFP_pi,     "f_int_pi_mfp[41]/D");
  fTree->Branch("f_int_n_mfp",      finuke_MFP_N,      "f_int_n_mfp[41]/D");
  fTree->Branch("f_int_pi_cex",     finuke_FrCEx_pi,   "f_int_pi_cex[41]/D");
  fTree->Branch("f_int_pi_elas",    finuke_FrElas_pi,  "f_int_pi_elas[41]/D");
  fTree->Branch("f_int_pi_inel",    finuke_FrInel_pi,  "f_int_pi_inel[41]/D");
  fTree->Branch("f_int_pi_abs",     finuke_FrAbs_pi,   "f_int_pi_abs[41]/D");
  fTree->Branch("f_int_pi_prod",    finuke_FrPiProd_pi,"f_int_pi_prod[41]/D");
  fTree->Branch("f_int_n_cex",      finuke_FrCEx_N,    "f_int_n_cex[41]/D");
  fTree->Branch("f_int_n_elas",     finuke_FrElas_N,   "f_int_n_elas[41]/D");
  fTree->Branch("f_int_n_inel",     finuke_FrInel_N,   "f_int_n_inel[41]/D");
  fTree->Branch("f_int_n_abs",      finuke_FrAbs_N,    "f_int_n_abs[41]/D");
  fTree->Branch("f_int_n_prod",     finuke_FrPiProd_N, "f_int_n_prod[41]/D");

  return;
}


//------------------------------------------------------------------------------
void MVASelect::beginSubRun(const art::SubRun& sr){
  art::Handle< sumdata::POTSummary > pots;
  if(sr.getByLabel("generator",pots))
    ftotPOT = pots->totpot;
  else
    ftotPOT = -1.;
}


//------------------------------------------------------------------------------
void MVASelect::analyze(art::Event const & evt)
{
  fRun = evt.id().run();
  fSubrun = evt.id().subRun();
  fEvent = evt.id().event();
  fMVAAlg.Run(evt,fMVAResult,fWeight);

  // fill branch variables
  fnMethods = fMVAResult.size();
  if(fnMethods==0)
    mf::LogWarning("MVASelect") << "Returning zero MVA methods";
  for(int i=0; i<fnMethods; i++){
    fmvaResults[i] = fMVAResult[i];
  }

  art::Handle< std::vector<simb::MCTruth> > mct;
  std::vector< art::Ptr<simb::MCTruth> > truth;
  if( evt.getByLabel("generator", mct) )
    art::fill_ptr_vector(truth, mct);
  else
    mf::LogWarning("MVASelect") << "No MCTruth.";

  art::Handle< std::vector<simb::MCFlux> > mcf;
  std::vector< art::Ptr<simb::MCFlux> > flux;
  if( evt.getByLabel("generator", mcf) )
    art::fill_ptr_vector(flux, mcf);
  else
    mf::LogWarning("MVASelect") << "No MCFlux.";

  art::Handle< std::vector<simb::GTruth> > gt;
  std::vector< art::Ptr<simb::GTruth> > gtru;
  if( evt.getByLabel("generator", gt) )
    art::fill_ptr_vector(gtru, gt);
  else
    mf::LogWarning("MVASelect") << "No GTruth.";


  for(size_t i=0; i<truth.size(); i++){
    
    if(i>1){
      mf::LogWarning("MVASelect") << "Skipping MC truth index " << i;
      continue;
    }

 
    fCCNC     = truth[i]->GetNeutrino().CCNC();  //0=CC 1=NC
    if(fCCNC==0){
      fCC = 1;
      fNC = 0;
    }
    else if(fCCNC==1){
      fCC = 0;
      fNC = 1;
    }
    fNuPdg    = truth[i]->GetNeutrino().Nu().PdgCode();
    fBeamPdg  = flux[i]->fntype;
    fMode     = truth[i]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    fEtrue    = truth[i]->GetNeutrino().Nu().E();
    fQ2       = truth[i]->GetNeutrino().QSqr();

    nuvtxx_truth = truth[i]->GetNeutrino().Nu().Vx();
    nuvtxy_truth = truth[i]->GetNeutrino().Nu().Vy();
    nuvtxz_truth = truth[i]->GetNeutrino().Nu().Vz();

    if(fMode==3)
      fIsCoh = 1;
    else
      fIsCoh = 0;
    
    // Cheat selection for now, will be refilled in script cutting on MVAResult
    fEvClass_reco = -1;
    if( fCCNC==1 )
      fEvClass_reco = 2;
    else if( fCCNC==0 && std::abs(fNuPdg)==12  )
      fEvClass_reco = 1;
    else if( fCCNC==0 && std::abs(fNuPdg)==14  )
      fEvClass_reco = 0;
    else
      mf::LogWarning("MVASelect") << "Bad EvClass_reco";

    nP        = 0;
    nN        = 0;
    nPip      = 0;
    nPim      = 0;
    nPi0      = 0;
    nKp       = 0;
    nKm       = 0;
    nK0       = 0;
    nEM       = 0;
    nOtherHad = 0;
    nNucleus  = 0;
    nUNKNOWN  = 0;

    for(int p=0; p<truth[i]->NParticles(); p++){
      if( truth[i]->GetParticle(p).StatusCode() == genie::kIStHadronInTheNucleus ){

	if      (truth[i]->GetParticle(p).PdgCode() == genie::kPdgProton  || 
		 truth[i]->GetParticle(p).PdgCode() == genie::kPdgAntiProton)          { nP++;       }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgNeutron || 
		 truth[i]->GetParticle(p).PdgCode() == genie::kPdgAntiNeutron)         { nN++;       }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgPiP)                 { nPip++;     }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgPiM)                 { nPim++;     }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgPi0)                 { nPi0++;     }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgKP)                  { nKp++;      }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgKM)                  { nKm++;      }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgK0    || 
		 truth[i]->GetParticle(p).PdgCode() == genie::kPdgAntiK0)              { nK0++;      }
	else if (truth[i]->GetParticle(p).PdgCode() == genie::kPdgGamma || 
		 truth[i]->GetParticle(p).PdgCode() == genie::kPdgElectron || 
		 truth[i]->GetParticle(p).PdgCode() == genie::kPdgPositron)            { nEM++;      }
	else if (genie::pdg::IsHadron(truth[i]->GetParticle(p).PdgCode()))             { nOtherHad++;}
	else if (genie::pdg::IsIon   (truth[i]->GetParticle(p).PdgCode()))             { nNucleus++; }
	else                                                                           { nUNKNOWN++; }
	
      }
    }

    unsigned int sigs_i = 0;
    double sig_step = 0.1;
    for(double sig=-2; sig<=2; sig+=sig_step){      
      if(knRwgts<=sigs_i)
	mf::LogError("MVASelect") << "too many sigma steps";
      fSigs[sigs_i] = sig;
      sigs_i++;
    }
    

    if(fReweight){ // takes a long time, only reweight if necessary
      rwgt::NuReweight *rwt;
      for(unsigned int r=0; r<fINukeLabel.size(); r++){
	
	for(double s=-2; s<=2; s+=sig_step){
	  
	  rwt = new rwgt::NuReweight();
	  rwt->ConfigureINuke();
	  rwt->ReweightIntraNuke(fINukeLabel[r],s);
	  double wgt = rwt->CalcWeight(*(truth[i]), *(gtru[i]));
	  if(wgt>10)
	    mf::LogVerbatim("MVASelect") << "High weight: " << wgt;
	  this->FillINuke(fINukeLabel[r],s,wgt);
	  
	} // all sigma steps s
      } // all reweights r
    } // if reweighting
    
  } // loop through MC truth i
  

  fTree->Fill();
  fMVAResult.clear();

  return;
}

//------------------------------------------------------------------------------
void MVASelect::FillINuke( rwgt::ReweightLabel_t label, double sig, double wgt){

  unsigned int sig_i = UINT_MAX;
  for(unsigned int a=0; a<knRwgts; a++)
    if( fSigs[a] == sig )
      sig_i = a;
	
  if(sig_i==UINT_MAX)
    mf::LogError("MVASelect") << "bad sigma index";

  if(      label == rwgt::fReweightMFP_pi )
    finuke_MFP_pi[sig_i] = wgt;
  else if( label == rwgt::fReweightMFP_N )
    finuke_MFP_N[sig_i] = wgt;
  else if( label == rwgt::fReweightFrCEx_pi )
    finuke_FrCEx_pi[sig_i] = wgt;
  else if( label == rwgt::fReweightFrElas_pi )
    finuke_FrElas_pi[sig_i] = wgt;
  else if( label == rwgt::fReweightFrInel_pi )
    finuke_FrInel_pi[sig_i] = wgt;
  else if( label == rwgt::fReweightFrAbs_pi )
    finuke_FrAbs_pi[sig_i] = wgt;
  else if( label == rwgt::fReweightFrPiProd_pi )
    finuke_FrPiProd_pi[sig_i] = wgt;
  else if( label == rwgt::fReweightFrCEx_N )
    finuke_FrCEx_N[sig_i] = wgt;
  else if( label == rwgt::fReweightFrElas_N )
    finuke_FrElas_N[sig_i] = wgt;
  else if( label == rwgt::fReweightFrInel_N )
    finuke_FrInel_N[sig_i] = wgt;
  else if( label == rwgt::fReweightFrAbs_N )
    finuke_FrAbs_N[sig_i] = wgt;
  else if( label == rwgt::fReweightFrPiProd_N )
    finuke_FrPiProd_N[sig_i] = wgt;

}

//------------------------------------------------------------------------------
void MVASelect::endSubRun(const art::SubRun& sr){
  fMVAAlg.endSubRun(sr);
}
  
DEFINE_ART_MODULE(MVASelect)

} // namespace dunemva

#endif // MVASelect_H


  /*

    ReweightLabel_t:

    fReweightMFP_pi = genie::rew::kINukeTwkDial_MFP_pi,           ///< tweak mean free path for pions
    fReweightMFP_N = genie::rew::kINukeTwkDial_MFP_N,             ///< tweak mean free path for nucleons
    fReweightFrCEx_pi = genie::rew::kINukeTwkDial_FrCEx_pi,       ///< tweak charge exchange probability for pions, for given total rescattering probability
    fReweightFrElas_pi = genie::rew::kINukeTwkDial_FrElas_pi,     ///< tweak elastic   probability for pions, for given total rescattering probability
    fReweightFrInel_pi = genie::rew::kINukeTwkDial_FrInel_pi,     ///< tweak inelastic probability for pions, for given total rescattering probability
    fReweightFrAbs_pi = genie::rew::kINukeTwkDial_FrAbs_pi,       ///< tweak absorption probability for pions, for given total rescattering probability
    fReweightFrPiProd_pi = genie::rew::kINukeTwkDial_FrPiProd_pi, ///< tweak pion production probability for pions, for given total rescattering probability
    fReweightFrCEx_N = genie::rew::kINukeTwkDial_FrCEx_N,         ///< tweak charge exchange probability for nucleons, for given total rescattering probability
    fReweightFrElas_N = genie::rew::kINukeTwkDial_FrElas_N,       ///< tweak elastic    probability for nucleons, for given total rescattering probability
    fReweightFrInel_N = genie::rew::kINukeTwkDial_FrInel_N,       ///< tweak inelastic  probability for nucleons, for given total rescattering probability
    fReweightFrAbs_N = genie::rew::kINukeTwkDial_FrAbs_N,         ///< tweak absorption probability for nucleons, for given total rescattering probability
    fReweightFrPiProd_N = genie::rew::kINukeTwkDial_FrPiProd_N,   ///< tweak pion production probability for nucleons, for given total rescattering probability

    kINukeTwkDial_MFP_pi            20%
    kINukeTwkDial_MFP_N             20%
    kINukeTwkDial_FrAbs_pi          20%
    kINukeTwkDial_FrCEx_pi          50%
    kINukeTwkDial_FrElas_pi         10%
    kINukeTwkDial_FrInel_pi         40%
    kINukeTwkDial_FrPiProd_pi       20%
    kINukeTwkDial_FrAbs_N           20%
    kINukeTwkDial_FrCEx_N           50%
    kINukeTwkDial_FrElas_N          30%
    kINukeTwkDial_FrInel_N          40%
    kINukeTwkDial_FrPiProd_N        20%

  */
