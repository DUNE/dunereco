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
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

#include "nutools/NuReweight/art/NuReweight.h"
#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "dune/FDSensOpt/MVAAlg/MVAAlg.h"

#include "dune/FDSensOpt/FDSensOptData/MVASelectPID.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"

namespace dunemva {

  class MVASelect : public art::EDProducer {

    public:

      explicit MVASelect(fhicl::ParameterSet const& pset);
      virtual ~MVASelect();
      void beginJob() override;
      void beginSubRun(art::SubRun& sr) override;
      void endSubRun(art::SubRun& sr) override;
      void reconfigure(fhicl::ParameterSet const& pset) ;
      void produce(art::Event& evt) override;


    private:

      MVAAlg fMVAAlg;
      std::string fMVAMethod;

      bool fSelNuE;
      bool fSelNuMu;

  }; // class MVASelect


  //------------------------------------------------------------------------------
  MVASelect::MVASelect(fhicl::ParameterSet const& pset)
    : fMVAAlg(pset)
  {
    produces<dunemva::MVASelectPID>();

    this->reconfigure(pset);
  }

  dunemva::MVASelect::~MVASelect(){}

  //------------------------------------------------------------------------------
  void MVASelect::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fMVAMethod=pset.get< std::string >("MVAMethod");
    fMVAAlg.reconfigure(pset);

    if(pset.get<std::string>("Select") == "nue"){
      fSelNuE  = true;
      fSelNuMu = false;
    } else if(pset.get<std::string>("Select") == "numu"){
      fSelNuE  = false;
      fSelNuMu = true;
    }

  }


  //------------------------------------------------------------------------------
  void MVASelect::beginJob()
  {
  }

  //------------------------------------------------------------------------------
  /*
     TH2D* MVASelect::MakeRespoHist(TString name){
     return fTfs->make<TH2D>( name, name,
     bins_true_E, bin_edges_true_E,
     n_respo_bins,bin_edges_respo);
     }
     */

  //------------------------------------------------------------------------------
  void MVASelect::beginSubRun(art::SubRun& sr){
  }


  //------------------------------------------------------------------------------
  void MVASelect::produce(art::Event& evt)
  {
    std::unique_ptr<dunemva::MVASelectPID> pidout = std::make_unique<MVASelectPID>();

    if(fSelNuE) pidout->selectMode = 12;
    if(fSelNuMu) pidout->selectMode = 14;

    // Weight seems to have osc probability etc in it. Not suitable for a reco
    // module, drop that output on the floor
    double weight;
    fMVAAlg.Run(evt, pidout->pid, weight);

    pidout->evtcharge = fMVAAlg.evtcharge;
    pidout->rawcharge = fMVAAlg.rawcharge;
    pidout->wirecharge = fMVAAlg.wirecharge;

    pidout->ntrack         = fMVAAlg.ntrack;
    pidout->avgtrklength   = fMVAAlg.avgtrklength;
    pidout->maxtrklength   = fMVAAlg.maxtrklength;
    pidout->trkdedx        = fMVAAlg.trkdedx;
    pidout->trkrch         = fMVAAlg.trkrch;
    pidout->trkrt          = fMVAAlg.trkrt;
    pidout->trkfr          = fMVAAlg.trkfr;
    pidout->trkpida_save   = fMVAAlg.trkpida_save;
    pidout->nshower        = fMVAAlg.nshower;
    pidout->showerdedx     = fMVAAlg.showerdedx;
    pidout->eshower        = fMVAAlg.eshower;
    pidout->frshower       = fMVAAlg.frshower;
    pidout->nhitspershw    = fMVAAlg.nhitspershw;
    pidout->shwlength      = fMVAAlg.shwlength;
    pidout->shwmax         = fMVAAlg.shwmax;
    pidout->fract_5_wires  = fMVAAlg.fract_5_wires;
    pidout->fract_10_wires = fMVAAlg.fract_10_wires;
    pidout->fract_50_wires = fMVAAlg.fract_50_wires;
    pidout->fract_100_wires= fMVAAlg.fract_100_wires;
    pidout->shwdis         = fMVAAlg.shwdis;
    pidout->shwdisx        = fMVAAlg.shwdisx;
    pidout->shwdisy        = fMVAAlg.shwdisy;
    pidout->shwdisz        = fMVAAlg.shwdisz;
    pidout->shwcosx        = fMVAAlg.shwcosx;
    pidout->shwcosy        = fMVAAlg.shwcosy;
    pidout->shwcosz        = fMVAAlg.shwcosz;
    pidout->trkcosx        = fMVAAlg.trkcosx;
    pidout->trkcosy        = fMVAAlg.trkcosy;
    pidout->trkcosz        = fMVAAlg.trkcosz;
    pidout->et             = fMVAAlg.ET;

    evt.put(std::move(pidout));
  }

  //------------------------------------------------------------------------------
  void MVASelect::endSubRun(art::SubRun& sr){
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
