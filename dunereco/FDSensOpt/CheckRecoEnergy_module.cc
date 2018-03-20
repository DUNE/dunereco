////////////////////////////////////////////////////////////////////////
//
// \file CheckRecoEnergy_module.cc
//
// n.grant.3@warwick.ac.uk
//
///////////////////////////////////////////////////////////////////////

#ifndef CheckRecoEnergy_H
#define CheckRecoEnergy_H

// Generic C++ includes
#include <iostream>
#include <vector>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

#include "nutools/NuReweight/art/NuReweight.h"
#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"

namespace dune {

  class CheckRecoEnergy : public art::EDAnalyzer {

    public:

      explicit CheckRecoEnergy(fhicl::ParameterSet const& pset);
      virtual ~CheckRecoEnergy();
      void beginJob() override;
      void endJob() override;
      void beginSubRun(const art::SubRun& sr) override;
      void endSubRun(const art::SubRun& sr) override;
      void reconfigure(fhicl::ParameterSet const& pset) ;
      void analyze(art::Event const & evt) override;


    private:

      std::string fEnergyRecoNueLabel;
      std::string fEnergyRecoNumuLabel;

      unsigned int fRun,fSubrun,fEvent;

      double fFidVolXMax;
      double fFidVolYMax;
      double fFidVolZMin;
      double fFidVolZMax;

      double fEtrue; 

      double fErecoNue;
      double fRecoLepEnNue;
      double fRecoHadEnNue;
      int fRecoMethodNue; // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
      double fErecoNumu; 
      double fRecoLepEnNumu;
      double fRecoHadEnNumu;
      int fRecoMethodNumu; // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
      int fLongestTrackContNumu; // 1 = contained, 0 = exiting, -1 = not set
      int fTrackMomMethodNumu; // 1 = range, 0 = MCS, -1 = not set

      double fNuVtxX;
      double fNuVtxY;
      double fNuVtxZ;
      int fNuPDG;
      int fMode; // 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      int fCCNC; // 0=CC, 1=NC
      double fNuMomX;
      double fNuMomY;
      double fNuMomZ;
      double fNuMomT;

      int fLepPDG;
      double fLepMomX;
      double fLepMomY;
      double fLepMomZ;
      double fLepMomT;
      double fLepNuAngle;

      TH1D* fLepMomResNumuCont;
      TH1D* fHadEnResNumuCont;
      TH1D* fEnergyResNumuCont;
      TH1D* fLepMomResNumuExit;
      TH1D* fHadEnResNumuExit;
      TH1D* fEnergyResNumuExit;
      TH1D* fLepEnResNue;
      TH1D* fHadEnResNue;
      TH1D* fEnergyResNue; 

      std::vector<double> fNumuContBinEdge;
      std::vector<double> fNumuExitBinEdge;
      std::vector<double> fNueBinEdge;    

      std::vector<TH1D*> fEnResEnNumuCont;
      std::vector<TH1D*> fEnResEnNumuExit;
      std::vector<TH1D*> fEnResEnNue;

      std::vector<TF1*> fFitEnResNumuCont;
      std::vector<TF1*> fFitEnResNumuExit;
      std::vector<TF1*> fFitEnResNue;

      TGraphErrors* fGrResEnNumuCont;
      TGraphErrors* fGrSigEnNumuCont;
      TGraphErrors* fGrResEnNumuExit;
      TGraphErrors* fGrSigEnNumuExit;
      TGraphErrors* fGrResEnNue;
      TGraphErrors* fGrSigEnNue;

  }; // class CheckRecoEnergy

  //------------------------------------------------------------------------------
  CheckRecoEnergy::CheckRecoEnergy(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  dune::CheckRecoEnergy::~CheckRecoEnergy(){}

  //------------------------------------------------------------------------------
  void CheckRecoEnergy::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
    fEnergyRecoNumuLabel = pset.get<std::string>("EnergyRecoNumuLabel");

    fFidVolXMax  = pset.get<double>("FidVolXMax");
    fFidVolYMax  = pset.get<double>("FidVolYMax");
    fFidVolZMin  = pset.get<double>("FidVolZMin");
    fFidVolZMax  = pset.get<double>("FidVolZMax");

    fNumuContBinEdge = pset.get<std::vector<double>>("NumuContBinEdge");
    fNumuExitBinEdge = pset.get<std::vector<double>>("NumuExitBinEdge");
    fNueBinEdge = pset.get<std::vector<double>>("NueBinEdge");

    return;
  }

  //------------------------------------------------------------------------------
  void CheckRecoEnergy::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    fLepMomResNumuCont = tfs->make<TH1D>("LepMomResNumuCont", "", 80, -2.0, 2.0);
    fHadEnResNumuCont = tfs->make<TH1D>("HadEnResNumuCont", "", 40, -2.0, 2.0);
    fEnergyResNumuCont = tfs->make<TH1D>("EnergyResNumuCont", "", 40, -2.0, 2.0);
    fLepMomResNumuExit = tfs->make<TH1D>("LepMomResNumuExit", "", 80, -2.0, 2.0);
    fHadEnResNumuExit = tfs->make<TH1D>("HadEnResNumuExit", "", 40, -2.0, 2.0);
    fEnergyResNumuExit = tfs->make<TH1D>("EnergyResNumuExit", "", 40, -2.0, 2.0);
    fLepEnResNue = tfs->make<TH1D>("LepEnResNue", "", 80, -2.0, 2.0);
    fHadEnResNue = tfs->make<TH1D>("HadEnResNue", "", 40, -2.0, 2.0);
    fEnergyResNue = tfs->make<TH1D>("EnergyResNue", "", 40, -2.0, 2.0);

    for(unsigned int h=0; h<fNumuContBinEdge.size()-1;h++)
      {
        TH1D *hEnResEnNumuCont = tfs->make<TH1D>(Form("hEnResEnNumuCont%d",h), "", 40, -2.0, 2.0);
        fEnResEnNumuCont.push_back(hEnResEnNumuCont);
        TF1 *gausEnResNumuCont = tfs->make<TF1>(Form("gausEnResNumuCont%d",h), "gaus", -1.0, 1.0);
        fFitEnResNumuCont.push_back(gausEnResNumuCont);
      }

    for(unsigned int h=0; h<fNumuExitBinEdge.size()-1;h++)
      {
        TH1D *hEnResEnNumuExit = tfs->make<TH1D>(Form("hEnResEnNumuExit%d",h), "", 40, -2.0, 2.0);
        fEnResEnNumuExit.push_back(hEnResEnNumuExit);
        TF1 *gausEnResNumuExit = tfs->make<TF1>(Form("gausEnResNumuExit%d",h), "gaus", -1.0, 1.0);
        fFitEnResNumuExit.push_back(gausEnResNumuExit);
      }

    for(unsigned int h=0; h<fNueBinEdge.size()-1;h++)
      {
        TH1D *hEnResEnNue = tfs->make<TH1D>(Form("hEnResEnNue%d",h), "", 40, -2.0, 2.0);
        fEnResEnNue.push_back(hEnResEnNue);
        TF1 *gausEnResNue = tfs->make<TF1>(Form("gausEnResNue%d",h), "gaus", -1.0, 1.0);
        fFitEnResNue.push_back(gausEnResNue);
      }

    fGrResEnNumuCont = tfs->makeAndRegister<TGraphErrors>("fGrResEnNumuCont", "fGrResEnNumuCont", (int)(fNumuContBinEdge.size())-1);
    fGrSigEnNumuCont = tfs->makeAndRegister<TGraphErrors>("fGrSigEnNumuCont", "fGrSigEnNumuCont", (int)(fNumuContBinEdge.size())-1);

    fGrResEnNumuExit = tfs->makeAndRegister<TGraphErrors>("fGrResEnNumuExit", "fGrResEnNumuExit", (int)(fNumuExitBinEdge.size())-1);
    fGrSigEnNumuExit = tfs->makeAndRegister<TGraphErrors>("fGrSigEnNumuExit", "fGrSigEnNumuExit", (int)(fNumuExitBinEdge.size())-1);

    fGrResEnNue = tfs->makeAndRegister<TGraphErrors>("fGrResEnNue", "fGrResEnNue", (int)(fNueBinEdge.size())-1);
    fGrSigEnNue = tfs->makeAndRegister<TGraphErrors>("fGrSigEnNue", "fGrSigEnNue", (int)(fNueBinEdge.size())-1);

    return;
  }

  //------------------------------------------------------------------------------
  void CheckRecoEnergy::beginSubRun(const art::SubRun& sr){
  }

  //------------------------------------------------------------------------------
  void CheckRecoEnergy::analyze(art::Event const & evt)
  {
    art::Handle<dune::EnergyRecoOutput> ereconuein;
    if (!evt.getByLabel(fEnergyRecoNueLabel, ereconuein)) return;

    art::Handle<dune::EnergyRecoOutput> ereconumuin;
    if (!evt.getByLabel(fEnergyRecoNumuLabel, ereconumuin)) return;

    fRun = evt.id().run();
    fSubrun = evt.id().subRun();
    fEvent = evt.id().event();

    fErecoNue          = ereconuein->fNuLorentzVector.E();
    fRecoLepEnNue      = ereconuein->fLepLorentzVector.E();
    fRecoHadEnNue      = ereconuein->fHadLorentzVector.E();
    fRecoMethodNue     = ereconuein->recoMethodUsed;
    fErecoNumu         = ereconumuin->fNuLorentzVector.E();
    fRecoLepEnNumu     = ereconumuin->fLepLorentzVector.E();
    fRecoHadEnNumu     = ereconumuin->fHadLorentzVector.E();
    fRecoMethodNumu    = ereconumuin->recoMethodUsed;
    fLongestTrackContNumu  = ereconumuin->longestTrackContained;
    fTrackMomMethodNumu    = ereconumuin->trackMomMethod;

    art::Handle< std::vector<simb::MCTruth> > mct;
    std::vector< art::Ptr<simb::MCTruth> > truth;
    if( evt.getByLabel("generator", mct) )
      art::fill_ptr_vector(truth, mct);
    else
      mf::LogWarning("CheckRecoEnergy") << "No MCTruth.";

    for(size_t i=0; i<truth.size(); i++){

      if(i>1){
        mf::LogWarning("CheckRecoEnergy") << "Skipping MC truth index " << i;
        continue;
      }

      fCCNC     = truth[i]->GetNeutrino().CCNC();  //0=CC 1=NC
      fNuPDG    = truth[i]->GetNeutrino().Nu().PdgCode();
      fMode     = truth[i]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      fNuVtxX = truth[i]->GetNeutrino().Nu().Vx();
      fNuVtxY = truth[i]->GetNeutrino().Nu().Vy();
      fNuVtxZ = truth[i]->GetNeutrino().Nu().Vz();
      fEtrue    = truth[i]->GetNeutrino().Nu().E();
      fNuMomX   = truth[i]->GetNeutrino().Nu().Momentum().X();
      fNuMomY   = truth[i]->GetNeutrino().Nu().Momentum().Y();
      fNuMomZ   = truth[i]->GetNeutrino().Nu().Momentum().Z();
      fNuMomT   = truth[i]->GetNeutrino().Nu().Momentum().T();
      fLepPDG     = truth[i]->GetNeutrino().Lepton().PdgCode();
      fLepMomX    = truth[i]->GetNeutrino().Lepton().Momentum().X();
      fLepMomY    = truth[i]->GetNeutrino().Lepton().Momentum().Y();
      fLepMomZ    = truth[i]->GetNeutrino().Lepton().Momentum().Z();
      fLepMomT    = truth[i]->GetNeutrino().Lepton().Momentum().T();
      fLepNuAngle = truth[i]->GetNeutrino().Nu().Momentum().Vect().Angle(truth[i]->GetNeutrino().Lepton().Momentum().Vect());

      //true CC event with true vertex in fiducial volume
      if(fCCNC == 0 && fabs(fNuVtxX) < fFidVolXMax && fabs(fNuVtxY) < fFidVolYMax && fabs(fNuVtxZ) > fFidVolZMin && fabs(fNuVtxZ) < fFidVolZMax) 
	{
	  if(fNuPDG == 14) //true numu
	    {
	      if(fLongestTrackContNumu == 1)  //longest reco track contained
		{
		  fLepMomResNumuCont->Fill((fRecoLepEnNumu - fLepMomT) / fLepMomT);
		  fHadEnResNumuCont->Fill((fRecoHadEnNumu - (fEtrue - fLepMomT)) / (fEtrue - fLepMomT));
		  fEnergyResNumuCont->Fill((fErecoNumu - fEtrue) / fEtrue);
		  for(unsigned int e=0; e<fNumuContBinEdge.size()-1; e++)
		    {
		      if(fEtrue > fNumuContBinEdge.at(e) && fEtrue <= fNumuContBinEdge.at(e+1))
			fEnResEnNumuCont.at(e)->Fill((fErecoNumu - fEtrue) / fEtrue);
		    }
		}
	      else if(fLongestTrackContNumu == 0)  //longest reco track exiting
		{
                  fLepMomResNumuExit->Fill((fRecoLepEnNumu - fLepMomT) / fLepMomT);
                  fHadEnResNumuExit->Fill((fRecoHadEnNumu - (fEtrue - fLepMomT)) / (fEtrue - fLepMomT));
                  fEnergyResNumuExit->Fill((fErecoNumu - fEtrue) / fEtrue);
		  for(unsigned int e=0; e<fNumuExitBinEdge.size()-1; e++)
                    {
                      if(fEtrue > fNumuExitBinEdge.at(e) && fEtrue <= fNumuExitBinEdge.at(e+1))
                        fEnResEnNumuExit.at(e)->Fill((fErecoNumu - fEtrue) / fEtrue);
                    }
		}
	    }
	  else if(fNuPDG == 12) //true nue
	    {
	      if(fRecoMethodNue == 2) //at least one reco shower
		{
		  fLepEnResNue->Fill((fRecoLepEnNue - fLepMomT) / fLepMomT);
		  fHadEnResNue->Fill((fRecoHadEnNue - (fEtrue - fLepMomT)) / (fEtrue - fLepMomT));
		  fEnergyResNue->Fill((fErecoNue - fEtrue) / fEtrue);
		  for(unsigned int e=0; e<fNueBinEdge.size()-1; e++)
                    {
                      if(fEtrue > fNueBinEdge.at(e) && fEtrue <= fNueBinEdge.at(e+1))
                        fEnResEnNue.at(e)->Fill((fErecoNue - fEtrue) / fEtrue);
                    }
		}
	    }
	}
    } //end of for(size_t i=0; i<truth.size(); i++)

    return;
  }
  
  void CheckRecoEnergy::endJob()
  {
    for(unsigned int h=0; h<fNumuContBinEdge.size()-1; h++)
      fEnResEnNumuCont.at(h)->Fit(fFitEnResNumuCont.at(h));

    for(unsigned int h=0; h<fNumuExitBinEdge.size()-1; h++)
      fEnResEnNumuExit.at(h)->Fit(fFitEnResNumuExit.at(h));

    for(unsigned int h=0; h<fNueBinEdge.size()-1; h++)
      fEnResEnNue.at(h)->Fit(fFitEnResNue.at(h));

    for(unsigned int h=0; h<fNumuContBinEdge.size()-1; h++)
      {
      fGrResEnNumuCont->SetPoint(h, 0.5 * (fNumuContBinEdge.at(h) + fNumuContBinEdge.at(h+1)), fFitEnResNumuCont.at(h)->GetParameter(1));
      fGrResEnNumuCont->SetPointError(h, 0.5 * (fNumuContBinEdge.at(h+1) - fNumuContBinEdge.at(h)), 0.0);
      fGrSigEnNumuCont->SetPoint(h, 0.5 * (fNumuContBinEdge.at(h) + fNumuContBinEdge.at(h+1)), fFitEnResNumuCont.at(h)->GetParameter(2));
      fGrSigEnNumuCont->SetPointError(h, 0.5 * (fNumuContBinEdge.at(h+1) - fNumuContBinEdge.at(h)), 0.0);
      }

    for(unsigned int h=0; h<fNumuExitBinEdge.size()-1; h++)
      {
	fGrResEnNumuExit->SetPoint(h, 0.5 * (fNumuExitBinEdge.at(h) + fNumuExitBinEdge.at(h+1)), fFitEnResNumuExit.at(h)->GetParameter(1));
	fGrResEnNumuExit->SetPointError(h, 0.5 * (fNumuExitBinEdge.at(h+1) - fNumuExitBinEdge.at(h)), 0.0);
	fGrSigEnNumuExit->SetPoint(h, 0.5 * (fNumuExitBinEdge.at(h) + fNumuExitBinEdge.at(h+1)), fFitEnResNumuExit.at(h)->GetParameter(2));
	fGrSigEnNumuExit->SetPointError(h, 0.5 * (fNumuExitBinEdge.at(h+1) - fNumuExitBinEdge.at(h)), 0.0);
      }

    for(unsigned int h=0; h<fNueBinEdge.size()-1; h++)
      {
	fGrResEnNue->SetPoint(h, 0.5 * (fNueBinEdge.at(h) + fNueBinEdge.at(h+1)), fFitEnResNue.at(h)->GetParameter(1));
	fGrResEnNue->SetPointError(h, 0.5 * (fNueBinEdge.at(h+1) - fNueBinEdge.at(h)), 0.0);
	fGrSigEnNue->SetPoint(h, 0.5 * (fNueBinEdge.at(h) + fNueBinEdge.at(h+1)), fFitEnResNue.at(h)->GetParameter(2));
	fGrSigEnNue->SetPointError(h, 0.5 * (fNueBinEdge.at(h+1) - fNueBinEdge.at(h)), 0.0);
      }
  }
  
  //------------------------------------------------------------------------------
  void CheckRecoEnergy::endSubRun(const art::SubRun& sr){
  }

  DEFINE_ART_MODULE(CheckRecoEnergy)

} // namespace dune

#endif // CheckRecoEnergy_H
