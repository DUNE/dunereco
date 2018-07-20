////////////////////////////////////////////////////////////////////////
//
// \file EnergyReco_module.cc
//
// n.grant.3@warwick.ac.uk
//
///////////////////////////////////////////////////////////////////////

#ifndef EnergyReco_H
#define EnergyReco_H

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
#include "art/Persistency/Common/PtrMaker.h"

#include "nutools/NuReweight/art/NuReweight.h"
#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

#include "TTimeStamp.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"

namespace dune {

  class EnergyReco : public art::EDProducer {

    public:

      explicit EnergyReco(fhicl::ParameterSet const& pset);
      virtual ~EnergyReco();
      void beginJob() override;
      void beginSubRun(art::SubRun& sr) override;
      void endSubRun(art::SubRun& sr) override;
      void reconfigure(fhicl::ParameterSet const& pset) ;
      void produce(art::Event& evt) override;

    private:

      void  PrepareEvent(const art::Event& event);
      bool insideContVol(const double posX, const double posY, const double posZ);

      //Run information
      int run;
      int subrun;
      int event;
      float evttime;
      float taulife;
      short isdata;

      std::string fWireModuleLabel;
      std::string fHitsModuleLabel;
      std::string fTrackModuleLabel;
      std::string fShowerModuleLabel;

      art::ServiceHandle<geo::Geometry> fGeom;

      calo::CalorimetryAlg fCaloAlg;

      int fRecoMethod;
      double fContVolCut;
      double fWirecharge;
      double fMaxTrackLength;
      double fLongestTrackMCSMom;
      double fLongestTrackCharge;
      double fTotalEventCharge;
      bool fLongestTrackContained;
      double fMaxShowerCharge;
      art::Ptr<recob::Track> fBestTrack;
      art::Ptr<recob::Shower> fBestShower;

      double fGradTrkMomRange;
      double fIntTrkMomRange;
      double fGradTrkMomMCS;
      double fIntTrkMomMCS;
      double fGradNuMuHadEnCont;
      double fIntNuMuHadEnCont;
      double fGradNuMuHadEnExit;
      double fIntNuMuHadEnExit;
      double fGradShwEnergy;
      double fIntShwEnergy;
      double fGradNuEHadEn; 
      double fIntNuEHadEn;
      double fRecombFactor;
  }; // class EnergyReco


  //------------------------------------------------------------------------------
  EnergyReco::EnergyReco(fhicl::ParameterSet const& pset)
    : fCaloAlg (pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  {
    produces<dune::EnergyRecoOutput>();
    produces<art::Assns<dune::EnergyRecoOutput, recob::Track>>();
    produces<art::Assns<dune::EnergyRecoOutput, recob::Shower>>();
    this->reconfigure(pset);
  }

  dune::EnergyReco::~EnergyReco(){}

  //------------------------------------------------------------------------------
  void EnergyReco::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fRecoMethod = pset.get<int>("RecoMethod");

    fWireModuleLabel = pset.get< std::string >("WireModuleLabel");
    fHitsModuleLabel = pset.get< std::string >("HitsModuleLabel");
    fTrackModuleLabel = pset.get< std::string >("TrackModuleLabel");
    fShowerModuleLabel = pset.get< std::string >("ShowerModuleLabel");

    fGradTrkMomRange = pset.get<double>("GradTrkMomRange");
    fIntTrkMomRange = pset.get<double>("IntTrkMomRange");
    fGradTrkMomMCS = pset.get<double>("GradTrkMomMCS");
    fIntTrkMomMCS = pset.get<double>("IntTrkMomMCS");
    fGradNuMuHadEnCont = pset.get<double>("GradNuMuHadEnCont");
    fIntNuMuHadEnCont = pset.get<double>("IntNuMuHadEnCont");
    fGradNuMuHadEnExit = pset.get<double>("GradNuMuHadEnExit");
    fIntNuMuHadEnExit = pset.get<double>("IntNuMuHadEnExit");
    fGradShwEnergy = pset.get<double>("GradShwEnergy");
    fIntShwEnergy = pset.get<double>("IntShwEnergy");
    fGradNuEHadEn = pset.get<double>("GradNuEHadEn");
    fIntNuEHadEn = pset.get<double>("IntNuEHadEn");
    fContVolCut = pset.get<double>("ContVolCut");
    fRecombFactor = pset.get<double>("RecombFactor");
  }


  //------------------------------------------------------------------------------
  void EnergyReco::beginJob()
  {
  }

  //------------------------------------------------------------------------------
  void EnergyReco::beginSubRun(art::SubRun& sr)
  {
  }

  //------------------------------------------------------------------------------
  void EnergyReco::produce(art::Event& evt)
  {
    auto erecoout = std::make_unique<dune::EnergyRecoOutput>();
    auto assnstrk = std::make_unique<art::Assns<dune::EnergyRecoOutput, recob::Track>>();
    auto assnsshw = std::make_unique<art::Assns<dune::EnergyRecoOutput, recob::Shower>>();
    //art::PtrMaker<dune::EnergyRecoOutput> makeEnergyRecoOutputPtr(evt, *this);
    this->PrepareEvent(evt);

    double longestTrackMom = 0.0;
    double maxShowerEnergy = 0.0;
    double corrHadEnergy = 0.0;

    erecoout->fRecoVertex.SetX(0.0);
    erecoout->fRecoVertex.SetY(0.0);
    erecoout->fRecoVertex.SetZ(0.0);
    erecoout->recoMethodUsed = -1;
    erecoout->longestTrackContained = -1;
    erecoout->trackMomMethod = -1;


    if (fRecoMethod == 1){//split event into longest reco track and hadronic part
      erecoout->recoMethodUsed = 1;
      //at least one reco track, longest reco track is either contained or is exiting with a defined value of MCS track momentum
      if (fMaxTrackLength >= 0.0 && 
          (fLongestTrackContained || (!fLongestTrackContained && fLongestTrackMCSMom >= 0.0))){

        erecoout->fRecoVertex.SetX(fBestTrack->Start().X());
	erecoout->fRecoVertex.SetY(fBestTrack->Start().Y());
	erecoout->fRecoVertex.SetZ(fBestTrack->Start().Z());

        if(fLongestTrackContained)
	  {
	    erecoout->longestTrackContained = 1;
	    //Some contained tracks can have reconstructed lengths that are too short for their momentum
	    //Use MCS momentum if ratio of range / MCS is < 0.7 (even though track is contained)
	    if(fLongestTrackMCSMom >= 0 
	       && ((fMaxTrackLength - fIntTrkMomRange) / fGradTrkMomRange) / ((fLongestTrackMCSMom - fIntTrkMomMCS) / fGradTrkMomMCS) < 0.7)           
	      {
		longestTrackMom = (fLongestTrackMCSMom - fIntTrkMomMCS) / fGradTrkMomMCS;
		corrHadEnergy = ((fCaloAlg.ElectronsFromADCArea(fTotalEventCharge - fLongestTrackCharge, 2) * (1.0 / fRecombFactor) / util::kGeVToElectrons) - fIntNuMuHadEnExit) / fGradNuMuHadEnExit;
                erecoout->trackMomMethod = 0;
                
	      }
	    else
	      {
                longestTrackMom = (fMaxTrackLength - fIntTrkMomRange) / fGradTrkMomRange;
                corrHadEnergy = ((fCaloAlg.ElectronsFromADCArea(fTotalEventCharge - fLongestTrackCharge, 2) * (1.0 / fRecombFactor) / util::kGeVToElectrons) - fIntNuMuHadEnCont) / fGradNuMuHadEnCont;
		erecoout->trackMomMethod = 1;
	      }
	  }
        else
	  {
	   erecoout->longestTrackContained = 0;
           longestTrackMom = (fLongestTrackMCSMom - fIntTrkMomMCS) / fGradTrkMomMCS;
	   corrHadEnergy = ((fCaloAlg.ElectronsFromADCArea(fTotalEventCharge - fLongestTrackCharge, 2) * (1.0 / fRecombFactor) / util::kGeVToElectrons) - fIntNuMuHadEnExit) / fGradNuMuHadEnExit;
	   erecoout->trackMomMethod = 0;
	  }
        erecoout->fLepLorentzVector.SetPx(longestTrackMom*fBestTrack->VertexDirection().X());
        erecoout->fLepLorentzVector.SetPy(longestTrackMom*fBestTrack->VertexDirection().Y());
        erecoout->fLepLorentzVector.SetPz(longestTrackMom*fBestTrack->VertexDirection().Z());
        double Emu = std::sqrt(longestTrackMom*longestTrackMom+0.1056583745*0.1056583745);
        erecoout->fLepLorentzVector.SetE(Emu);

        erecoout->fHadLorentzVector.SetE(corrHadEnergy);
        erecoout->fNuLorentzVector.SetE(Emu + corrHadEnergy);
      }
      else{
	erecoout->recoMethodUsed = 3;
        erecoout->fNuLorentzVector.SetE(fCaloAlg.ElectronsFromADCArea(fWirecharge, 2)/fRecombFactor / util::kGeVToElectrons);
      }
    }
    else if (fRecoMethod == 2){//split event into reco shower with highest charge and hadronic part
      erecoout->recoMethodUsed = 2;
      //at least one reco shower
      if(fMaxShowerCharge >= 0.0){

	erecoout->fRecoVertex.SetX(fBestShower->ShowerStart().X());
        erecoout->fRecoVertex.SetY(fBestShower->ShowerStart().Y());
        erecoout->fRecoVertex.SetZ(fBestShower->ShowerStart().Z());

        maxShowerEnergy = ((fCaloAlg.ElectronsFromADCArea(fMaxShowerCharge, 2) * (1.0 / fRecombFactor) / util::kGeVToElectrons) - fIntShwEnergy) / fGradShwEnergy;
        corrHadEnergy = ((fCaloAlg.ElectronsFromADCArea(fTotalEventCharge - fMaxShowerCharge, 2) * (1.0 / fRecombFactor) / util::kGeVToElectrons) - fIntNuEHadEn) / fGradNuEHadEn;

        erecoout->fLepLorentzVector.SetPx(maxShowerEnergy*fBestShower->Direction().X());
        erecoout->fLepLorentzVector.SetPy(maxShowerEnergy*fBestShower->Direction().Y());
        erecoout->fLepLorentzVector.SetPz(maxShowerEnergy*fBestShower->Direction().Z());
        double Eel = std::sqrt(maxShowerEnergy*maxShowerEnergy+0.0005109989461*0.0005109989461);
        erecoout->fLepLorentzVector.SetE(Eel);
	erecoout->fHadLorentzVector.SetE(corrHadEnergy);
        erecoout->fNuLorentzVector.SetE(Eel + corrHadEnergy);
      }
      else{
	erecoout->recoMethodUsed = 3;
        erecoout->fNuLorentzVector.SetE(fCaloAlg.ElectronsFromADCArea(fWirecharge, 2)/fRecombFactor / util::kGeVToElectrons);
        //0.63: recombination factor, 1/4.966e-3: calorimetry constant to convert ADC to number of electrons, Wion = 23.6 eV
      }
    }
    else if (fRecoMethod == 3){//use charges of all hits and convert to energy
      erecoout->recoMethodUsed = 3;
      erecoout->fNuLorentzVector.SetE(fCaloAlg.ElectronsFromADCArea(fWirecharge, 2)/fRecombFactor / util::kGeVToElectrons);
    }

    art::ProductID const prodId = getProductID<dune::EnergyRecoOutput>();
    art::EDProductGetter const* prodGetter = evt.productGetter(prodId);
    art::Ptr<dune::EnergyRecoOutput> EnergyRecoOutputPtr{ prodId, 0U, prodGetter };
    //art::Ptr<dune::EnergyRecoOutput> EnergyRecoOutputPtr = makeEnergyRecoOutputPtr(0);
    if (fBestTrack.isAvailable()) assnstrk->addSingle(EnergyRecoOutputPtr, fBestTrack);
    if (fBestShower.isAvailable()) assnsshw->addSingle(EnergyRecoOutputPtr, fBestShower);
    evt.put(std::move(erecoout));
    evt.put(std::move(assnstrk));
    evt.put(std::move(assnsshw));
  }

  //------------------------------------------------------------------------------
  void EnergyReco::PrepareEvent(const art::Event& evt){

    auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();
    taulife = detprop->ElectronLifetime();
    isdata = evt.isRealData();

    double t0 = detprop->TriggerOffset();

    // Wires
    art::Handle< std::vector<recob::Wire>> wireListHandle;
    std::vector<art::Ptr<recob::Wire>> wirelist;
    if (evt.getByLabel(fWireModuleLabel, wireListHandle))
      art::fill_ptr_vector(wirelist, wireListHandle);

    // Hits
//    art::Handle< std::vector<recob::Hit> > hitListHandle;
//    std::vector<art::Ptr<recob::Hit> > hitlist;
//    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
//      art::fill_ptr_vector(hitlist, hitListHandle);
    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    auto const& hitlist = *hitListHandle;

    // Tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);

    // Showers
    art::Handle<std::vector<recob::Shower>> shwListHandle;
    std::vector<art::Ptr<recob::Shower>> shwlist;
    if (evt.getByLabel(fShowerModuleLabel,shwListHandle))
      art::fill_ptr_vector(shwlist, shwListHandle);

    // Associations
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::SpacePoint> fmhs(hitListHandle, evt, fTrackModuleLabel);

    fWirecharge = 0;
    for (size_t i = 0; i<wirelist.size(); ++i){
      if (fGeom->SignalType(wirelist[i]->Channel()) == geo::kCollection){
	const recob::Wire::RegionsOfInterest_t& signalROI = wirelist[i]->SignalROI();
	for(const auto& range : signalROI.get_ranges()){
	  const std::vector<float>& signal = range.data();
	  raw::TDCtick_t roiFirstBinTick = range.begin_index();
	  for (size_t j = 0; j<signal.size(); ++j){
	    fWirecharge += signal[j]*exp((j+roiFirstBinTick)*0.5/taulife);
	  }
	}
      }
    }

    int ntracks = tracklist.size();

    fTotalEventCharge = 0.0;

    for (recob::Hit const& hit: hitlist) {
      if (hit.WireID().Plane == 2){
        fTotalEventCharge += hit.Integral() * fCaloAlg.LifetimeCorrection(hit.PeakTime(), t0);
      }
    }

    fMaxTrackLength = -1.0;
    int iLongestTrack = -1;

    for (int i = 0; i < ntracks; ++i){
      if(tracklist[i]->Length() > fMaxTrackLength){
	fMaxTrackLength = tracklist[i]->Length();
	iLongestTrack = i;
        fBestTrack = tracklist[i];
      }
    }

    fLongestTrackCharge = 0.0;
    fLongestTrackMCSMom = -1.0;
    fLongestTrackContained = true;

    trkf::TrackMomentumCalculator TrkMomCalc;

    if(iLongestTrack >= 0 && iLongestTrack <= ntracks-1 ){
      if (fmth.isValid()){
	std::vector< art::Ptr<recob::Hit> > vhit = fmth.at(iLongestTrack);
	for (size_t h = 0; h < vhit.size(); ++h){
          if (vhit[h]->WireID().Plane == 2){
            fLongestTrackCharge += vhit[h]->Integral() * fCaloAlg.LifetimeCorrection(vhit[h]->PeakTime(), t0);
            std::vector<art::Ptr<recob::SpacePoint> > spts = fmhs.at(vhit[h].key());
            if (spts.size()){
              if (!insideContVol(spts[0]->XYZ()[0], spts[0]->XYZ()[1], spts[0]->XYZ()[2]))
                fLongestTrackContained = false;
            }
	  }
	}
      }
      fLongestTrackMCSMom = TrkMomCalc.GetMomentumMultiScatterChi2(tracklist[iLongestTrack]);
    }

    fMaxShowerCharge = -1.0;
    double showerCharge;

    if (shwListHandle.isValid()){
      art::FindManyP<recob::Hit> fmsh(shwListHandle, evt, fShowerModuleLabel);
      for (size_t i = 0; i<shwlist.size(); ++i){
	showerCharge = 0.0;
	if (fmsh.isValid()){
	  std::vector< art::Ptr<recob::Hit> > vhit = fmsh.at(i);
	  for (size_t h = 0; h < vhit.size(); ++h){
            if (vhit[h]->WireID().Plane == 2)
              showerCharge += vhit[h]->Integral() * fCaloAlg.LifetimeCorrection(vhit[h]->PeakTime(), t0);
          }
	}
	if(showerCharge > fMaxShowerCharge){
          fBestShower = shwlist[i];
	  fMaxShowerCharge = showerCharge;
        }
      }
    }
  }

  //------------------------------------------------------------------------------
  bool EnergyReco::insideContVol(const double posX, const double posY, const double posZ)
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
	    (dista > fContVolCut) && (distb > fContVolCut)) inside = true;
	//y
	dista = fabs(maxy - posY);
	distb = fabs(posY - miny);
	if (inside && (posY > miny) && (posY < maxy) &&
	    (dista > fContVolCut) && (distb > fContVolCut)) inside = true;
	else inside = false;
	//z
	dista = fabs(maxz - posZ);
	distb = fabs(posZ - minz);
	if (inside && (posZ > minz) && (posZ < maxz) &&
	    (dista > fContVolCut) && (distb > fContVolCut)) inside = true;
	else inside = false;
      }

    return inside;

  }

  //------------------------------------------------------------------------------
  void EnergyReco::endSubRun(art::SubRun& sr){
  }

  DEFINE_ART_MODULE(EnergyReco)

} // namespace dunemva

#endif // EnergyReco_H


