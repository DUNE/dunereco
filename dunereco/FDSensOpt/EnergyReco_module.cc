////////////////////////////////////////////////////////////////////////
//
// \file EnergyReco_module.cc
//
// updated by Dom Brailsford (d.brailsford@lancaster.ac.uk)
// original module by Nick Grant (n.grant.3@warwick.ac.uk)
//
///////////////////////////////////////////////////////////////////////

#ifndef EnergyReco_H
#define EnergyReco_H

//STL
#include <limits>
//ROOT
//ART
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//DUNE
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/FDSensOpt/ParticleSelectionAlg/ParticleSelectionAlg.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

namespace dune {

    class EnergyReco : public art::EDProducer {

        public:

            explicit EnergyReco(fhicl::ParameterSet const& pset);
            void produce(art::Event& evt) override;

        private:
            std::string fWireLabel;
            std::string fHitLabel;
            std::string fTrackLabel;
            std::string fShowerLabel;
            std::string fTrackToHitLabel;
            std::string fShowerToHitLabel;
            std::string fHitToSpacePointLabel;
            std::string fPFParticleLabel;

            int fRecoMethod;
            int fLongestTrackMethod;
            bool fApplyMuFilters;
            bool fApplyElectronFilters;

            NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
            ParticleSelectionAlg fParticleSelectionAlg;

            ///The energy reconstruction method
            enum EnergyRecoMethod
            {
                kMuonAndHadronic = 1,                                 ///< muon momentum and hadronic deposited energy method
                kElectronAndHadronic,                                 ///< electron deposited energy and hadronic deposited energy method
                kAllCharges,                                          ///< summed wire charges
                kMuonProtonPionHadronic,                              ///< muon, protons and pions momenta + all other hadronic deposited energy (PID)
                kElectronProtonPionHadronic,                          ///< protons and pions momenta + electron and all other hadronic deposited energy (PID)
                kProtonPionHadronic,                                  ///< protons and pions momenta + summed wire charges
            };

            ///The muon momentum reconstruction method
            enum LongestTrackMethod
            {
                kTrackMethodNotSet = 0,                               ///< method not set
                kbyRange,                                             ///< muon momentum by range
                kbyMCS                                                ///< muon momentum by multiple scattering
            };

    }; // class EnergyReco

//-----------------------------------------------------------------------------------------------------------------------------------------

EnergyReco::EnergyReco(fhicl::ParameterSet const& pset) :
    EDProducer(pset),
    fWireLabel(pset.get<std::string>("WireLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fTrackLabel(pset.get<std::string>("TrackLabel")),
    fShowerLabel(pset.get<std::string>("ShowerLabel")),
    fTrackToHitLabel(pset.get<std::string>("TrackToHitLabel")),
    fShowerToHitLabel(pset.get<std::string>("ShowerToHitLabel")),
    fHitToSpacePointLabel(pset.get<std::string>("HitToSpacePointLabel")),
    fPFParticleLabel(pset.get<std::string>("PFParticleLabel")),
    fRecoMethod(pset.get<int>("RecoMethod")),
    fLongestTrackMethod(pset.get<int>("LongestTrackMethod")),
    fApplyMuFilters(pset.get<bool>("ApplyMuFilters")),
    fApplyElectronFilters(pset.get<bool>("ApplyElectronFilters")),
    fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"),fTrackLabel,fShowerLabel,
        fHitLabel,fWireLabel,fTrackToHitLabel,fShowerToHitLabel,fHitToSpacePointLabel),
    fParticleSelectionAlg(pset.get<fhicl::ParameterSet>("ParticleSelectionAlg"),fPFParticleLabel,fTrackLabel,fShowerLabel,
        fHitLabel,fTrackToHitLabel,fShowerToHitLabel,fHitToSpacePointLabel)
{
    produces<dune::EnergyRecoOutput>();
    produces<art::Assns<dune::EnergyRecoOutput, recob::Track>>();
    produces<art::Assns<dune::EnergyRecoOutput, recob::Shower>>();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void EnergyReco::produce(art::Event& evt)
{
    std::unique_ptr<dune::EnergyRecoOutput> energyRecoOutput;
    auto assnstrk = std::make_unique<art::Assns<dune::EnergyRecoOutput, recob::Track>>();
    auto assnsshw = std::make_unique<art::Assns<dune::EnergyRecoOutput, recob::Shower>>();


    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    art::Ptr<recob::Track> longestTrack(fParticleSelectionAlg.GetLongestTrackPID(evt, fApplyMuFilters));
    art::Ptr<recob::Shower> highestChargeShower(fParticleSelectionAlg.GetHighestChargeShowerPID(evt, fApplyElectronFilters));

    if (fRecoMethod == EnergyRecoMethod::kMuonAndHadronic)
        energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(longestTrack, evt, fLongestTrackMethod));
    else if (fRecoMethod == EnergyRecoMethod::kElectronAndHadronic)
        energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(highestChargeShower, evt));
    else if (fRecoMethod == EnergyRecoMethod::kAllCharges)
        energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(evt));
    else
    {
        std::vector<art::Ptr<recob::Track>> protonTracks(fParticleSelectionAlg.GenProtonCandidates(evt));
        std::vector<art::Ptr<recob::Track>> pionTracks(fParticleSelectionAlg.GenPionCandidates(evt));
        if (fRecoMethod == EnergyRecoMethod::kMuonProtonPionHadronic)
            energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergyPID(longestTrack, evt, fLongestTrackMethod, protonTracks, pionTracks));
        else if (fRecoMethod == EnergyRecoMethod::kElectronProtonPionHadronic)
            energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergyPID(highestChargeShower, evt, protonTracks, pionTracks));
        else if (fRecoMethod == EnergyRecoMethod::kProtonPionHadronic)
            energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergyPID(evt, protonTracks, pionTracks));
    }

    art::ProductID const prodId = evt.getProductID<dune::EnergyRecoOutput>();
    art::EDProductGetter const* prodGetter = evt.productGetter(prodId);
    art::Ptr<dune::EnergyRecoOutput> EnergyRecoOutputPtr{ prodId, 0U, prodGetter };
    if (longestTrack.isAvailable()) assnstrk->addSingle(EnergyRecoOutputPtr, longestTrack);
    if (highestChargeShower.isAvailable()) assnsshw->addSingle(EnergyRecoOutputPtr, highestChargeShower);
    evt.put(std::move(energyRecoOutput));
    evt.put(std::move(assnstrk));
    evt.put(std::move(assnsshw));
  }



DEFINE_ART_MODULE(EnergyReco)

} // namespace dune

#endif // EnergyReco_H
