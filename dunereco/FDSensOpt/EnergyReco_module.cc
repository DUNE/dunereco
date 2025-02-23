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
#pragma GCC diagnostic ignored "-Wunused-variable"

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

            int fRecoMethod;
            int fLongestTrackMethod;



            NeutrinoEnergyRecoAlg fNeutrinoEnergyRecoAlg;
            ParticleSelectionAlg fParticleSelectionAlg;
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
    fRecoMethod(pset.get<int>("RecoMethod")),
    fLongestTrackMethod(pset.get<int>("LongestTrackMethod")),
    fNeutrinoEnergyRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"),fTrackLabel,fShowerLabel,
        fHitLabel,fWireLabel,fTrackToHitLabel,fShowerToHitLabel,fHitToSpacePointLabel),
    fParticleSelectionAlg(pset.get<fhicl::ParameterSet>("ParticleSelectionAlg"),fTrackLabel,fShowerLabel,
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
    art::Ptr<recob::Track> longestTrack(fParticleSelectionAlg.GetLongestTrack(evt, true));
    art::Ptr<recob::Shower> highestChargeShower(fParticleSelectionAlg.GetHighestChargeShower(clockData, detProp, evt));

    if (fRecoMethod == 1)
    {
        if (fLongestTrackMethod == 0 || !longestTrack.isAvailable() || longestTrack.isNull())
            energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(longestTrack, evt));
        else if (fLongestTrackMethod == 1)
            energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergyViaMuonRanging(longestTrack, evt));
        else if (fLongestTrackMethod == 2)
        {
            energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergyViaMuonMCS(longestTrack, evt));
        }
    }
    else if (fRecoMethod == 2)
        energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(highestChargeShower, evt));
    else if (fRecoMethod == 3)
        energyRecoOutput = std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(evt));

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
