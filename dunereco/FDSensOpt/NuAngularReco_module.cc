////////////////////////////////////////////////////////////////////////
//
// \file NuAngularReco_module.cc
//
// Generated at Mon Oct 30 08:54:19 2023 by Henrique Souza using cetskelgen
//
////////////////////////////////////////////////////////////////////////

#ifndef NuAngularReco_H
#define NuAngularReco_H

//ART
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
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
#include "dunereco/FDSensOpt/FDSensOptData/AngularRecoOutput.h"
#include "dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include <memory>

namespace dune {


  class NuAngularReco : public art::EDProducer {
    public:
      explicit NuAngularReco(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      NuAngularReco(NuAngularReco const&) = delete;
      NuAngularReco(NuAngularReco&&) = delete;
      NuAngularReco& operator=(NuAngularReco const&) = delete;
      NuAngularReco& operator=(NuAngularReco&&) = delete;

      // Required functions.
      void produce(art::Event& evt) override;

      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

    private:
            art::Ptr<recob::Track> GetLongestTrack(const art::Event& event);
            art::Ptr<recob::Shower> GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                           detinfo::DetectorPropertiesData const& detProp,
                                                           const art::Event& event);
      std::string fWireLabel;
      std::string fHitLabel;
      std::string fTrackLabel;
      std::string fShowerLabel;
      std::string fTrackToHitLabel;
      std::string fShowerToHitLabel;
      std::string fHitToSpacePointLabel;

      int fRecoMethod;

      NeutrinoAngularRecoAlg fNeutrinoAngularRecoAlg;
  }; // class NuAngularReco


//-----------------------------------------------------------------------------------------------------------------------------------------

NuAngularReco::NuAngularReco(fhicl::ParameterSet const& pset) :
  EDProducer(pset),
    fWireLabel(pset.get<std::string>("WireLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fTrackLabel(pset.get<std::string>("TrackLabel")),
    fShowerLabel(pset.get<std::string>("ShowerLabel")),
    fTrackToHitLabel(pset.get<std::string>("TrackToHitLabel")),
    fShowerToHitLabel(pset.get<std::string>("ShowerToHitLabel")),
    fHitToSpacePointLabel(pset.get<std::string>("HitToSpacePointLabel")),
    fRecoMethod(pset.get<int>("RecoMethod")),
    fNeutrinoAngularRecoAlg(pset.get<fhicl::ParameterSet>("NeutrinoAngularRecoAlg"),fTrackLabel,fShowerLabel,
        fHitLabel,fWireLabel,fTrackToHitLabel,fShowerToHitLabel,fHitToSpacePointLabel)
{
  produces<dune::AngularRecoOutput>();
  produces<art::Assns<dune::AngularRecoOutput, recob::Track>>();
  produces<art::Assns<dune::AngularRecoOutput, recob::Shower>>();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void NuAngularReco::produce(art::Event& evt)
{
  std::unique_ptr<dune::AngularRecoOutput> angularRecoOutput;
  auto assnstrk = std::make_unique<art::Assns<dune::AngularRecoOutput, recob::Track>>();
  auto assnsshw = std::make_unique<art::Assns<dune::AngularRecoOutput, recob::Shower>>();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  art::Ptr<recob::Track> longestTrack(this->GetLongestTrack(evt));
  art::Ptr<recob::Shower> highestChargeShower(this->GetHighestChargeShower(clockData, detProp, evt));

  if (fRecoMethod == 1)
      angularRecoOutput = std::make_unique<dune::AngularRecoOutput>(fNeutrinoAngularRecoAlg.CalculateNeutrinoAngle(longestTrack, evt));
  else if (fRecoMethod == 2)
      angularRecoOutput = std::make_unique<dune::AngularRecoOutput>(fNeutrinoAngularRecoAlg.CalculateNeutrinoAngle(highestChargeShower, evt));

  art::ProductID const prodId = evt.getProductID<dune::AngularRecoOutput>();
  art::EDProductGetter const* prodGetter = evt.productGetter(prodId);
  art::Ptr<dune::AngularRecoOutput> AngularRecoOutputPtr{ prodId, 0U, prodGetter };
  if (longestTrack.isAvailable()) assnstrk->addSingle(AngularRecoOutputPtr, longestTrack);
  if (highestChargeShower.isAvailable()) assnsshw->addSingle(AngularRecoOutputPtr, highestChargeShower);
  evt.put(std::move(angularRecoOutput));
  evt.put(std::move(assnstrk));
  evt.put(std::move(assnsshw));
}


//-----------------------------------------------------------------------------------------------------------------------------------------

void NuAngularReco::beginJob()
{
  // Implementation of optional member function here.
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void NuAngularReco::endJob()
{
  // Implementation of optional member function here.
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Track> NuAngularReco::GetLongestTrack(const art::Event &event)
{
    art::Ptr<recob::Track> pTrack{};
    const std::vector<art::Ptr<recob::Track> > tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
    art::FindManyP<recob::PFParticle> fmPFParticle(tracks, event, fTrackLabel);
    
    if (0 == tracks.size())
        return pTrack;

    double longestLength(std::numeric_limits<double>::lowest());
    for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
    {
        if (fmPFParticle.isValid()){
            std::vector<art::Ptr<recob::PFParticle>> pfp = fmPFParticle.at(iTrack);
            if (!lar_pandora::LArPandoraHelper::IsTrack(pfp[0]))
                continue;
        }
        const double length(tracks[iTrack]->Length());
        if (length-longestLength > std::numeric_limits<double>::epsilon())
        {
            longestLength = length;
            pTrack = tracks[iTrack];
        }
    }
    return pTrack;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Shower> NuAngularReco::GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                           detinfo::DetectorPropertiesData const& detProp,
                                                           const art::Event &event)
{
    art::Ptr<recob::Shower> pShower{};
    const std::vector<art::Ptr<recob::Shower> > showers(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));
    if (0 == showers.size())
        return pShower;

    double maxCharge(std::numeric_limits<double>::lowest());
    for (unsigned int iShower = 0; iShower < showers.size(); ++iShower)
    {
        const std::vector<art::Ptr<recob::Hit> > showerHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(showers[iShower],
            event,fShowerToHitLabel),2));
        const double showerCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, showerHits));
        if (showerCharge-maxCharge > std::numeric_limits<double>::epsilon())
        {
            maxCharge = showerCharge;
            pShower = showers[iShower];
        }
    }
    return pShower;
}

DEFINE_ART_MODULE(NuAngularReco)

} // namespace dune

#endif // AngularReco_H