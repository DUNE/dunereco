#include "HighestEnergyRecoVertexShowerSelector.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


FDSelectionTools::HighestEnergyRecoVertexShowerSelector::HighestEnergyRecoVertexShowerSelector(fhicl::ParameterSet const& ps) :
    fShowerModuleLabel(ps.get< std::string>("ModuleLabels.ShowerModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string>("ModuleLabels.PFParticleModuleLabel")),
    fShowerEnergyAlg(ps.get<fhicl::ParameterSet>("ShowerEnergyAlg"))
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Shower> FDSelectionTools::HighestEnergyRecoVertexShowerSelector::SelectShower(art::Event const & evt)
{
  art::Ptr<recob::Shower> selShower;

  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return selShower;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel);

  art::ServiceHandle<geo::Geometry> geom;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  double highestEnergy = -999.0;

  for (art::Ptr<recob::PFParticle> childPFP : nuChildren) 
  {
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(childPFP, evt, fPFParticleModuleLabel, fShowerModuleLabel))
      continue;

    art::Ptr<recob::Shower> childShower = dune_ana::DUNEAnaPFParticleUtils::GetShower(childPFP, evt, fPFParticleModuleLabel, fShowerModuleLabel);

    std::map<int,double> showerEnergy;

    if (childShower->Energy().size() > 0)
    {
        for (unsigned int plane = 0; plane < geom->MaxPlanes(); ++plane)
          showerEnergy[plane] = childShower->Energy()[plane];
    }
    else 
    {
      std::vector<art::Ptr<recob::Hit>> childHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(childPFP, evt, fPFParticleModuleLabel);

      for (unsigned int plane = 0; plane < geom->MaxPlanes(); ++plane)
        showerEnergy[plane] = fShowerEnergyAlg.ShowerEnergy(clockData, detProp, childHits, plane);
    }

    int bestPlane = -1;

    for (auto &entry : showerEnergy)
    {
      if (entry.second > highestEnergy)
      {
        highestEnergy = entry.second;
        bestPlane = entry.first;
      }
    }

    if (bestPlane < 0)
      continue;

    double currentEnergy = showerEnergy.at(bestPlane);
    if (currentEnergy > highestEnergy) 
    {
      highestEnergy = currentEnergy;
      selShower = childShower;
    }
  }

  return selShower;
}
