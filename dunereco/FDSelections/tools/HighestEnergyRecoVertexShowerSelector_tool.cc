#include "HighestEnergyRecoVertexShowerSelector.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

art::Ptr<recob::Shower> FDSelectionTools::HighestEnergyRecoVertexShowerSelector::SelectShower(art::Event const & evt){
  art::Ptr<recob::Shower> myshower;


  std::cout << "fShowerModuleLabel: " << fShowerModuleLabel << std::endl;

  //Need the associated hits for the showers
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  if (!evt.getByLabel(fShowerModuleLabel, showerListHandle)){
    return myshower;
  }
  //art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);


  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)){
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);
  }

  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);
  //Loop over the neutrinos
  if (nu_pfps.size() != 1){
    std::cout<<"HighestEnergyRecoVertexShowerSelector: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return myshower; //empty track
  }
  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];
  //Loop over each PFP, find the associated showers and then find which one is the highest energy
  //int i_energyist_shower = -1;
  double energy_energyist_shower = -999.;
  art::ServiceHandle<geo::Geometry> geom;

  for (int i_child = 0; i_child < nu_pfp->NumDaughters(); i_child++){
    //Use the PFParticle map to get the child PFPs
    int id_child = nu_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[id_child];
    //Now get the associated shower 
    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, fShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(child_pfp.key());
    if (pfp_shower_vector.size() > 1){ //Found a PFP with more than one shower matched.  Complain and exit
      std::cout<<"Number of associated showers to a PFP is greater than 1: " << pfp_shower_vector.size() << std::endl;
      continue;
    } 
    else if (pfp_shower_vector.size() == 0){ //Don't bother complaining of no shower was found.  It's either missing or its a shower.  No biggie here
      continue;
    }
    const art::Ptr<recob::Shower> shower = pfp_shower_vector[0];
    //const std::vector<art::Ptr<recob::Hit> > showerHits = fmhs.at(shower.key());
    std::map<int,double> showerEnergy;
    if (shower->Energy().size() > 0)
    {
        for (unsigned int plane = 0; plane < geom->MaxPlanes(); ++plane)
          showerEnergy[plane] = shower->Energy()[plane];
    }
    else 
        showerEnergy[2] = dune_ana::DUNEAnaShowerUtils::GetHits(shower, evt, "pandoraShower").size();
      //showerEnergy[plane] = fShowerEnergyAlg.ShowerEnergy(showerHits, plane);
    int best_plane = -1;
    double highest_energy_plane = 0;
    for (std::map<int,double>::const_iterator showerEnergyIt = showerEnergy.begin(); showerEnergyIt != showerEnergy.end(); ++showerEnergyIt) {
      if (showerEnergyIt->second > highest_energy_plane) {
        highest_energy_plane = showerEnergyIt->second;
        best_plane = showerEnergyIt->first;
      }
    }
    if (best_plane < 0){
      continue;
    }

    double current_energy = showerEnergy.at(best_plane);
    if (current_energy > energy_energyist_shower) {
      energy_energyist_shower = current_energy;
      myshower = shower;
    }
  }

  return myshower;


}
