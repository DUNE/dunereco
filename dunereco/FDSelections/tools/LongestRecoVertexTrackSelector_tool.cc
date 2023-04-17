#include "LongestRecoVertexTrackSelector.h"

art::Ptr<recob::Track> FDSelectionTools::LongestRecoVertexTrackSelector::SelectTrack(art::Event const & evt){
  art::Ptr<recob::Track> mytrack;

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
    std::cout<<"LongestRecoVertexTrackSelector: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return mytrack; //empty track
  }
  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];
  //We are going to loop over the PFP children, get the associated track and then find the longest one.  So, default initialise the track length tracker to something grossly negative
  double longest_track_length = -999;
  //Loop over the children of the neutrino
  for (int i_child = 0; i_child < nu_pfp->NumDaughters(); i_child++){
    //Use the PFParticle map to get the child PFPs
    int id_child = nu_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[id_child];
    //Now get the associated track 
    art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fTrackModuleLabel);
    const std::vector<art::Ptr<recob::Track> > pfp_track_vector = fmtpfp.at(child_pfp.key());
    if (pfp_track_vector.size() > 1){ //Found a PFP with more than one track matched.  Complain and move on
      std::cout<<"Number of associated tracks to a PFP is greater than 1: " << pfp_track_vector.size() << std::endl;
      continue; //empty
    } 
    else if (pfp_track_vector.size() == 0){ //Don't bother complaining of no track was found.  It's either missing or its a shower.  No biggie here
      continue; 
    }
    const art::Ptr<recob::Track> track = pfp_track_vector[0];
    double current_track_length = track->Length();
    if (current_track_length > longest_track_length){
      longest_track_length = current_track_length;
      mytrack = track;
    }
  }
  return mytrack;
}
