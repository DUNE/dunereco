#include "HighestPandizzleScoreRecoVertexTrackSelector.h"

FDSelectionTools::HighestPandizzleScoreRecoVertexTrackSelector::HighestPandizzleScoreRecoVertexTrackSelector(fhicl::ParameterSet const& ps) :
    fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fPandizzleWeightFileName(ps.get< std::string> ("PandizzleWeightFileName")),
    fPandizzleAlg(ps),
    fPandizzleReader("")
{
  fPandizzleReader.AddVariable("PFPMichelNHits",&fTMVAPFPMichelNHits);
  fPandizzleReader.AddVariable("PFPMichelElectronMVA",&fTMVAPFPMichelElectronMVA);
  fPandizzleReader.AddVariable("PFPMichelRecoEnergyPlane2",&fTMVAPFPMichelRecoEnergyPlane2);
  fPandizzleReader.AddVariable("PFPTrackDeflecAngleSD",&fTMVAPFPTrackDeflecAngleSD);
  fPandizzleReader.AddVariable("PFPTrackLength",&fTMVAPFPTrackLength);
  fPandizzleReader.AddVariable("PFPTrackEvalRatio",&fTMVAPFPTrackEvalRatio);
  fPandizzleReader.AddVariable("PFPTrackConcentration",&fTMVAPFPTrackConcentration);
  fPandizzleReader.AddVariable("PFPTrackCoreHaloRatio",&fTMVAPFPTrackCoreHaloRatio);
  fPandizzleReader.AddVariable("PFPTrackConicalness",&fTMVAPFPTrackConicalness);
  fPandizzleReader.AddVariable("PFPTrackdEdxStart",&fTMVAPFPTrackdEdxStart);
  fPandizzleReader.AddVariable("PFPTrackdEdxEnd",&fTMVAPFPTrackdEdxEnd);
  fPandizzleReader.AddVariable("PFPTrackdEdxEndRatio",&fTMVAPFPTrackdEdxEndRatio);
  std::string weight_file_name = fPandizzleWeightFileName;
  std::string weight_file_path;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(weight_file_name, weight_file_path);
  fPandizzleReader.BookMVA("BDTG",weight_file_path);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Track> FDSelectionTools::HighestPandizzleScoreRecoVertexTrackSelector::SelectTrack(art::Event const & evt){
    /*
    std::cout << "selector fTrackModuleLabel: " << fTrackModuleLabel << std::endl;
    std::cout << "selector fPFParticleModuleLabel: " << fPFParticleModuleLabel << std::endl;
    std::cout << "selector fPandizzleWeightFileName: " << fPandizzleWeightFileName << std::endl; 
    */
  art::Ptr<recob::Track> selTrack;

  // Get PFParticles from event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;

  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  // Build PFParticle map
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  // Select Neutrino PFParticles
  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

  // Loop over the neutrinos
  if (nu_pfps.size() != 1)
  {
    std::cout << "HighestPandizzleScoreRecoVertexTrackSelector: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return selTrack; //empty track
  }

  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  // Loop over primaries, get the associated track and then find that with the highest pandizzle score
  double highestPandizzleScore(std::numeric_limits<double>::lowest());

  for (int i_child = 0; i_child < nu_pfp->NumDaughters(); i_child++)
  {
    // Use the PFParticle map to get the child PFPs
    int id_child = nu_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[id_child];

    // Get the associated track 
    art::FindManyP<recob::Track> fmtpfp(pfparticleListHandle, evt, fTrackModuleLabel);
    const std::vector<art::Ptr<recob::Track> > pfp_track_vector = fmtpfp.at(child_pfp.key());

    if (pfp_track_vector.size() > 1)
    { 
      // Found a PFP with more than one track matched.  Complain and move on
      std::cout<< "Number of associated tracks to a PFP is greater than 1: " << pfp_track_vector.size() << std::endl;
      continue; // empty
    } 
    else if (pfp_track_vector.size() == 0)
    { 
      // Don't bother complaining of no track was found.  It's either missing or its a shower.  No biggie here
      continue; 
    }
    const art::Ptr<recob::Track> track = pfp_track_vector[0];

    // Get Pandizzle score
    fPandizzleAlg.ProcessPFParticle(child_pfp, evt);

    fTMVAPFPMichelNHits = (float)(fPandizzleAlg.GetIntVar("PFPMichelNHits"));
    fTMVAPFPMichelElectronMVA = fPandizzleAlg.GetFloatVar("PFPMichelElectronMVA");
    fTMVAPFPMichelRecoEnergyPlane2 = fPandizzleAlg.GetFloatVar("PFPMichelRecoEnergyPlane2");
    fTMVAPFPTrackDeflecAngleSD = fPandizzleAlg.GetFloatVar("PFPTrackDeflecAngleSD");
    fTMVAPFPTrackLength = fPandizzleAlg.GetFloatVar("PFPTrackLength");
    fTMVAPFPTrackEvalRatio = fPandizzleAlg.GetFloatVar("PFPTrackEvalRatio");
    fTMVAPFPTrackConcentration = fPandizzleAlg.GetFloatVar("PFPTrackConcentration");
    fTMVAPFPTrackCoreHaloRatio = fPandizzleAlg.GetFloatVar("PFPTrackCoreHaloRatio");
    fTMVAPFPTrackConicalness = fPandizzleAlg.GetFloatVar("PFPTrackConicalness");
    fTMVAPFPTrackdEdxStart = fPandizzleAlg.GetFloatVar("PFPTrackdEdxStart");
    fTMVAPFPTrackdEdxEnd = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEnd");
    fTMVAPFPTrackdEdxEndRatio = fPandizzleAlg.GetFloatVar("PFPTrackdEdxEndRatio");
    fTMVAPFPTrackPIDA = fPandizzleAlg.GetFloatVar("PFPTrackPIDA");
    /*
    std::cout << "fTMVAPFPMichelNHits: " << fTMVAPFPMichelNHits << std::endl;
    std::cout << "fTMVAPFPMichelElectronMVA: " << fTMVAPFPMichelElectronMVA << std::endl;
    std::cout << "fTMVAPFPMichelRecoEnergyPlane2: " << fTMVAPFPMichelRecoEnergyPlane2 << std::endl;
    std::cout << "fTMVAPFPTrackDeflecAngleSD: " << fTMVAPFPTrackDeflecAngleSD << std::endl;
    std::cout << "fTMVAPFPTrackLength" << fTMVAPFPTrackLength << std::endl;
    std::cout << "fTMVAPFPTrackEvalRatio: " << fTMVAPFPTrackEvalRatio << std::endl;
    std::cout << "fTMVAPFPTrackConcentration: " << fTMVAPFPTrackConcentration << std::endl;
    std::cout << "fTMVAPFPTrackCoreHaloRatio: " << fTMVAPFPTrackCoreHaloRatio << std::endl;
    std::cout << "fTMVAPFPTrackConicalness: " << fTMVAPFPTrackConicalness << std::endl;
    std::cout << "fTMVAPFPTrackdEdxStart: " << fTMVAPFPTrackdEdxStart << std::endl;
    std::cout << "fTMVAPFPTrackdEdxEnd: " << fTMVAPFPTrackdEdxEnd << std::endl;
    std::cout << "fTMVAPFPTrackdEdxEndRatio: " << fTMVAPFPTrackdEdxEndRatio << std::endl;
    std::cout << "fTMVAPFPTrackPIDA: " << fTMVAPFPTrackPIDA << std::endl;
    */
    const double pandizzleScore = fPandizzleReader.EvaluateMVA("BDTG");

    //std::cout << "pandizzle score: " << pandizzleScore << std::endl;

    if (pandizzleScore > highestPandizzleScore)
    {
        highestPandizzleScore = pandizzleScore;
        selTrack = track;
    }
  }

  //std::cout << "highest score: " << highestPandizzleScore << std::endl;

  return selTrack;
}
