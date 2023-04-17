#include "HighestPandrizzleScoreRecoVertexShowerSelector.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector::HighestPandrizzleScoreRecoVertexShowerSelector(fhicl::ParameterSet const& ps) :
    fShowerModuleLabel(ps.get< std::string> ("ModuleLabels.ShowerModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fCheatCharacterisation(ps.get<bool>("CheatCharacterisation", false)),
    fPFParticleHitCut(ps.get<unsigned int>("PFParticleHitCut", 0)),
    fDemandBDTScore(ps.get<bool>("DemandBDTScore", false)),
    fPandrizzleAlg(ps),
    fRecoNuVtxX(-9999),
    fRecoNuVtxY(-9999),
    fRecoNuVtxZ(-9999)
{
    if (fCheatCharacterisation)
    {
        //std::cout << "HighestPandrizzleScoreRecoVertexShowerSelector - fCheatCharacterisation is true" << std::endl; 

        fCheatShowerModuleLabel = ps.get< std::string >("ModuleLabels.CheatShowerModuleLabel");
        fNuGenModuleLabel = ps.get<std::string>("ModuleLabels.NuGenModuleLabel");
        fClusterModuleLabel = ps.get<std::string>("ModuleLabels.ClusterModuleLabel");
        fShowerPDGToCheat = ps.get<std::vector<int>> ("ShowerPDGToCheat");
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Shower> FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector::SelectShower(art::Event const & evt){

  art::Ptr<recob::Shower> selShower;
  art::Ptr<recob::Shower> backupSelShower;

  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
    art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  lar_pandora::PFParticleVector nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

  // Check we have one neutrino...
  if (nu_pfps.size() != 1)
  {
    std::cout << "HighestEnergyRecoVertexShowerSelector: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
    return selShower;
  }

  art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

  // Get neutrino vertex for pandrizzle MVA
  art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Vertex> > sel_pfp_vertices = fmvpfp.at(nu_pfp.key());

  if (sel_pfp_vertices.size() == 0)
  {
    return selShower;
  }
  else if (sel_pfp_vertices.size() > 1)
  {
    std::cout<< "CCNuSelection::FillVertexInformation Number of matched vertices bigger than 1: " << sel_pfp_vertices.size() << std::endl;
  }

  //always take the first vertex, even if there's more than one
  art::Ptr<recob::Vertex> matched_vertex = sel_pfp_vertices[0];
  fRecoNuVtxX = matched_vertex->position().X();
  fRecoNuVtxY = matched_vertex->position().Y();
  fRecoNuVtxZ = matched_vertex->position().Z();

  //Loop over each PFP, find the associated showers and then find which one is the highest energy
  double highestJamPandrizzleScore = std::numeric_limits<double>::lowest();
  double highestPandrizzleScore = std::numeric_limits<double>::lowest();
  bool found = false;

  for (int i_child = 0; i_child < nu_pfp->NumDaughters(); i_child++)
  {
    //Use the PFParticle map to get the child PFPs
    int id_child = nu_pfp->Daughter(i_child);
    art::Ptr<recob::PFParticle> child_pfp = pfparticleMap[id_child];

    /////////////////////
    bool cheatCharacterisation(false);

    if (fCheatCharacterisation && !fShowerPDGToCheat.empty())
    {
        simb::MCParticle* matched_mcparticle(nullptr);
        std::vector<art::Ptr<recob::Hit> > pfp_hits = GetPFPHits(child_pfp, evt);
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
        int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfp_hits, 1);
        if (g4id > 0)
        {
            art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
            matched_mcparticle = pi_serv->ParticleList().at(g4id);
        }

        if (matched_mcparticle)
        {
            const int absPdg(std::abs(matched_mcparticle->PdgCode()));
            //std::cout << "HighestEnergyRecoVertexShowerSelector - absPdg: " << absPdg << std::endl; 

            for (int cheatPDG : fShowerPDGToCheat)
            {
                const int absCheatPdg(std::abs(cheatPDG));

                if (absPdg != absCheatPdg)
                    continue;

                if ((absCheatPdg == 11 || absCheatPdg == 13))
                    continue;

                cheatCharacterisation = true;
                break;
            }

            // Get the generator record
            art::Handle<std::vector<simb::MCTruth> > mcTruthListHandle;
            std::vector<art::Ptr<simb::MCTruth> > mcList;

            if (evt.getByLabel(fNuGenModuleLabel, mcTruthListHandle))
                art::fill_ptr_vector(mcList, mcTruthListHandle);

            if (mcList.size() == 1)
            {
                const bool isNC = mcList[0]->GetNeutrino().CCNC();
                const bool isNue = (std::abs(mcList[0]->GetNeutrino().Nu().PdgCode()) == 12);
                const bool isNumu = (std::abs(mcList[0]->GetNeutrino().Nu().PdgCode()) == 14);
                const bool isCCNue = (isNue && !isNC);
                const bool isCCNumu = (isNumu && !isNC);

                //std::cout << "isNC ? " << (isNC ? "yes" : "no") << std::endl;
                //std::cout << "isNue ? " << (isNue ? "yes" : "no") <<std::endl;
                //std::cout << "isNumu ? " << (isNumu ? "yes" : "no") <<std::endl;
                //std::cout << "isCCNue ? " << (isCCNue ? "yes" : "no") <<std::endl;
                //std::cout << "isCCNumu ? " << (isCCNumu ? "yes" : "no") <<std::endl;

                for (int cheatPDG : fShowerPDGToCheat)
                {
                    const int absCheatPdg(std::abs(cheatPDG));

                    if (absPdg != absCheatPdg)
                        continue;

                    if ((absPdg == 11 && isCCNue) || (absPdg == 13 && isCCNumu))
                    {
                        cheatCharacterisation = true;
                        break;
                    }
                }
            }
            //if (cheatCharacterisation)
            //std::cout << "HighestPandrizzleScore - Am cheating shower with abs(PDG): " << absPdg << std::endl;
        }
    }

    //std::cout << "showerModuleLabel: " << (cheatCharacterisation ? fCheatShowerModuleLabel : fShowerModuleLabel) << std::endl;
    ////////////////////////////////

    //Now get the associated shower 
    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, cheatCharacterisation ? fCheatShowerModuleLabel : fShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(child_pfp.key());

    if (pfp_shower_vector.size() > 1)
    { 
      //Found a PFP with more than one shower matched.  Complain and exit
        std::cout<< "Number of associated showers to a PFP is greater than 1: " << pfp_shower_vector.size() << std::endl;
        continue;
    } 
    else if (pfp_shower_vector.size() == 0)
    { 
      //Don't bother complaining of no shower was found.  It's either missing or its a track.  No biggie here
      continue;
    }

    const art::Ptr<recob::Shower> shower = pfp_shower_vector[0];

    art::Handle< std::vector<recob::Shower> > showerListHandle;
    if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle)))
    {
        std::cout<<"Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel << std::endl;
        continue;
    }
    art::FindManyP<recob::Hit> fmhs(showerListHandle, evt, fShowerModuleLabel);
    const std::vector<art::Ptr<recob::Hit> > current_hits = fmhs.at(shower.key());

    if (current_hits.size() < fPFParticleHitCut)
        continue;

    // get the metadata here...
    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn(pfparticleListHandle, evt, "pandoraSel");
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(child_pfp.key());

    if (fDemandBDTScore)
    {
        if (pfpMetadata.size() != 1)
            continue;

        larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetadata[0]->GetPropertiesMap());

        if (propertiesMap.find("ElectronConnectionPathwayScore") == propertiesMap.end())
            continue;
    }

    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(shower, TVector3(fRecoNuVtxX, fRecoNuVtxY, fRecoNuVtxZ), evt));
    float pandrizzleBDTMethod = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);
    double pandrizzleScore = pandrizzleRecord.GetMVAScore();

    if ((std::fabs(pandrizzleBDTMethod - 2.0) < std::numeric_limits<float>::epsilon()) && (pandrizzleScore > highestJamPandrizzleScore))
    {
        highestJamPandrizzleScore = pandrizzleScore;
        selShower = shower;
        found = true;
    }

    if ((std::fabs(pandrizzleBDTMethod - 1.0) < std::numeric_limits<float>::epsilon()) && (pandrizzleScore > highestPandrizzleScore))
    {
        highestPandrizzleScore = pandrizzleScore;
        backupSelShower = shower;
    }
  }

  //return ((highestJamPandrizzleScore > highestPandrizzleScore) ? selShower : backupSelShower);

  return (found ? selShower : backupSelShower);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<art::Ptr<recob::Hit> > FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector::GetPFPHits(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
  std::vector<art::Ptr<recob::Hit> > pfp_hits;

  //Get the PFP handle out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)))
  {
      mf::LogWarning("HighestPandrizzleScoreRecoVertexShowerSelector") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
      return pfp_hits;
  }

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  if (!(evt.getByLabel(fClusterModuleLabel, clusterListHandle)))
  {
      mf::LogWarning("HighestPandrizzleScoreRecoVertexShowerSelector") << "Unable to find std::vector<recob::Cluster> with module label: " << fClusterModuleLabel;
      return pfp_hits;
  }

  art::FindManyP<recob::Cluster> fmcpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
  const std::vector<art::Ptr<recob::Cluster> > sel_pfp_clusters = fmcpfp.at(pfp.key());
  art::FindManyP<recob::Hit> fmhc(clusterListHandle, evt, fClusterModuleLabel);

  //Loop over the clusters and retrieve the hits
  for (unsigned i_cluster = 0; i_cluster < sel_pfp_clusters.size(); i_cluster++)
  {
      art::Ptr<recob::Cluster> cluster = sel_pfp_clusters[i_cluster];
      const std::vector<art::Ptr<recob::Hit> > sel_pfp_hits = fmhc.at(cluster.key());

    for (unsigned i_hit = 0; i_hit < sel_pfp_hits.size(); i_hit++)
        pfp_hits.push_back(sel_pfp_hits[i_hit]);
  }

  return pfp_hits;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

