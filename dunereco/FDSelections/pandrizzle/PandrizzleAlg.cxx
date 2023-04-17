///////////////////////////////////////////////
// PandrizzleAlg.cxx
//
// D Brailsford
///////////////////////////////////////////////

//STL
#include <limits>

//ART
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"

//Custom
#include "PandrizzleAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Helpers to convert 
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"


namespace
{
    constexpr Float_t kDefValue(std::numeric_limits<Float_t>::lowest());

    using namespace FDSelection;
    void Reset(PandrizzleAlg::InputVarsToReader &inputVarsToReader)
    {
        for (PandrizzleAlg::Vars var = PandrizzleAlg::kEvalRatio; var < PandrizzleAlg::kTerminatingValue; var=static_cast<PandrizzleAlg::Vars>(static_cast<int>(var)+1)) 
        {
            auto [itr, inserted] = inputVarsToReader.try_emplace(var, std::make_unique<Float_t>(kDefValue));
            if (!inserted)
                *(itr->second.get()) = kDefValue;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandrizzleAlg::Record::Record(const InputVarsToReader &inputVarsToReader, const Float_t mvaScore, const bool isFilled) :
    fMVAScore(mvaScore),
    fIsFilled(isFilled)
{
    for (PandrizzleAlg::Vars var = PandrizzleAlg::kEvalRatio; var < PandrizzleAlg::kTerminatingValue; var=static_cast<PandrizzleAlg::Vars>(static_cast<int>(var)+1)) 
        fInputs.try_emplace(var, *(inputVarsToReader.at(var)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandrizzleAlg::PandrizzleAlg(const fhicl::ParameterSet& pset) :
    fPFParticleModuleLabel(pset.get<std::string>("ModuleLabels.PFParticleModuleLabel")),
    fShowerModuleLabel(pset.get<std::string>("ModuleLabels.ShowerModuleLabel")),
    fClusterModuleLabel(pset.get<std::string>("ModuleLabels.ClusterModuleLabel")),
    fPIDModuleLabel(pset.get<std::string>("ModuleLabels.PIDModuleLabel")),
    fPandrizzleWeightFileName(pset.get< std::string > ("PandrizzleWeightFileName")),
    fJamPandrizzleWeightFileName(pset.get< std::string > ("JamPandrizzleWeightFileName", "")),
    fCheatCharacterisation(pset.get<bool>("CheatCharacterisation", false)),
    fShiftDisplacement(pset.get<bool>("ShiftDisplacement", false)),
    fShiftdEdX(pset.get<bool>("ShiftdEdX", false)),
    fReader("", 0),
    fJamReader("", 0),
    fMakeTree(pset.get<bool>("MakeTree", false)),
    fUseConcentration(pset.get<bool>("UseConcentration", true)),
    fUseDisplacement(pset.get<bool>("UseDisplacement", true)),
    fUseDCA(pset.get<bool>("UseDCA", true)),
    fUseBDTVariables(pset.get<bool>("UseBDTVariables", false)),
    fUseModularShowerVariables(pset.get<bool>("UseModularShowerVariables", false)),
    fEnhancedPandrizzleHitCut(pset.get<int>("EnhancedPandrizzleHitCut", 100)),
    fModularShowerPandrizzleHitCut(pset.get<int>("ModularShowerPandrizzleHitCut", 25))
{
    if (fCheatCharacterisation)
    {
        //std::cout << "PandrizzleAlg - fCheatCharacterisation is true" << std::endl; 

        fCheatShowerModuleLabel = pset.get< std::string >("ModuleLabels.CheatShowerModuleLabel");
        fCheatPIDModuleLabel = pset.get<std::string>("ModuleLabels.CheatPIDModuleLabel");
        fNuGenModuleLabel = pset.get<std::string>("ModuleLabels.NuGenModuleLabel");
        fShowerPDGToCheat = pset.get<std::vector<int>> ("ShowerPDGToCheat");
    }

    Reset(fInputsToReader);

    fReader.AddVariable("EvalRatio",GetVarPtr(kEvalRatio));

    if (fUseConcentration)
    {
        fReader.AddVariable("Concentration",GetVarPtr(kConcentration));
    }

    fReader.AddVariable("CoreHaloRatio",GetVarPtr(kCoreHaloRatio));
    fReader.AddVariable("Conicalness",GetVarPtr(kConicalness));
    fReader.AddVariable("dEdxBestPlane",GetVarPtr(kdEdxBestPlane));

    if (fUseDisplacement)
    {
        fReader.AddVariable("Displacement",GetVarPtr(kDisplacement));
    }

    if (fUseDCA)
    {
        fReader.AddVariable("DCA",GetVarPtr(kDCA));
    }

    fReader.AddVariable("Wideness",GetVarPtr(kWideness));
    fReader.AddVariable("EnergyDensity",GetVarPtr(kEnergyDensity));


    fJamReader.AddVariable("EvalRatio",GetVarPtr(kEvalRatio));
    fJamReader.AddVariable("Concentration",GetVarPtr(kConcentration));
    fJamReader.AddVariable("CoreHaloRatio",GetVarPtr(kCoreHaloRatio));
    fJamReader.AddVariable("Conicalness",GetVarPtr(kConicalness));
    fJamReader.AddVariable("dEdxBestPlane",GetVarPtr(kdEdxBestPlane));
    fJamReader.AddVariable("Displacement",GetVarPtr(kDisplacement));

    if (fUseDCA)
    {
        fJamReader.AddVariable("DCA",GetVarPtr(kDCA));
    }


    fJamReader.AddVariable("Wideness",GetVarPtr(kWideness));
    fJamReader.AddVariable("EnergyDensity",GetVarPtr(kEnergyDensity));

    if (fUseBDTVariables)
    {
        fJamReader.AddVariable("PathwayLengthMin", GetVarPtr(kPathwayLengthMin));
        fJamReader.AddVariable("MaxShowerStartPathwayScatteringAngle2D", GetVarPtr(kMaxShowerStartPathwayScatteringAngle2D));
        fJamReader.AddVariable("MaxNPostShowerStartHits", GetVarPtr(kMaxNPostShowerStartHits));
        fJamReader.AddVariable("MaxPostShowerStartScatterAngle", GetVarPtr(kMaxPostShowerStartScatterAngle));
        fJamReader.AddVariable("MaxPostShowerStartNuVertexEnergyAsymmetry", GetVarPtr(kMaxPostShowerStartNuVertexEnergyAsymmetry));
        fJamReader.AddVariable("MaxPostShowerStartShowerStartEnergyAsymmetry", GetVarPtr(kMaxPostShowerStartShowerStartEnergyAsymmetry));
        fJamReader.AddVariable("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", GetVarPtr(kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance));
        fJamReader.AddVariable("MinPostShowerStartShowerStartMoliereRadius", GetVarPtr(kMinPostShowerStartShowerStartMoliereRadius));
        fJamReader.AddVariable("MaxPostShowerStartOpeningAngle", GetVarPtr(kMaxPostShowerStartOpeningAngle));
        fJamReader.AddVariable("MaxFoundHitRatio", GetVarPtr(kMaxFoundHitRatio));
        fJamReader.AddVariable("MaxInitialGapSize", GetVarPtr(kMaxInitialGapSize));
        fJamReader.AddVariable("MinLargestProjectedGapSize", GetVarPtr(kMinLargestProjectedGapSize));
        fJamReader.AddVariable("NViewsWithAmbiguousHits", GetVarPtr(kNViewsWithAmbiguousHits));
        fJamReader.AddVariable("AmbiguousHitMaxUnaccountedEnergy", GetVarPtr(kAmbiguousHitMaxUnaccountedEnergy));

        const std::string jamWeightFileName(fJamPandrizzleWeightFileName);
        std::string jamWeightFilePath;
        cet::search_path jSP("FW_SEARCH_PATH");
        jSP.find_file(jamWeightFileName, jamWeightFilePath);

        fJamReader.BookMVA("BDTG", jamWeightFilePath);
    }

    if (fUseModularShowerVariables)
    {
        fReader.AddVariable("ModularShowerPathwayLengthMin",GetVarPtr(kModularShowerPathwayLengthMin));
        fReader.AddVariable("ModularShowerMaxNShowerHits",GetVarPtr(kModularShowerMaxNShowerHits));
        fReader.AddVariable("ModularShowerMaxNuVertexChargeWeightedMeanRadialDistance",GetVarPtr(kModularShowerMaxNuVertexChargeWeightedMeanRadialDistance));
    }

    const std::string weightFileName(fPandrizzleWeightFileName);
    std::string weightFilePath;

    cet::search_path sP("FW_SEARCH_PATH");
    sP.find_file(weightFileName, weightFilePath);

    std::cout << "weightFileName: " << weightFileName << std::endl;
    std::cout << "weightFilePath: " << weightFilePath << std::endl;

    fReader.BookMVA("BDTG",weightFilePath);

    if (fMakeTree)
    {
        InitialiseTrees();
        ResetTreeVariables();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::InitialiseTrees() {
  fSignalShowerTree = tfs->make<TTree>("DrizzleSigShowerTree","Pandrizzle Signal Shower Tree");
  fBackgroundShowerTree = tfs->make<TTree>("DrizzleBgShowerTree","Pandrizzle Background Shower Tree");

  std::map<std::string, TTree*> treeMap;
  treeMap["signalShower"] = fSignalShowerTree;
  treeMap["backgroundShower"] = fBackgroundShowerTree;
  for (std::map<std::string, TTree*>::iterator mapIt = treeMap.begin(); mapIt != treeMap.end(); mapIt++){
    TTree *tree = mapIt->second;
    BookTreeInt(tree, "Event");
    BookTreeInt(tree, "Run");
    BookTreeInt(tree, "SubRun");
    BookTreeInt(tree, "TruePDG");
    BookTreeFloat(tree, "TrueEnergy");
    BookTreeInt(tree, "PFPPDG");
    BookTreeInt(tree, "PFPNHits");
    BookTreeFloat(tree, "EvalRatio");
    BookTreeFloat(tree, "Concentration");
    BookTreeFloat(tree, "CoreHaloRatio");
    BookTreeFloat(tree, "Conicalness");
    BookTreeFloat(tree, "dEdxBestPlane");
    BookTreeFloat(tree, "Displacement");
    BookTreeFloat(tree, "DCA");
    BookTreeFloat(tree, "Wideness");
    BookTreeFloat(tree, "EnergyDensity");

    // Enhancement variables
    BookTreeFloat(tree, "FoundConnectionPathway");
    BookTreeFloat(tree, "ConnectionBDTScore");
    BookTreeFloat(tree, "PathwayLengthMin");
    BookTreeFloat(tree, "MaxShowerStartPathwayScatteringAngle2D");
    BookTreeFloat(tree, "MaxNPostShowerStartHits");
    BookTreeFloat(tree, "MaxPostShowerStartScatterAngle");
    BookTreeFloat(tree, "MaxPostShowerStartNuVertexEnergyAsymmetry");
    BookTreeFloat(tree, "MaxPostShowerStartShowerStartEnergyAsymmetry");
    BookTreeFloat(tree, "MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance");
    BookTreeFloat(tree, "MinPostShowerStartShowerStartMoliereRadius");
    BookTreeFloat(tree, "MaxPostShowerStartOpeningAngle");
    BookTreeFloat(tree, "MaxFoundHitRatio");
    BookTreeFloat(tree, "MaxInitialGapSize");
    BookTreeFloat(tree, "MinLargestProjectedGapSize");
    BookTreeFloat(tree, "NViewsWithAmbiguousHits");
    BookTreeFloat(tree, "AmbiguousHitMaxUnaccountedEnergy");

    // Backup variables
    BookTreeFloat(tree, "FoundTrackStub"); // 1.0 for yes -1.0 for no (why? because root TMVA is trash)
    BookTreeFloat(tree, "ModularShowerPathwayLengthMin");
    BookTreeFloat(tree, "ModularShowerPathwayKink3D");
    BookTreeFloat(tree, "ModularShowerMaxNShowerHits");
    BookTreeFloat(tree, "ModularShowerMaxScatterAngle"); //looks weird
    BookTreeFloat(tree, "ModularShowerMaxNuVertexChargeAsymmetry");
    BookTreeFloat(tree, "ModularShowerMaxShowerStartChargeAsymmetry");
    BookTreeFloat(tree, "ModularShowerMaxNuVertexChargeWeightedMeanRadialDistance");
    BookTreeFloat(tree, "ModularShowerMinShowerStartMoliereRadius");
    BookTreeFloat(tree, "ModularShowerMaxOpeningAngle");
    BookTreeFloat(tree, "ModularShowerMaxFoundHitRatio");
    BookTreeFloat(tree, "ModularShowerProjectedInitialGapSize");
    BookTreeFloat(tree, "ModularShowerGlobalInitialGapSize");
    BookTreeFloat(tree, "ModularShowerLargestProjectedGapSize");

    // For drawing purposes
    BookTreeFloat(tree, "StartX");
    BookTreeFloat(tree, "StartY");
    BookTreeFloat(tree, "StartZ");
    BookTreeFloat(tree, "EndX");
    BookTreeFloat(tree, "EndY");
    BookTreeFloat(tree, "EndZ");
    BookTreeFloat(tree, "StartPX");
    BookTreeFloat(tree, "StartPY");
    BookTreeFloat(tree, "StartPZ");
    BookTreeFloat(tree, "EndPX");
    BookTreeFloat(tree, "EndPY");
    BookTreeFloat(tree, "EndPZ");

    fVarHolder.m_trajPositionX = std::vector<float>();
    fVarHolder.m_trajPositionY = std::vector<float>();
    fVarHolder.m_trajPositionZ = std::vector<float>();
    tree->Branch("TrajPositionX", &fVarHolder.m_trajPositionX);
    tree->Branch("TrajPositionY", &fVarHolder.m_trajPositionY);
    tree->Branch("TrajPositionZ", &fVarHolder.m_trajPositionZ);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::BookTreeInt(TTree *tree, std::string branch_name)
{
    fVarHolder.IntVars[branch_name] = -9999;
    tree->Branch(branch_name.c_str(), &(fVarHolder.IntVars[branch_name]));
    return;
}

void FDSelection::PandrizzleAlg::BookTreeFloat(TTree *tree, std::string branch_name)
{
    fVarHolder.FloatVars[branch_name] = -9999.f;
    tree->Branch(branch_name.c_str(), &(fVarHolder.FloatVars[branch_name]));
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::Run(const art::Event& evt) {

  //Get the PFPs out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)))
  {
      mf::LogWarning("PandrizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
      return;
  }
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

  //Ceate the full PFP map
  lar_pandora::PFParticleMap pfparticleMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleList, pfparticleMap);

  //Grab the primary PFPs (the neutrinos) from the 
  std::vector<art::Ptr<recob::PFParticle> > nu_pfps;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

  //Now grab the primary children of these PFP
  for (unsigned int i_nupfp = 0; i_nupfp < nu_pfps.size(); i_nupfp++)
  {
      art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[i_nupfp];
      std::vector<art::Ptr<recob::PFParticle> > child_pfps = SelectChildPFParticles(nu_pfp, pfparticleMap);

      //Assess each child pfp
      for (unsigned int i_childpfp = 0; i_childpfp < child_pfps.size(); i_childpfp++)
      {
          art::Ptr<recob::PFParticle> child_pfp = child_pfps[i_childpfp];

          //Process the child PFP
          ProcessPFParticle(child_pfp, evt);
          FillTree();
          ResetTreeVariables();
      }
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<art::Ptr<recob::PFParticle> > FDSelection::PandrizzleAlg::SelectChildPFParticles(const art::Ptr<recob::PFParticle> parent_pfp, const lar_pandora::PFParticleMap & pfp_map)
{
  std::vector<art::Ptr<recob::PFParticle> > child_pfps;

  for (int i_child = 0; i_child < parent_pfp->NumDaughters(); i_child++)
  {
      int child_id = parent_pfp->Daughter(i_child);
      art::Ptr<recob::PFParticle> child_pfp = pfp_map.at(child_id);
      child_pfps.push_back(child_pfp);
  }

  return child_pfps;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
    art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
    if (evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle))
        art::fill_ptr_vector(pfparticleList, pfparticleListHandle);

    /*
    art::Handle< std::vector<larpandoraobj::PFParticleMetadata> > pfpMetadataListHandle;
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetadataList; 
    if (evt.getByLabel("pandoraSel", pfpMetadataListHandle)){
      art::fill_ptr_vector(pfpMetadataList, pfpMetadataListHandle);
    }
    */

    //Get the matched MCParticle
    simb::MCParticle* matched_mcparticle(nullptr);
    std::vector<art::Ptr<recob::Hit> > pfp_hits = GetPFPHits(pfp, evt);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfp_hits, 1);
    if (g4id > 0)
    {
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        matched_mcparticle = pi_serv->ParticleList().at(g4id);
    }

    std::string showerModuleLabel(fShowerModuleLabel);
    std::string pidModuleLabel(fPIDModuleLabel);


    // for cheating characterisation studies
    if (matched_mcparticle && !fShowerPDGToCheat.empty())
    {
        const int absPdg(std::abs(matched_mcparticle->PdgCode()));

        for (int cheatPDG : fShowerPDGToCheat)
        {
            const int absCheatPdg(std::abs(cheatPDG));

            if (absPdg != absCheatPdg)
                continue;

            if ((absCheatPdg == 11 || absCheatPdg == 13))
                continue;

            showerModuleLabel = fCheatShowerModuleLabel;
            pidModuleLabel = fCheatPIDModuleLabel;
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

            for (int cheatPDG : fShowerPDGToCheat)
            {
                const int absCheatPdg(std::abs(cheatPDG));

                if (absPdg != absCheatPdg)
                    continue;

                if ((absPdg == 11 && isCCNue) || (absPdg == 13 && isCCNumu))
                {
                    showerModuleLabel = fCheatShowerModuleLabel;
                    pidModuleLabel = fCheatPIDModuleLabel;
                    break;
                }
            }
        }

        //if (showerModuleLabel == "pandoraShowerSelCheat")
        //std::cout << "PandrizzleAlg - Am cheating shower with abs(PDG): " << absPdg << std::endl;
    }

    //std::cout << "PandrizzleAlg - showerModuleLabel: " << showerModuleLabel << std::endl;
    //std::cout << "PandrizzleAlg - pidModuleLabel: " << pidModuleLabel << std::endl;

    art::FindManyP<recob::Shower> fmspfp(pfparticleListHandle, evt, showerModuleLabel);
    const std::vector<art::Ptr<recob::Shower> > pfp_shower_vector = fmspfp.at(pfp.key());

    if (pfp_shower_vector.size() > 1)
    { 
        //Found a PFP with more than one shower matched.  Complain and exit
        std::cout<< "Number of associated showers to a PFP is greater than 1: " << pfp_shower_vector.size() << std::endl;
        return;
    } 
    else if (pfp_shower_vector.size() == 0)
    { 
        //Don't bother complaining of no shower was found.  It's either missing or its a track.  No biggie here
        return;
    }

    const art::Ptr<recob::Shower> pShower = pfp_shower_vector[0];

    lar_pandora::PFParticleVector nu_pfps;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticleList, nu_pfps);

    // Check we have one neutrino...
    if (nu_pfps.size() != 1)
    {
        std::cout << "PandrizzleAlg: Number of Nu PFPs does not equal 1: " << nu_pfps.size() << std::endl;
        return;
    }

    art::Ptr<recob::PFParticle> nu_pfp = nu_pfps[0];

    art::FindManyP<recob::Vertex> fmvpfp(pfparticleListHandle, evt, fPFParticleModuleLabel);
    const std::vector<art::Ptr<recob::Vertex> > nu_vertices = fmvpfp.at(nu_pfp.key());

    if (nu_vertices.size() == 0)
    {
        return;
    }
    else if (nu_vertices.size() > 1)
    {
        std::cout << "PandrizzleAlg: Number of matched neutrino vertices bigger than 1: " << nu_vertices.size() << std::endl;
    }

    //always take the first vertex, even if there's more than one
    TVector3 nuVertex = TVector3(nu_vertices[0]->position().X(), nu_vertices[0]->position().Y(), nu_vertices[0]->position().Z());


    // SIGNED STUFF
    /*
    double displacementX = pShower->ShowerStart().X() - nu_vertices[0]->position().X();
    double displacementY = pShower->ShowerStart().Y() - nu_vertices[0]->position().Y();
    double displacementZ = pShower->ShowerStart().Z() - nu_vertices[0]->position().Z();

    double displacement = std::sqrt((displacementX * displacementX) + (displacementY * displacementY) + (displacementZ * displacementZ));



    double dot = (pShower->Direction().X() * displacementX) + (pShower->Direction().Y() * displacementY) + (pShower->Direction().Z() * displacementZ);
    double showerDirectionMag = std::sqrt((pShower->Direction().X() * pShower->Direction().X()) + (pShower->Direction().Y() * pShower->Direction().Y()) +
                                          (pShower->Direction().Z() * pShower->Direction().Z()));

    if (std::fabs(showerDirectionMag) < std::numeric_limits<float>::epsilon())
        std::cout << "displacement: (" << displacementX << ", " << displacementY << ", " << displacementZ << ")" << std::endl;

    double cosOpeningAngle = dot / (showerDirectionMag * displacement);
    bool nuVertexDownstream(false);

    if (cosOpeningAngle < 0.f)
        nuVertexDownstream = true;

    if (nuVertexDownstream)
        displacement = displacement * (-1.f);
    */

    art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower> >{pShower}, evt, pidModuleLabel);
    art::Ptr<anab::MVAPIDResult> mvaPIDResult(findPIDResult.at(0));

    //Get event,subrun,run
    fVarHolder.IntVars["Run"] = evt.run();
    fVarHolder.IntVars["SubRun"] = evt.subRun();
    fVarHolder.IntVars["Event"] = evt.event();

    //Fill the PFP hits
    fVarHolder.IntVars["PFPPDG"] = pfp->PdgCode();
    fVarHolder.IntVars["PFPNHits"] = pfp_hits.size();

    if (matched_mcparticle)
    {
        fVarHolder.IntVars["TruePDG"] = matched_mcparticle->PdgCode();
        fVarHolder.FloatVars["TrueEnergy"] = matched_mcparticle->Momentum().T();
    }

    //MVAPID vars
    if (mvaPIDResult.isAvailable())
    {
        fVarHolder.FloatVars["EvalRatio"] = (isnan(mvaPIDResult->evalRatio) ? -0.5f : static_cast<Float_t>(mvaPIDResult->evalRatio));
        fVarHolder.FloatVars["Concentration"] = (isnan(mvaPIDResult->concentration) ? -2.f : std::min(static_cast<Float_t>(mvaPIDResult->concentration), 50.f));
        fVarHolder.FloatVars["CoreHaloRatio"] = static_cast<Float_t>(mvaPIDResult->coreHaloRatio);
        fVarHolder.FloatVars["Conicalness"] = std::min(static_cast<Float_t>(mvaPIDResult->conicalness), 100.f);
    }

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn(pfparticleListHandle, evt, "pandoraSel");
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(pfp.key());

    //dEdx
    if (pShower->dEdx().size() > 0)
        fVarHolder.FloatVars["dEdxBestPlane"] = std::max(std::min(static_cast<Float_t>(pShower->dEdx().at(pShower->best_plane())), 20.f), -2.f);

    if (fShiftdEdX && (pfpMetadata.size() == 1))
    {
        larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetadata[0]->GetPropertiesMap());

        if (propertiesMap.find("dEdX") != propertiesMap.end())
        {
            fVarHolder.FloatVars["dEdxBestPlane"] = propertiesMap.at("dEdX");
        }
    }

    //Displacement
    TVector3 showerVertex = pShower->ShowerStart();
    if (fShiftDisplacement && (pfpMetadata.size() == 1))
    {
      larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetadata[0]->GetPropertiesMap());

      if ((propertiesMap.find("ShowerVertexX") != propertiesMap.end()) && (propertiesMap.find("ShowerVertexY") != propertiesMap.end()) && 
          (propertiesMap.find("ShowerVertexZ") != propertiesMap.end()))
      {
          showerVertex = TVector3(propertiesMap.at("ShowerVertexX"), propertiesMap.at("ShowerVertexY"), propertiesMap.at("ShowerVertexZ"));
      }
    }

    fVarHolder.FloatVars["Displacement"] = std::min(static_cast<Float_t>((showerVertex - nuVertex).Mag()), 100.f);

    //Distance of closest approach
    double alpha((pShower->ShowerStart() - nuVertex).Dot(pShower->Direction()));
    TVector3 r(pShower->ShowerStart() + alpha*pShower->Direction());

    fVarHolder.FloatVars["DCA"] = std::min(static_cast<Float_t>((r-nuVertex).Mag()), 50.f);

    //Wideness
    Float_t wideness(static_cast<Float_t>(pShower->OpenAngle()/pShower->Length()));

    fVarHolder.FloatVars["Wideness"] = isnan(wideness) ? -0.01f : std::min(wideness, 0.1f);

    //Energy density
    if (pShower->Energy().size() > 0)
    {
        Float_t volume(static_cast<Float_t>((M_PI * pShower->Length() * pShower->Length() * pShower->Length() * std::tan(pShower->OpenAngle()))/3.f));
        Float_t energyDensity(std::min(std::max(static_cast<Float_t>(pShower->Energy().at(2))/volume, -0.1f), 5.f));

        fVarHolder.FloatVars["EnergyDensity"] = isnan(energyDensity) ? -0.1f : energyDensity;
    }

    if (pfpMetadata.size() == 1)
    {
        larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetadata[0]->GetPropertiesMap());

        fVarHolder.FloatVars["FoundConnectionPathway"] = propertiesMap.find("ElectronConnectionPathwayScore") != propertiesMap.end() ? 1.0 : 0.0;
        fVarHolder.FloatVars["ConnectionBDTScore"] = propertiesMap.find("ElectronConnectionPathwayScore") != propertiesMap.end() ? propertiesMap.at("ElectronConnectionPathwayScore") : -100.0;
        fVarHolder.FloatVars["PathwayLengthMin"] = propertiesMap.find("PathwayLengthMin") != propertiesMap.end() ? propertiesMap.at("PathwayLengthMin") : -100.0;
        fVarHolder.FloatVars["MaxShowerStartPathwayScatteringAngle2D"] = propertiesMap.find("MaxShowerStartPathwayScatteringAngle2D") != propertiesMap.end() ? 
          propertiesMap.at("MaxShowerStartPathwayScatteringAngle2D") : -100.0;
        fVarHolder.FloatVars["MaxNPostShowerStartHits"] = propertiesMap.find("MaxNPostShowerStartHits") != propertiesMap.end() ? propertiesMap.at("MaxNPostShowerStartHits") : -100.0;
        fVarHolder.FloatVars["MaxPostShowerStartScatterAngle"] = propertiesMap.find("MaxPostShowerStartScatterAngle") != propertiesMap.end() ? propertiesMap.at("MaxPostShowerStartScatterAngle") : -100.0;
        fVarHolder.FloatVars["MaxPostShowerStartNuVertexEnergyAsymmetry"] = propertiesMap.find("MaxPostShowerStartNuVertexEnergyAsymmetry") != propertiesMap.end() ? 
          propertiesMap.at("MaxPostShowerStartNuVertexEnergyAsymmetry") : -100.0;
        fVarHolder.FloatVars["MaxPostShowerStartShowerStartEnergyAsymmetry"] = propertiesMap.find("MaxPostShowerStartShowerStartEnergyAsymmetry") != propertiesMap.end() ? 
          propertiesMap.at("MaxPostShowerStartShowerStartEnergyAsymmetry") : -100.0;
        fVarHolder.FloatVars["MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance"] = propertiesMap.find("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance") != propertiesMap.end() ? 
          propertiesMap.at("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance") : -100.0;
        fVarHolder.FloatVars["MinPostShowerStartShowerStartMoliereRadius"] = propertiesMap.find("MinPostShowerStartShowerStartMoliereRadius") != propertiesMap.end() ? 
          propertiesMap.at("MinPostShowerStartShowerStartMoliereRadius") : -100.0;
        fVarHolder.FloatVars["MaxPostShowerStartOpeningAngle"] = propertiesMap.find("MaxPostShowerStartOpeningAngle") != propertiesMap.end() ? 
            propertiesMap.at("MaxPostShowerStartOpeningAngle") : -100.0;
        fVarHolder.FloatVars["MaxFoundHitRatio"] = propertiesMap.find("MaxFoundHitRatio") != propertiesMap.end() ? propertiesMap.at("MaxFoundHitRatio") : -100.0;
        fVarHolder.FloatVars["MaxInitialGapSize"] = propertiesMap.find("MaxInitialGapSize") != propertiesMap.end() ? propertiesMap.at("MaxInitialGapSize") : -100.0;
        fVarHolder.FloatVars["MinLargestProjectedGapSize"] = propertiesMap.find("MinLargestProjectedGapSize") != propertiesMap.end() ?
            propertiesMap.at("MinLargestProjectedGapSize") : -100.0;
        fVarHolder.FloatVars["NViewsWithAmbiguousHits"] = propertiesMap.find("NViewsWithAmbiguousHits") != propertiesMap.end() ? propertiesMap.at("NViewsWithAmbiguousHits") : -100.0;
        fVarHolder.FloatVars["AmbiguousHitMaxUnaccountedEnergy"] = propertiesMap.find("AmbiguousHitMaxUnaccountedEnergy") != propertiesMap.end() ? 
          propertiesMap.at("AmbiguousHitMaxUnaccountedEnergy") : -100.0;
    }
    else
    {
        fVarHolder.FloatVars["FoundConnectionPathway"] = 0.0;
        fVarHolder.FloatVars["ConnectionBDTScore"] = -100.0; 
        fVarHolder.FloatVars["PathwayLengthMin"] = -100.0; 
        fVarHolder.FloatVars["MaxShowerStartPathwayScatteringAngle2D"] = -100.0; 
        fVarHolder.FloatVars["MaxNPostShowerStartHits"] = -100.0; 
        fVarHolder.FloatVars["MaxPostShowerStartScatterAngle"] = -100.0; 
        fVarHolder.FloatVars["MaxPostShowerStartNuVertexEnergyAsymmetry"] = -100.0; 
        fVarHolder.FloatVars["MaxPostShowerStartShowerStartEnergyAsymmetry"] = -100.0; 
        fVarHolder.FloatVars["MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance"] = -100.0; 
        fVarHolder.FloatVars["MinPostShowerStartShowerStartMoliereRadius"] = -100.0; 
        fVarHolder.FloatVars["MaxPostShowerStartOpeningAngle"] = -100.0;
        fVarHolder.FloatVars["MaxFoundHitRatio"] = -100.0;
        fVarHolder.FloatVars["MaxInitialGapSize"] = -100.0;
        fVarHolder.FloatVars["MinLargestProjectedGapSize"] = -100.0;
        fVarHolder.FloatVars["NViewsWithAmbiguousHits"] = -100.0; 
        fVarHolder.FloatVars["AmbiguousHitMaxUnaccountedEnergy"] = -100.0; 
    }

    art::Handle< std::vector<recob::Shower> > showerListHandle;
    evt.getByLabel(fShowerModuleLabel, showerListHandle);
    art::FindManyP<recob::Track> initialTrackAssn(showerListHandle, evt, fShowerModuleLabel);
    std::vector<art::Ptr<recob::Track>> initialTrackStub = initialTrackAssn.at(pShower.key());

    if (initialTrackStub.size() == 1)
    {
        art::Ptr<recob::Track> trackStub = initialTrackStub.at(0);

        /*
        // Wanna see what the bad boi track stubs look like?
        fVarHolder.FloatVars["StartX"] = initialTrackStub.at(0)->Start().X();
        fVarHolder.FloatVars["StartY"] = initialTrackStub.at(0)->Start().Y();
        fVarHolder.FloatVars["StartZ"] = initialTrackStub.at(0)->Start().Z();
        fVarHolder.FloatVars["EndX"] = initialTrackStub.at(0)->End().X();
        fVarHolder.FloatVars["EndY"] = initialTrackStub.at(0)->End().Y();
        fVarHolder.FloatVars["EndZ"] = initialTrackStub.at(0)->End().Z();
        fVarHolder.FloatVars["StartPX"] = initialTrackStub.at(0)->StartDirection().X();
        fVarHolder.FloatVars["StartPY"] = initialTrackStub.at(0)->StartDirection().Y();
        fVarHolder.FloatVars["StartPZ"] = initialTrackStub.at(0)->StartDirection().Z();
        fVarHolder.FloatVars["EndPX"] = initialTrackStub.at(0)->EndDirection().X();
        fVarHolder.FloatVars["EndPY"] = initialTrackStub.at(0)->EndDirection().Y();
        fVarHolder.FloatVars["EndPZ"] = initialTrackStub.at(0)->EndDirection().Z();

        int nTrajectoryPoints(initialTrackStub.at(0)->NumberTrajectoryPoints());

        for (int i = 0; i < nTrajectoryPoints; ++i)
        {
            fVarHolder.m_trajPositionX.push_back(initialTrackStub.at(0)->LocationAtPoint(i).X());
            fVarHolder.m_trajPositionY.push_back(initialTrackStub.at(0)->LocationAtPoint(i).Y());
            fVarHolder.m_trajPositionZ.push_back(initialTrackStub.at(0)->LocationAtPoint(i).Z());
        }
        */

        fVarHolder.FloatVars["FoundTrackStub"] = 1.0;

        float modularShowerPathwayLengthMin(-10.f), modularShowerPathwayKink3D(-10.f);

        GetPathwayVariables(trackStub, modularShowerPathwayLengthMin, modularShowerPathwayKink3D);

        fVarHolder.FloatVars["ModularShowerPathwayLengthMin"] = std::min(modularShowerPathwayLengthMin, 30.f);
        fVarHolder.FloatVars["ModularShowerPathwayKink3D"] = std::min(modularShowerPathwayKink3D, 20.f);

        SetInitialRegionVariables(nuVertex, trackStub);

        float modularShowerMaxNShowerHits(-10.f), modularShowerMaxFoundHitRatio(-10.f), modularShowerMaxOpeningAngle(-10.f), modularShowerMaxScatterAngle(-10.f), 
            modularShowerMaxNuVertexChargeAsymmetry(-10.f), modularShowerMaxShowerStartChargeAsymmetry(-10.f), modularShowerMaxNuVertexChargeWeightedMeanRadialDistance(-10.f), 
            modularShowerMinShowerStartMoliereRadius(-10.f);

        GetShowerRegionVariables(nuVertex, pfp, trackStub, modularShowerMaxNShowerHits, modularShowerMaxFoundHitRatio, modularShowerMaxOpeningAngle, modularShowerMaxScatterAngle, 
            modularShowerMaxNuVertexChargeAsymmetry, modularShowerMaxShowerStartChargeAsymmetry, modularShowerMaxNuVertexChargeWeightedMeanRadialDistance, 
            modularShowerMinShowerStartMoliereRadius, evt);

        fVarHolder.FloatVars["ModularShowerMaxNShowerHits"] = std::min(modularShowerMaxNShowerHits, 2000.f);
        fVarHolder.FloatVars["ModularShowerMaxFoundHitRatio"] = std::min(modularShowerMaxFoundHitRatio, 1.5f);
        fVarHolder.FloatVars["ModularShowerMaxOpeningAngle"] = std::min(modularShowerMaxOpeningAngle, 20.f);
        fVarHolder.FloatVars["ModularShowerMaxScatterAngle"] = std::min(modularShowerMaxScatterAngle, 40.f);
        fVarHolder.FloatVars["ModularShowerMaxNuVertexChargeAsymmetry"] = modularShowerMaxNuVertexChargeAsymmetry;
        fVarHolder.FloatVars["ModularShowerMaxShowerStartChargeAsymmetry"] = modularShowerMaxShowerStartChargeAsymmetry;
        fVarHolder.FloatVars["ModularShowerMaxNuVertexChargeWeightedMeanRadialDistance"] = std::min(modularShowerMaxNuVertexChargeWeightedMeanRadialDistance, 20.f);
        fVarHolder.FloatVars["ModularShowerMinShowerStartMoliereRadius"] = std::min(modularShowerMinShowerStartMoliereRadius, 10.f);
    }
    else
    {
        fVarHolder.FloatVars["FoundTrackStub"] = -1.0;
        fVarHolder.FloatVars["ModularShowerPathwayLengthMin"] = -100.0;
        fVarHolder.FloatVars["ModularShowerPathwayKink3D"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxNShowerHits"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxFoundHitRatio"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxOpeningAngle"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxScatterAngle"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxNuVertexChargeAsymmetry"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxShowerStartChargeAsymmetry"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMaxNuVertexChargeWeightedMeanRadialDistance"] = -100.0;
        fVarHolder.FloatVars["ModularShowerMinShowerStartMoliereRadius"] = -100.0;
        fVarHolder.FloatVars["ModularShowerProjectedInitialGapSize"] = -100.0;
        fVarHolder.FloatVars["ModularShowerGlobalInitialGapSize"] = -100.0;
        fVarHolder.FloatVars["ModularShowerLargestProjectedGapSize"] = -100.0;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::GetPathwayVariables(art::Ptr<recob::Track> &trackStub, float &modularShowerPathwayLengthMin, float &modularShowerPathwayKink3D)
{
    modularShowerPathwayLengthMin = std::min(std::sqrt((trackStub->Start() - trackStub->End()).Mag2()), 30.0);
    modularShowerPathwayKink3D = std::min(GetLargest3DPathwayKink(trackStub), 20.f);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

float FDSelection::PandrizzleAlg::GetLargest3DPathwayKink(art::Ptr<recob::Track> &trackStub)
{
    const int nLayersHalfWindow(5);
    const int nLayersSpanned(trackStub->NumberTrajectoryPoints());

    if (nLayersSpanned <= 3 * nLayersHalfWindow)
    {
        std::cout << "NOT ENOUGH LAYERS SPANNED" << std::endl;
        return -10.f;
    }

    const unsigned int maxCentralLayer(nLayersSpanned - (2 * nLayersHalfWindow) - 1);
    const unsigned int minCentralLayer(nLayersHalfWindow);

    float highestOpeningAngle(-10.f);

    for (unsigned int index = minCentralLayer; index <= maxCentralLayer; ++index)
    {
        bool found(false);
        float thisOpeningAngle(std::numeric_limits<float>::max());

        for (int i = 0; i < nLayersHalfWindow; ++i)
        {
            unsigned int testIndex = index + i;
            const unsigned int firstIndex(testIndex - nLayersHalfWindow);
            const unsigned int centralIndex(testIndex);
            const unsigned int secondIndex(testIndex + nLayersHalfWindow);

            const geo::Point_t firstPosition(trackStub->LocationAtPoint(firstIndex));
            const geo::Point_t centralPosition(trackStub->LocationAtPoint(centralIndex));
            const geo::Point_t secondPosition(trackStub->LocationAtPoint(secondIndex));

            geo::Vector_t firstDirection(centralPosition - firstPosition);
            firstDirection /= std::sqrt(firstDirection.Mag2());

            geo::Vector_t secondDirection(secondPosition - centralPosition);
            secondDirection /= std::sqrt(secondDirection.Mag2());

            const float openingAngle(acos(firstDirection.Dot(secondDirection)) * 180.0 / 3.14);

            if (openingAngle < thisOpeningAngle)
            {
                found = true;
                thisOpeningAngle = openingAngle;
            }
        }

        if (found && (thisOpeningAngle > highestOpeningAngle))
            highestOpeningAngle = thisOpeningAngle;
    }

    return highestOpeningAngle;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::SetInitialRegionVariables(TVector3 &nuVertex, art::Ptr<recob::Track> &trackStub)
{
    float maxProjectedGapSize(-10.f);
    float projectedInitialGapSize(-10.f);
    float globalInitialGapSize(std::numeric_limits<float>::max());

    std::vector<float> longitudinalProjections;
    const geo::Vector_t &startDirection(trackStub->StartDirection());

    geo::Vector_t nuVertexPosition(nuVertex.X(), nuVertex.Y(), nuVertex.Z());

    bool found = false;
    for (unsigned int i = 0; i < trackStub->NumberTrajectoryPoints(); ++i)
    {
        const geo::Vector_t hitPosition(trackStub->LocationAtPoint(i));
        const geo::Vector_t displacement(hitPosition - nuVertexPosition);
        const float distanceToNuVertex(std::sqrt(displacement.Mag2()));
        const float longitudinalProjection(startDirection.Dot(displacement));

        longitudinalProjections.push_back(longitudinalProjection);

        if (distanceToNuVertex < globalInitialGapSize)
        {
            found = true;
            globalInitialGapSize = distanceToNuVertex;
        }
    }

    if (!found)
        globalInitialGapSize = -10.f;

    std::sort(longitudinalProjections.begin(), longitudinalProjections.end());

    projectedInitialGapSize = std::fabs(longitudinalProjections[0]);

    const long unsigned int nSampleHits(10);
    const unsigned int nIterations(std::min(longitudinalProjections.size(), nSampleHits) - 1); // basically how many hits you want to consider

    for (unsigned int i = 0; i < nIterations; ++i)
        maxProjectedGapSize = std::max(std::fabs(longitudinalProjections[i] - longitudinalProjections[i + 1]), maxProjectedGapSize);

    fVarHolder.FloatVars["ModularShowerProjectedInitialGapSize"] = std::min(projectedInitialGapSize, 4.f);
    fVarHolder.FloatVars["ModularShowerGlobalInitialGapSize"] = std::min(globalInitialGapSize, 4.f);
    fVarHolder.FloatVars["ModularShowerLargestProjectedGapSize"] = std::min(maxProjectedGapSize, 2.f);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::GetShowerRegionVariables(const TVector3 &nuVertex, const art::Ptr<recob::PFParticle> &pfp, art::Ptr<recob::Track> &trackStub, float &modularShowerMaxNShowerHits,
    float &modularShowerMaxFoundHitRatio, float &modularShowerMaxOpeningAngle, float &modularShowerMaxScatterAngle, float &modularShowerMaxNuVertexChargeAsymmetry,
    float &modularShowerMaxShowerStartChargeAsymmetry, float &modularShowerMaxNuVertexChargeWeightedMeanRadialDistance, float &modularShowerMinShowerStartMoliereRadius, const art::Event& evt)
{
    const geo::Point_t nuVertexPosition(nuVertex.X(), nuVertex.Y(), nuVertex.Z());
    const bool isDownstream((trackStub->Start() - nuVertexPosition).Mag2() < (trackStub->End() - nuVertexPosition).Mag2());

    const pandora::CartesianVector nuVertexU(nuVertexPosition.X(), 0.0, YZtoU(nuVertexPosition.Y(), nuVertexPosition.Z()));
    const pandora::CartesianVector nuVertexV(nuVertexPosition.X(), 0.0, YZtoV(nuVertexPosition.Y(), nuVertexPosition.Z()));
    const pandora::CartesianVector nuVertexW(nuVertexPosition.X(), 0.0, YZtoW(nuVertexPosition.Y(), nuVertexPosition.Z()));

    const pandora::CartesianVector pandoraNuVertex(nuVertex.X(), nuVertex.Y(), nuVertex.Z());
    const pandora::CartesianVector connectionPathwaySeed(pandoraNuVertex + 
    (pandora::CartesianVector(trackStub->StartDirection().X(), trackStub->StartDirection().Y(), trackStub->StartDirection().Z()) * 1.0));

    //const pandora::CartesianVector connectionPathwaySeed(trackStub->Start().X(), trackStub->Start().Y(), trackStub->Start().Z());

    const pandora::CartesianVector connectionPathwaySeedU(connectionPathwaySeed.GetX(), 0.0, YZtoU(connectionPathwaySeed.GetY(), connectionPathwaySeed.GetZ()));
    const pandora::CartesianVector connectionPathwaySeedV(connectionPathwaySeed.GetX(), 0.0, YZtoV(connectionPathwaySeed.GetY(), connectionPathwaySeed.GetZ()));
    const pandora::CartesianVector connectionPathwaySeedW(connectionPathwaySeed.GetX(), 0.0, YZtoW(connectionPathwaySeed.GetY(), connectionPathwaySeed.GetZ()));

    pandora::CartesianVector connectionPathwayDirectionU((connectionPathwaySeedU - nuVertexU).GetUnitVector());
    pandora::CartesianVector connectionPathwayDirectionV((connectionPathwaySeedV - nuVertexV).GetUnitVector());
    pandora::CartesianVector connectionPathwayDirectionW((connectionPathwaySeedW - nuVertexW).GetUnitVector());

    const geo::Vector_t showerStart(trackStub->End());
    const geo::Vector_t showerStartU(showerStart.X(), 0.0, YZtoU(showerStart.Y(), showerStart.Z()));
    const geo::Vector_t showerStartV(showerStart.X(), 0.0, YZtoV(showerStart.Y(), showerStart.Z()));
    const geo::Vector_t showerStartW(showerStart.X(), 0.0, YZtoW(showerStart.Y(), showerStart.Z()));

    const geo::Vector_t showerDirectionSeed(showerStart + (trackStub->EndDirection() * 1.0));
    const geo::Vector_t showerDirectionSeedU(showerDirectionSeed.X(), 0.0, YZtoU(showerDirectionSeed.Y(), showerDirectionSeed.Z()));
    const geo::Vector_t showerDirectionSeedV(showerDirectionSeed.X(), 0.0, YZtoV(showerDirectionSeed.Y(), showerDirectionSeed.Z()));
    const geo::Vector_t showerDirectionSeedW(showerDirectionSeed.X(), 0.0, YZtoW(showerDirectionSeed.Y(), showerDirectionSeed.Z()));

    geo::Vector_t showerDirectionU(showerDirectionSeedU - showerStartU);
    geo::Vector_t showerDirectionV(showerDirectionSeedV - showerStartV);
    geo::Vector_t showerDirectionW(showerDirectionSeedW - showerStartW);

    showerDirectionU = showerDirectionU / std::sqrt(showerDirectionU.Mag2());
    showerDirectionV = showerDirectionV / std::sqrt(showerDirectionV.Mag2());
    showerDirectionW = showerDirectionW / std::sqrt(showerDirectionW.Mag2());

    // Find 2D showers

    art::ServiceHandle<geo::Geometry const> theGeometry;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    std::vector<art::Ptr<recob::Hit>> allShowerHits(dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, "pandoraSel"));
    std::vector<art::Ptr<recob::Hit>> showerHitsU, showerHitsV, showerHitsW;
    pandora::CartesianPointVector cartesianPointVectorU, cartesianPointVectorV, cartesianPointVectorW;

    int nHitsU(0), nHitsV(0), nHitsW(0);

    for (auto pHit : allShowerHits)
    {
        const geo::Vector_t hitPosition(GetPandoraHitPosition(pHit, evt));
        const geo::View_t pandoraView(GetPandoraHitView(pHit));

        if (pandoraView == geo::kW)
        {
            ++nHitsW;

            const geo::Vector_t displacement(hitPosition - showerStartW);
            const float l(showerDirectionW.Dot(displacement));
            const float t(std::sqrt(showerDirectionW.Cross(displacement).Mag2()));

            if ((l > 0.f) && (t < 14.0))
            {
                showerHitsW.push_back(pHit);
                cartesianPointVectorW.push_back(pandora::CartesianVector(hitPosition.X(), hitPosition.Y(), hitPosition.Z()));
            }
        }
        else if (pandoraView == geo::kU)
        {
            ++nHitsU;

            const geo::Vector_t displacement(hitPosition - showerStartU);
            const float l(showerDirectionU.Dot(displacement));
            const float t(std::sqrt(showerDirectionU.Cross(displacement).Mag2()));

            if ((l > 0.f) && (t < 14.0))
            {
                showerHitsU.push_back(pHit);
                cartesianPointVectorU.push_back(pandora::CartesianVector(hitPosition.X(), hitPosition.Y(), hitPosition.Z()));
            }
        }
        else if (pandoraView == geo::kV) 
        {
            ++nHitsV;

            const geo::Vector_t displacement(hitPosition - showerStartV);
            const float l(showerDirectionV.Dot(displacement));
            const float t(std::sqrt(showerDirectionV.Cross(displacement).Mag2()));

            if ((l > 0.f) && (t < 14.0))
            {
                showerHitsV.push_back(pHit);
                cartesianPointVectorV.push_back(pandora::CartesianVector(hitPosition.X(), hitPosition.Y(), hitPosition.Z()));
            }
        }
        else 
        {
            throw cet::exception("LArPandora")
                << "CreatePandoraHits2D - this wire view not recognised (View=" << pandoraView << ") ";
        }
    }

    std::cout << "showerHitsU.size(): " << showerHitsU.size() << std::endl;
    std::cout << "showerHitsV.size(): " << showerHitsV.size() << std::endl;
    std::cout << "showerHitsW.size(): " << showerHitsW.size() << std::endl;

    modularShowerMaxNShowerHits = std::min(std::max(std::max(showerHitsU.size(), showerHitsV.size()), showerHitsW.size()), static_cast<long unsigned int>(2000));

    float foundHitRatioU = showerHitsU.size() == 0 ? 0 : static_cast<float>(showerHitsU.size()) / static_cast<float>(nHitsU);
    float foundHitRatioV = showerHitsV.size() == 0 ? 0 : static_cast<float>(showerHitsV.size()) / static_cast<float>(nHitsV);
    float foundHitRatioW = showerHitsW.size() == 0 ? 0 : static_cast<float>(showerHitsW.size()) / static_cast<float>(nHitsW);

    modularShowerMaxFoundHitRatio = std::min(std::max(std::max(foundHitRatioU, foundHitRatioV), foundHitRatioW), 1.5f);

    // Need to make fits etc...
    try
    {
        lar_content::TwoDSlidingFitResult slidingFitResultU(&cartesianPointVectorU, 1000, theGeometry->WirePitch(geo::kW));
        lar_content::TwoDSlidingFitResult slidingFitResultV(&cartesianPointVectorV, 1000, theGeometry->WirePitch(geo::kW));
        lar_content::TwoDSlidingFitResult slidingFitResultW(&cartesianPointVectorW, 1000, theGeometry->WirePitch(geo::kW));

        pandora::CartesianVector fittedShowerDirectionU(isDownstream ? slidingFitResultU.GetGlobalMinLayerDirection() : slidingFitResultU.GetGlobalMaxLayerDirection() * -1.0);
        pandora::CartesianVector fittedShowerDirectionV(isDownstream ? slidingFitResultV.GetGlobalMinLayerDirection() : slidingFitResultV.GetGlobalMaxLayerDirection() * -1.0);
        pandora::CartesianVector fittedShowerDirectionW(isDownstream ? slidingFitResultW.GetGlobalMinLayerDirection() : slidingFitResultW.GetGlobalMaxLayerDirection() * -1.0);

        // now update...
        for (auto pHit : allShowerHits)
        {
            const geo::Vector_t hitPosition(GetPandoraHitPosition(pHit, evt));
            const geo::View_t pandoraView(GetPandoraHitView(pHit));

            if (pandoraView == geo::kW)
            {
                if (std::find(showerHitsW.begin(), showerHitsW.end(), pHit) != showerHitsW.end())
                    continue;

                const geo::Vector_t displacement(hitPosition - showerStartW);
                const float l(fittedShowerDirectionW.GetDotProduct(pandora::CartesianVector(displacement.X(), displacement.Y(), displacement.Z())));
                const float t(fittedShowerDirectionW.GetCrossProduct(pandora::CartesianVector(displacement.X(), displacement.Y(), displacement.Z())).GetMagnitude());

                if ((l > 0.f) && (t < 14.0))
                {
                    showerHitsW.push_back(pHit);
                    cartesianPointVectorW.push_back(pandora::CartesianVector(hitPosition.X(), hitPosition.Y(), hitPosition.Z()));
                }
            }
            else if (pandoraView == geo::kU)
            {
                if (std::find(showerHitsU.begin(), showerHitsU.end(), pHit) != showerHitsU.end())
                    continue;

                const geo::Vector_t displacement(hitPosition - showerStartU);
                const float l(fittedShowerDirectionU.GetDotProduct(pandora::CartesianVector(displacement.X(), displacement.Y(), displacement.Z())));
                const float t(fittedShowerDirectionU.GetCrossProduct(pandora::CartesianVector(displacement.X(), displacement.Y(), displacement.Z())).GetMagnitude());

                if ((l > 0.f) && (t < 14.0))
                {
                    showerHitsU.push_back(pHit);
                    cartesianPointVectorU.push_back(pandora::CartesianVector(hitPosition.X(), hitPosition.Y(), hitPosition.Z()));
                }
            }
            else if (pandoraView == geo::kV)
            {
                if (std::find(showerHitsV.begin(), showerHitsV.end(), pHit) != showerHitsV.end())
                    continue;

                const geo::Vector_t displacement(hitPosition - showerStartV);
                const float l(fittedShowerDirectionV.GetDotProduct(pandora::CartesianVector(displacement.X(), displacement.Y(), displacement.Z())));
                const float t(fittedShowerDirectionV.GetCrossProduct(pandora::CartesianVector(displacement.X(), displacement.Y(), displacement.Z())).GetMagnitude());

                if ((l > 0.f) && (t < 14.0))
                {
                    showerHitsV.push_back(pHit);
                    cartesianPointVectorV.push_back(pandora::CartesianVector(hitPosition.X(), hitPosition.Y(), hitPosition.Z()));
                }
            }
            else 
            {
                throw cet::exception("LArPandora")
                    << "CreatePandoraHits2D - this wire view not recognised (View=" << pandoraView << ") ";
            }
        }

        std::cout << "showerHitsU.size(): " << showerHitsU.size() << std::endl;
        std::cout << "showerHitsV.size(): " << showerHitsV.size() << std::endl;
        std::cout << "showerHitsW.size(): " << showerHitsW.size() << std::endl;

        modularShowerMaxNShowerHits = std::min(std::max(std::max(showerHitsU.size(), showerHitsV.size()), showerHitsW.size()), static_cast<long unsigned int>(2000));

        foundHitRatioU = showerHitsU.size() == 0 ? 0 : static_cast<float>(showerHitsU.size()) / static_cast<float>(nHitsU);
        foundHitRatioV = showerHitsV.size() == 0 ? 0 : static_cast<float>(showerHitsV.size()) / static_cast<float>(nHitsV);
        foundHitRatioW = showerHitsW.size() == 0 ? 0 : static_cast<float>(showerHitsW.size()) / static_cast<float>(nHitsW);

        modularShowerMaxFoundHitRatio = std::min(std::max(std::max(foundHitRatioU, foundHitRatioV), foundHitRatioW), 1.5f);
        
        // this doesn't need to happen... o.O (i'm not sure if it is a mistake in the original code - it seems to work fine?)
        lar_content::TwoDSlidingFitResult updatedSlidingFitResultU(&cartesianPointVectorU, 1000, theGeometry->WirePitch(geo::kW));
        lar_content::TwoDSlidingFitResult updatedSlidingFitResultV(&cartesianPointVectorV, 1000, theGeometry->WirePitch(geo::kW));
        lar_content::TwoDSlidingFitResult updatedSlidingFitResultW(&cartesianPointVectorW, 1000, theGeometry->WirePitch(geo::kW));

        const pandora::CartesianVector &updatedFittedShowerDirectionU(isDownstream ? updatedSlidingFitResultU.GetGlobalMinLayerDirection() : 
            updatedSlidingFitResultU.GetGlobalMaxLayerDirection() * -1.0);
        const pandora::CartesianVector &updatedFittedShowerDirectionV(isDownstream ? updatedSlidingFitResultV.GetGlobalMinLayerDirection() : 
            updatedSlidingFitResultV.GetGlobalMaxLayerDirection() * -1.0);
        const pandora::CartesianVector &updatedFittedShowerDirectionW(isDownstream ? updatedSlidingFitResultW.GetGlobalMinLayerDirection() : 
            updatedSlidingFitResultW.GetGlobalMaxLayerDirection() * -1.0);
            /*
        const pandora::CartesianVector &fittedShowerStartU(isDownstream ? updatedSlidingFitResultU.GetGlobalMinLayerPosition() : 
            updatedSlidingFitResultU.GetGlobalMaxLayerPosition());
        const pandora::CartesianVector &fittedShowerStartV(isDownstream ? updatedSlidingFitResultV.GetGlobalMinLayerPosition() : 
            updatedSlidingFitResultV.GetGlobalMaxLayerPosition());
        const pandora::CartesianVector &fittedShowerStartW(isDownstream ? updatedSlidingFitResultW.GetGlobalMinLayerPosition() : 
            updatedSlidingFitResultW.GetGlobalMaxLayerPosition());
        */

        // Use the shower start here incase the shower has a gap (i think this makes the most sense)
        const float openingAngleU(GetShowerOpeningAngle(showerStartU, fittedShowerDirectionU, cartesianPointVectorU));
        const float openingAngleV(GetShowerOpeningAngle(showerStartV, fittedShowerDirectionV, cartesianPointVectorV));
        const float openingAngleW(GetShowerOpeningAngle(showerStartW, fittedShowerDirectionW, cartesianPointVectorW));

        modularShowerMaxOpeningAngle = std::min(std::max(std::max(openingAngleU, openingAngleV), openingAngleW), 20.f);

        const float scatterAngleU = std::acos(updatedFittedShowerDirectionU.GetDotProduct(connectionPathwayDirectionU)) * 180.0 /3.14;
        const float scatterAngleV = std::acos(updatedFittedShowerDirectionV.GetDotProduct(connectionPathwayDirectionV)) * 180.0 /3.14;
        const float scatterAngleW = std::acos(updatedFittedShowerDirectionW.GetDotProduct(connectionPathwayDirectionW)) * 180.0 /3.14;

        modularShowerMaxScatterAngle = std::min(std::max(std::max(scatterAngleU, scatterAngleV), scatterAngleW), 40.f);

        const pandora::CartesianVector &fittedShowerStartU(isDownstream ? slidingFitResultU.GetGlobalMinLayerPosition() : 
            slidingFitResultU.GetGlobalMaxLayerPosition());
        const pandora::CartesianVector &fittedShowerStartV(isDownstream ? slidingFitResultV.GetGlobalMinLayerPosition() : 
            slidingFitResultV.GetGlobalMaxLayerPosition());
        const pandora::CartesianVector &fittedShowerStartW(isDownstream ? slidingFitResultW.GetGlobalMinLayerPosition() : 
            slidingFitResultW.GetGlobalMaxLayerPosition());

        float nuVertexChargeAsymmetryU(-0.5f), showerStartChargeAsymmetryU(-0.5f), nuVertexChargeWeightedMeanRadialDistanceU(-10.f), showerStartMoliereRadiusU(-10.f);

        GetShowerChargeDistributionVariables(nuVertexU, connectionPathwayDirectionU, fittedShowerStartU, fittedShowerDirectionU, showerHitsU, nuVertexChargeAsymmetryU,
            showerStartChargeAsymmetryU, nuVertexChargeWeightedMeanRadialDistanceU, showerStartMoliereRadiusU, evt);

        float nuVertexChargeAsymmetryV(-0.5f), showerStartChargeAsymmetryV(-0.5f), nuVertexChargeWeightedMeanRadialDistanceV(-10.f), showerStartMoliereRadiusV(-10.f);

        GetShowerChargeDistributionVariables(nuVertexV, connectionPathwayDirectionV, fittedShowerStartV, fittedShowerDirectionV, showerHitsV, nuVertexChargeAsymmetryV,
            showerStartChargeAsymmetryV, nuVertexChargeWeightedMeanRadialDistanceV, showerStartMoliereRadiusV, evt);

        float nuVertexChargeAsymmetryW(-0.5f), showerStartChargeAsymmetryW(-0.5f), nuVertexChargeWeightedMeanRadialDistanceW(-10.f), showerStartMoliereRadiusW(-10.f);

        GetShowerChargeDistributionVariables(nuVertexW, connectionPathwayDirectionW, fittedShowerStartW, fittedShowerDirectionW, showerHitsW, nuVertexChargeAsymmetryW,
            showerStartChargeAsymmetryW, nuVertexChargeWeightedMeanRadialDistanceW, showerStartMoliereRadiusW, evt);

        modularShowerMaxNuVertexChargeAsymmetry = std::max(std::max(nuVertexChargeAsymmetryU, nuVertexChargeAsymmetryV), nuVertexChargeAsymmetryW);
        modularShowerMaxShowerStartChargeAsymmetry = std::max(std::max(showerStartChargeAsymmetryU, showerStartChargeAsymmetryV), showerStartChargeAsymmetryW);
        modularShowerMaxNuVertexChargeWeightedMeanRadialDistance = std::min(std::max(std::max(nuVertexChargeWeightedMeanRadialDistanceU, 
            nuVertexChargeWeightedMeanRadialDistanceV), nuVertexChargeWeightedMeanRadialDistanceW), 20.f);
        modularShowerMinShowerStartMoliereRadius = std::min(std::min(std::min(showerStartMoliereRadiusU, showerStartMoliereRadiusV), showerStartMoliereRadiusW), 10.f);
    }
    catch (...)
    {
        modularShowerMaxOpeningAngle = -10.f;
        modularShowerMaxScatterAngle = -10.f;
        modularShowerMaxNuVertexChargeAsymmetry = -10.f;
        modularShowerMaxShowerStartChargeAsymmetry = -10.f;
        modularShowerMaxNuVertexChargeWeightedMeanRadialDistance = -10.f;
        modularShowerMinShowerStartMoliereRadius = -10.f;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

float FDSelection::PandrizzleAlg::GetShowerOpeningAngle(const geo::Vector_t &showerStart, const pandora::CartesianVector &fittedShowerDirection, 
    pandora::CartesianPointVector &cartesianPointVector)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    pandora::CartesianVector pandoraShowerStart(showerStart.X(), showerStart.Y(), showerStart.Z());

    float openingAngle(-10.f);

    try
    {
        pandora::CartesianVector orthoAxis = fittedShowerDirection.GetCrossProduct(pandora::CartesianVector(0.f, 1.f, 0.f));

        std::map<int, float> positiveEdges, negativeEdges;

        for (const pandora::CartesianVector &position : cartesianPointVector)
        {
            const pandora::CartesianVector displacement(position - pandoraShowerStart);
            const float thisT(fittedShowerDirection.GetCrossProduct(displacement).GetMagnitude());
            const float thisL(fittedShowerDirection.GetDotProduct(displacement));
            const float orthoL(orthoAxis.GetDotProduct(displacement));

            std::map<int, float> &edgeMap(orthoL > 0.f ? positiveEdges : negativeEdges);

            const int lIndex(std::floor(thisL / 2.f));

            edgeMap[lIndex] = (edgeMap.find(lIndex) == edgeMap.end() ? thisT : std::max(edgeMap[lIndex] , thisT));
        }

        pandora::CartesianPointVector positiveEdgePositions, negativeEdgePositions;

        for (auto &entry : positiveEdges)
            positiveEdgePositions.push_back(pandora::CartesianVector(entry.second, 0.f, entry.first));

        for (auto &entry : negativeEdges)
            negativeEdgePositions.push_back(pandora::CartesianVector(entry.second, 0.f, entry.first));

        const lar_content::TwoDSlidingFitResult positiveEdgeFit(&positiveEdgePositions, 1000, theGeometry->WirePitch(geo::kW));
        const lar_content::TwoDSlidingFitResult negativeEdgeFit(&negativeEdgePositions, 1000, theGeometry->WirePitch(geo::kW));

        const pandora::CartesianVector positiveMinLayer(positiveEdgeFit.GetGlobalMinLayerPosition());
        const pandora::CartesianVector positiveMaxLayer(positiveEdgeFit.GetGlobalMaxLayerPosition());
        const pandora::CartesianVector negativeMinLayer(negativeEdgeFit.GetGlobalMinLayerPosition());
        const pandora::CartesianVector negativeMaxLayer(negativeEdgeFit.GetGlobalMaxLayerPosition());

        const pandora::CartesianVector globalPositiveMinLayer(pandoraShowerStart + (fittedShowerDirection * positiveMinLayer.GetZ()) + (orthoAxis * positiveMinLayer.GetX()));
        const pandora::CartesianVector globalPositiveMaxLayer(pandoraShowerStart + (fittedShowerDirection * positiveMaxLayer.GetZ()) + (orthoAxis * positiveMaxLayer.GetX()));
        const pandora::CartesianVector globalNegativeMinLayer(pandoraShowerStart + (fittedShowerDirection * negativeMinLayer.GetZ()) - (orthoAxis * negativeMinLayer.GetX()));
        const pandora::CartesianVector globalNegativeMaxLayer(pandoraShowerStart + (fittedShowerDirection * negativeMaxLayer.GetZ()) - (orthoAxis * negativeMaxLayer.GetX()));

        // need to make sure that they're pointing in the same direction......
        pandora::CartesianVector positiveEdgeVector((globalPositiveMaxLayer - globalPositiveMinLayer).GetUnitVector());
        positiveEdgeVector *= ((globalPositiveMinLayer - pandoraShowerStart).GetMagnitudeSquared() < (globalPositiveMaxLayer - pandoraShowerStart).GetMagnitudeSquared() ? 1.0 : -1.0);

        pandora::CartesianVector negativeEdgeVector((globalNegativeMaxLayer - globalNegativeMinLayer).GetUnitVector());
        negativeEdgeVector *= ((globalNegativeMinLayer - pandoraShowerStart).GetMagnitudeSquared() < (globalNegativeMaxLayer - pandoraShowerStart).GetMagnitudeSquared() ? 1.0 : -1.0);

        const float positiveOpeningAngle = fittedShowerDirection.GetOpeningAngle(positiveEdgeVector) * 180.f / M_PI;
        const float negativeOpeningAngle = fittedShowerDirection.GetOpeningAngle(negativeEdgeVector) * 180.f / M_PI;
        openingAngle = std::max(positiveOpeningAngle, negativeOpeningAngle);
    }
    catch(...)
    {
        openingAngle = -10.f;
    }

    return openingAngle;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::GetShowerChargeDistributionVariables(const pandora::CartesianVector nuVertexPosition, const pandora::CartesianVector connectionPathwayDirection,
    const pandora::CartesianVector &fittedShowerStart, const pandora::CartesianVector &fittedShowerDirection, std::vector<art::Ptr<recob::Hit>> &showerHits, float &nuVertexChargeAsymmetry,
    float &showerStartChargeAsymmetry, float &nuVertexChargeWeightedMeanRadialDistance, float &showerStartMoliereRadius, const art::Event& evt)
{
    pandora::CartesianVector pathwayOrthoAxis = connectionPathwayDirection.GetCrossProduct(pandora::CartesianVector(0.f, 1.f, 0.f));
    pandora::CartesianVector showerOrthoAxis = fittedShowerDirection.GetCrossProduct(pandora::CartesianVector(0.f, 1.f, 0.f));

    // Nu vertex energy asymmetry
    float totalCharge(0.f);
    nuVertexChargeAsymmetry = 0.f;

    for (auto pHit : showerHits)
    {
        const float hitCharge(std::fabs(pHit->Integral()));
        totalCharge += hitCharge;

        const geo::Vector_t hitPosition(GetPandoraHitPosition(pHit, evt));
        const pandora::CartesianVector pandoraHitPosition(hitPosition.X(), hitPosition.Y(), hitPosition.Z());
        const pandora::CartesianVector displacement(pandoraHitPosition - nuVertexPosition);
        const float thisL(pathwayOrthoAxis.GetDotProduct(displacement));

        nuVertexChargeAsymmetry += (thisL < 0.f) ? (-1.f * hitCharge) : hitCharge;
    }

    nuVertexChargeAsymmetry = (totalCharge < std::numeric_limits<float>::epsilon()) ? -0.5f : (nuVertexChargeAsymmetry / totalCharge);
    nuVertexChargeAsymmetry = std::fabs(nuVertexChargeAsymmetry);

    // Shower start energy asymmetry
    showerStartChargeAsymmetry = 0.f;

    for (auto pHit : showerHits)
    {
        const float hitCharge(std::fabs(pHit->Integral()));
        const geo::Vector_t hitPosition(GetPandoraHitPosition(pHit, evt));
        const pandora::CartesianVector pandoraHitPosition(hitPosition.X(), hitPosition.Y(), hitPosition.Z());
        const pandora::CartesianVector displacement(pandoraHitPosition - fittedShowerStart);
        const float thisL(showerOrthoAxis.GetDotProduct(displacement));

        showerStartChargeAsymmetry += (thisL < 0.f) ? (-1.f * hitCharge) : hitCharge;
    }

    showerStartChargeAsymmetry = (totalCharge < std::numeric_limits<float>::epsilon()) ? -0.5f : (showerStartChargeAsymmetry / totalCharge);
    showerStartChargeAsymmetry = std::fabs(showerStartChargeAsymmetry);

    // Mean radial distance
    nuVertexChargeWeightedMeanRadialDistance = 0.f;

    for (auto pHit : showerHits)
    {
        const float hitCharge(std::fabs(pHit->Integral()));
        const geo::Vector_t hitPosition(GetPandoraHitPosition(pHit, evt));
        const pandora::CartesianVector pandoraHitPosition(hitPosition.X(), hitPosition.Y(), hitPosition.Z());
        const pandora::CartesianVector displacement(pandoraHitPosition - nuVertexPosition);
        const float thisT(connectionPathwayDirection.GetCrossProduct(displacement).GetMagnitude());

        nuVertexChargeWeightedMeanRadialDistance += (thisT * hitCharge);
    }

    nuVertexChargeWeightedMeanRadialDistance = (totalCharge < std::numeric_limits<float>::epsilon()) ? -10.f : nuVertexChargeWeightedMeanRadialDistance / totalCharge;

    // Molliere radius
    std::sort(showerHits.begin(), showerHits.end(),
        [this, &fittedShowerStart, &fittedShowerDirection, &evt](auto pHitA, auto pHitB) -> bool {
            const geo::Vector_t hitPositionA(GetPandoraHitPosition(pHitA, evt));
            const pandora::CartesianVector pandoraHitPositionA(hitPositionA.X(), hitPositionA.Y(), hitPositionA.Z());
            const pandora::CartesianVector displacementA(pandoraHitPositionA - fittedShowerStart);

            const geo::Vector_t hitPositionB(GetPandoraHitPosition(pHitB, evt));
            const pandora::CartesianVector pandoraHitPositionB(hitPositionB.X(), hitPositionB.Y(), hitPositionB.Z());
            const pandora::CartesianVector displacementB(pandoraHitPositionB - fittedShowerStart);

            const float tA(fittedShowerDirection.GetCrossProduct(displacementA).GetMagnitude());
            const float tB(fittedShowerDirection.GetCrossProduct(displacementB).GetMagnitude());

            return tA < tB;
        });

    float showerStartRunningChargeSum(0.f);
    showerStartMoliereRadius = -10.f;

    for (auto pHit : showerHits)
    {
        const float hitCharge(std::fabs(pHit->Integral()));

        showerStartRunningChargeSum += hitCharge;

        if ((showerStartRunningChargeSum / totalCharge) > 0.9f)
        {
            const geo::Vector_t hitPosition(GetPandoraHitPosition(pHit, evt));
            const pandora::CartesianVector pandoraHitPosition(hitPosition.X(), hitPosition.Y(), hitPosition.Z());
            const pandora::CartesianVector displacement(pandoraHitPosition - fittedShowerStart);
            showerStartMoliereRadius = fittedShowerDirection.GetCrossProduct(displacement).GetMagnitude();

            break;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::FillTree(){
  if (std::abs(fVarHolder.IntVars["TruePDG"]) == 11)
  { 
      if (fVarHolder.IntVars["PFPPDG"] == 11)
          fSignalShowerTree->Fill();
  }
  else 
  { 
      if (fVarHolder.IntVars["PFPPDG"] == 11)
          fBackgroundShowerTree->Fill();
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FDSelection::PandrizzleAlg::ResetTreeVariables()
{
    for (std::map<std::string, int>::iterator mapIt = fVarHolder.IntVars.begin(); mapIt != fVarHolder.IntVars.end(); mapIt++)
        mapIt->second = -9999;

    for (std::map<std::string, float>::iterator mapIt = fVarHolder.FloatVars.begin(); mapIt != fVarHolder.FloatVars.end(); mapIt++)
        mapIt->second = -9999.f;

    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<art::Ptr<recob::Hit> > FDSelection::PandrizzleAlg::GetPFPHits(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
  std::vector<art::Ptr<recob::Hit> > pfp_hits;

  //Get the PFP handle out of the event
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  if (!(evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle)))
  {
      mf::LogWarning("PandrizzleAlg") << "Unable to find std::vector<recob::PFParticle> with module label: " << fPFParticleModuleLabel;
      return pfp_hits;
  }

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  if (!(evt.getByLabel(fClusterModuleLabel, clusterListHandle)))
  {
      mf::LogWarning("PandrizzleAlg") << "Unable to find std::vector<recob::Cluster> with module label: " << fClusterModuleLabel;
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

art::Ptr<recob::PFParticle> FDSelection::PandrizzleAlg::GetPFParticle(const art::Ptr<recob::Shower> pShower, const art::Event& evt)
{
    art::Ptr<recob::PFParticle> matchedPFParticle;

    art::Handle< std::vector<recob::Shower> > showerListHandle;
    if (!evt.getByLabel(fShowerModuleLabel, showerListHandle))
        return matchedPFParticle;

    art::FindManyP<recob::PFParticle> showerAssn(showerListHandle, evt, fShowerModuleLabel);
    std::vector<art::Ptr<recob::PFParticle>> matchedPFParticleVector = showerAssn.at(pShower.key());

    if (matchedPFParticleVector.size() != 1)
        return matchedPFParticle;

    return matchedPFParticleVector[0];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandrizzleAlg::Record FDSelection::PandrizzleAlg::RunPID(const art::Ptr<recob::Shower> pShower, const TVector3& nuVertex, const art::Event& evt) 
{
  art::Handle< std::vector<recob::PFParticle> > pfparticleListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfparticleList;
  evt.getByLabel(fPFParticleModuleLabel, pfparticleListHandle);

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector <art::Ptr<recob::Shower>> showerList;
  if (evt.getByLabel(fShowerModuleLabel, showerListHandle))
      art::fill_ptr_vector(showerList, showerListHandle);

  bool cheat(false);
  art::Handle< std::vector<recob::Shower> > cheatShowerListHandle;
  if (fCheatCharacterisation && (std::find(showerList.begin(), showerList.end(), pShower) == showerList.end()))
    cheat = true;

  //if (cheat)
  //std::cout << "Run PID, selected shower is in the cheating file" << std::endl;

  art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Shower> >{pShower}, evt, cheat ? fCheatPIDModuleLabel : fPIDModuleLabel);
  art::Ptr<anab::MVAPIDResult> mvaPIDResult(findPIDResult.at(0));
    //MVAPID vars
    if (mvaPIDResult.isAvailable())
    {
        if (isnan(mvaPIDResult->evalRatio))
            SetVar(kEvalRatio, -0.5f);
        else
            SetVar(kEvalRatio, static_cast<Float_t>(mvaPIDResult->evalRatio));

        if (isnan(mvaPIDResult->concentration))
            SetVar(kConcentration, -2.f);
        else
            SetVar(kConcentration, std::min(static_cast<Float_t>(mvaPIDResult->concentration), 50.f));
        SetVar(kCoreHaloRatio, static_cast<Float_t>(mvaPIDResult->coreHaloRatio));
        SetVar(kConicalness, std::min(static_cast<Float_t>(mvaPIDResult->conicalness), 100.f));
    }
    else 
        ReturnEmptyRecord();

    art::Ptr<recob::PFParticle> matchedPFParticle = GetPFParticle(pShower, evt);
    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn(pfparticleListHandle, evt, "pandoraSel");
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(matchedPFParticle.key());
    larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetadata[0]->GetPropertiesMap());

    if (fShiftdEdX && (pfpMetadata.size() == 1))
    {
        if (propertiesMap.find("dEdX") != propertiesMap.end())
        {
            SetVar(kdEdxBestPlane, propertiesMap.at("dEdX"));
        }
        else
        {
            if (pShower->dEdx().size() > 0)
            {
                SetVar(kdEdxBestPlane, std::max(std::min(static_cast<Float_t>(pShower->dEdx().at(pShower->best_plane())), 20.f), -2.f));
            }
            else
            {
                ReturnEmptyRecord();
            }
        }
    }
    else
    {
        if (pShower->dEdx().size() > 0)
        {
            SetVar(kdEdxBestPlane, std::max(std::min(static_cast<Float_t>(pShower->dEdx().at(pShower->best_plane())), 20.f), -2.f));
        }
        else
        {
            ReturnEmptyRecord();
        }
    }

    //Displacement
    TVector3 showerVertex = pShower->ShowerStart();

    if (fShiftDisplacement && (pfpMetadata.size() == 1))
    {
      if ((propertiesMap.find("ShowerVertexX") != propertiesMap.end()) && (propertiesMap.find("ShowerVertexY") != propertiesMap.end()) && 
          (propertiesMap.find("ShowerVertexZ") != propertiesMap.end()))
      {
          showerVertex = TVector3(propertiesMap.at("ShowerVertexX"), propertiesMap.at("ShowerVertexY"), propertiesMap.at("ShowerVertexZ"));
      }
    }

    SetVar(kDisplacement, std::min(static_cast<Float_t>((showerVertex - nuVertex).Mag()), 100.f));

    //Distance of closest approach
    double alpha((pShower->ShowerStart() - nuVertex).Dot(pShower->Direction()));
    TVector3 r(pShower->ShowerStart() + alpha*pShower->Direction());
    SetVar(kDCA, std::min(static_cast<Float_t>((r-nuVertex).Mag()), 50.f));

    //Wideness
    Float_t wideness(static_cast<Float_t>(pShower->OpenAngle()/pShower->Length()));
    if (isnan(wideness))
        SetVar(kWideness, -0.01f);
    else
        SetVar(kWideness, std::min( wideness, 0.1f));

    //Energy density
    if (pShower->Energy().size() > 0)
    {
        Float_t volume(static_cast<Float_t>((M_PI * pShower->Length() * pShower->Length() * pShower->Length() * std::tan(pShower->OpenAngle()))/3.));
        Float_t energyDensity(std::min(std::max(static_cast<Float_t>(pShower->Energy().at(2))/volume, -0.1f), 5.f));
        if (isnan(energyDensity))
            energyDensity = -0.1f;
        SetVar(kEnergyDensity, energyDensity);
    }
    else
        ReturnEmptyRecord();

    // For hit cuts
    std::vector<art::Ptr<recob::Hit>> allShowerHits(dune_ana::DUNEAnaPFParticleUtils::GetHits(matchedPFParticle, evt, "pandoraSel"));
    const int nShowerHits = allShowerHits.size(); 

    // Connecting pathway BDT variables
    if (fUseBDTVariables && (pfpMetadata.size() == 1) && (propertiesMap.find("PathwayLengthMin") != propertiesMap.end()) && (nShowerHits > fEnhancedPandrizzleHitCut))
    {
        SetVar(kPathwayLengthMin, propertiesMap.find("PathwayLengthMin") != propertiesMap.end() ? propertiesMap.at("PathwayLengthMin") : -100.0);
        SetVar(kMaxShowerStartPathwayScatteringAngle2D, propertiesMap.find("MaxShowerStartPathwayScatteringAngle2D") != propertiesMap.end() ?
          propertiesMap.at("MaxShowerStartPathwayScatteringAngle2D") : -100.0);
        SetVar(kMaxNPostShowerStartHits, propertiesMap.find("MaxNPostShowerStartHits") != propertiesMap.end() ? propertiesMap.at("MaxNPostShowerStartHits") : -100.0);
        SetVar(kMaxPostShowerStartScatterAngle, propertiesMap.find("MaxPostShowerStartScatterAngle") != propertiesMap.end() ? propertiesMap.at("MaxPostShowerStartScatterAngle") : -100.0);
        SetVar(kMaxPostShowerStartNuVertexEnergyAsymmetry, propertiesMap.find("MaxPostShowerStartNuVertexEnergyAsymmetry") != propertiesMap.end() ?
          propertiesMap.at("MaxPostShowerStartNuVertexEnergyAsymmetry") : -100.0);
        SetVar(kMaxPostShowerStartShowerStartEnergyAsymmetry, propertiesMap.find("MaxPostShowerStartShowerStartEnergyAsymmetry") != propertiesMap.end() ?
          propertiesMap.at("MaxPostShowerStartShowerStartEnergyAsymmetry") : -100.0);
        SetVar(kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance, propertiesMap.find("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance") != propertiesMap.end() ?
          propertiesMap.at("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance") : -100.0);
        SetVar(kMinPostShowerStartShowerStartMoliereRadius, propertiesMap.find("MinPostShowerStartShowerStartMoliereRadius") != propertiesMap.end() ?
          propertiesMap.at("MinPostShowerStartShowerStartMoliereRadius") : -100.0);
        SetVar(kMaxPostShowerStartOpeningAngle, propertiesMap.find("MaxPostShowerStartOpeningAngle") != propertiesMap.end() ? propertiesMap.at("MaxPostShowerStartOpeningAngle") : -100.0);
        SetVar(kMaxFoundHitRatio, propertiesMap.find("MaxFoundHitRatio") != propertiesMap.end() ? propertiesMap.at("MaxFoundHitRatio") : -100.0);
        SetVar(kMaxInitialGapSize, propertiesMap.find("MaxInitialGapSize") != propertiesMap.end() ? propertiesMap.at("MaxInitialGapSize") : -100.0);
        SetVar(kMinLargestProjectedGapSize, propertiesMap.find("MinLargestProjectedGapSize") != propertiesMap.end() ? propertiesMap.at("MinLargestProjectedGapSize") : -100.0);
        SetVar(kNViewsWithAmbiguousHits, propertiesMap.find("NViewsWithAmbiguousHits") != propertiesMap.end() ? propertiesMap.at("NViewsWithAmbiguousHits") : -100.0);
        SetVar(kAmbiguousHitMaxUnaccountedEnergy, propertiesMap.find("AmbiguousHitMaxUnaccountedEnergy") != propertiesMap.end() ? propertiesMap.at("AmbiguousHitMaxUnaccountedEnergy") : -100.0);
        SetVar(kBDTMethod, 2);

        return Record(fInputsToReader, fJamReader.EvaluateMVA("BDTG"), true);
    }
    else
    {
        SetVar(kPathwayLengthMin, -100.0);
        SetVar(kMaxShowerStartPathwayScatteringAngle2D, -100.0);
        SetVar(kMaxNPostShowerStartHits, -100.0);
        SetVar(kMaxPostShowerStartScatterAngle, -100.0);
        SetVar(kMaxPostShowerStartNuVertexEnergyAsymmetry, -100.0);
        SetVar(kMaxPostShowerStartShowerStartEnergyAsymmetry, -100.0);
        SetVar(kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance, -100.0);
        SetVar(kMinPostShowerStartShowerStartMoliereRadius, -100.0);
        SetVar(kMaxPostShowerStartOpeningAngle, -100.0);
        SetVar(kMaxFoundHitRatio, -100.0);
        SetVar(kMaxInitialGapSize, -100.0);
        SetVar(kMinLargestProjectedGapSize, -100.0);
        SetVar(kNViewsWithAmbiguousHits, -100.0);
        SetVar(kAmbiguousHitMaxUnaccountedEnergy, -100.0);
    }

    if (fUseModularShowerVariables && (nShowerHits > fModularShowerPandrizzleHitCut))
    {
        art::FindManyP<recob::Track> initialTrackAssn(showerListHandle, evt, fShowerModuleLabel);
        std::vector<art::Ptr<recob::Track>> initialTrackStub = initialTrackAssn.at(pShower.key());

        if (initialTrackStub.size() == 1)
        {
            art::Ptr<recob::Track> trackStub = initialTrackStub.at(0);

            float modularShowerPathwayLengthMin(-10.f), modularShowerPathwayKink3D(-10.f);
            GetPathwayVariables(trackStub, modularShowerPathwayLengthMin, modularShowerPathwayKink3D);

            SetVar(kModularShowerPathwayLengthMin, std::min(modularShowerPathwayLengthMin, 30.f));

            float modularShowerMaxNShowerHits(-10.f), modularShowerMaxFoundHitRatio(-10.f), modularShowerMaxOpeningAngle(-10.f), modularShowerMaxScatterAngle(-10.f), 
                modularShowerMaxNuVertexChargeAsymmetry(-10.f), modularShowerMaxShowerStartChargeAsymmetry(-10.f), modularShowerMaxNuVertexChargeWeightedMeanRadialDistance(-10.f), 
                modularShowerMinShowerStartMoliereRadius(-10.f);
            GetShowerRegionVariables(nuVertex, matchedPFParticle, trackStub, modularShowerMaxNShowerHits, modularShowerMaxFoundHitRatio, modularShowerMaxOpeningAngle, modularShowerMaxScatterAngle, 
                modularShowerMaxNuVertexChargeAsymmetry, modularShowerMaxShowerStartChargeAsymmetry, modularShowerMaxNuVertexChargeWeightedMeanRadialDistance, 
                modularShowerMinShowerStartMoliereRadius, evt);

            SetVar(kModularShowerMaxNShowerHits, std::min(modularShowerMaxNShowerHits, 2000.f));
            SetVar(kModularShowerMaxNuVertexChargeWeightedMeanRadialDistance, std::min(modularShowerMaxNuVertexChargeWeightedMeanRadialDistance, 20.f));
        }
        else
        {
            ReturnEmptyRecord(); // if no track stub found then backup...

            SetVar(kModularShowerPathwayLengthMin, -100.0);
            SetVar(kModularShowerMaxNuVertexChargeWeightedMeanRadialDistance, -100.0);
            SetVar(kModularShowerMaxNShowerHits, -100.0);
        }
    }

    SetVar(kBDTMethod, 1);
    return Record(fInputsToReader, fReader.EvaluateMVA("BDTG"), true);
}

FDSelection::PandrizzleAlg::Record FDSelection::PandrizzleAlg::ReturnEmptyRecord()
{
    Reset(fInputsToReader);
    return Record(fInputsToReader, kDefValue, false);
}

////////////////////////////////////////////////////////////////////////////////////////////////

double FDSelection::PandrizzleAlg::YZtoU(double y, double z)
{
    const double wireAngleU(0.623257);
    return (z * std::cos(wireAngleU) - y * std::sin(wireAngleU));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double FDSelection::PandrizzleAlg::YZtoV(double y, double z)
{
    const double wireAngleV(-0.623257);
    return (z * std::cos(wireAngleV) - y * std::sin(wireAngleV));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double FDSelection::PandrizzleAlg::YZtoW(double y, double z)
{
    const double wireAngleW(0.0);
    return (z * std::cos(wireAngleW) - y * std::sin(wireAngleW));
}

////////////////////////////////////////////////////////////////////////////////////////////////

const geo::Vector_t FDSelection::PandrizzleAlg::GetPandoraHitPosition(art::Ptr<recob::Hit> pHit, const art::Event& evt)
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

   const geo::WireID hit_WireID(pHit->WireID());
   const double hit_Time(pHit->PeakTime());
   const geo::View_t hit_View(pHit->View());
   const geo::CryostatID cryostatID(hit_WireID.Cryostat);

   const geo::View_t pandora_View(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));

   geo::Point_t hitXYZ = theGeometry->Cryostat(cryostatID).TPC(hit_WireID.TPC).Plane(hit_WireID.Plane).Wire(hit_WireID.Wire).GetCenter();
   const double hitY(hitXYZ.Y());
   const double hitZ(hitXYZ.Z());

   geo::Vector_t hitPosition(detProp.ConvertTicksToX(hit_Time, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat), 0.0, 0.0);

   if (pandora_View == geo::kW)
   {
       hitPosition.SetZ(YZtoW(hitY, hitZ));
   }
   else if (pandora_View == geo::kU)
   {
       hitPosition.SetZ(YZtoU(hitY, hitZ));
   }
   else if (pandora_View == geo::kV) 
   {
       hitPosition.SetZ(YZtoV(hitY, hitZ));
   }
   else 
   {
       throw cet::exception("LArPandora")
           << "CreatePandoraHits2D - this wire view not recognised (View=" << hit_View << ") ";
   }

   return hitPosition;
}

////////////////////////////////////////////////////////////////////////////////////////////////

const geo::View_t FDSelection::PandrizzleAlg::GetPandoraHitView(art::Ptr<recob::Hit> pHit)
{
    const geo::WireID hit_WireID(pHit->WireID());
    const geo::View_t hit_View(pHit->View());
    const geo::View_t pandora_View(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));

    return pandora_View;
}

////////////////////////////////////////////////////////////////////////////////////////////////
