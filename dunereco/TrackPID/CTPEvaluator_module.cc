/**
 *  @file   larpandora/DUNEAnaAnalysis/CTPEvaluator_module.cc
 *
 *  @brief  This module uses the analysis utilities to demonstrate 
 *          some of their usage. This can be used as a basis for 
 *          writing analysis code using these tools
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"

#include "dune/TrackPID/CTPHelper.h"
#include "dune/TrackPID/CTPResult.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <fstream>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace ctp
{

/**
 *  @brief  CTPEvaluator class
 */
class CTPEvaluator : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     CTPEvaluator(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~CTPEvaluator();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);

private:

  CTPHelper fConvTrackPID;
  std::string fParticleLabel;
};

DEFINE_ART_MODULE(CTPEvaluator)

} // namespace ctp

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "TRandom3.h"

#include <iostream>
#include <random>

namespace ctp
{

CTPEvaluator::CTPEvaluator(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset),
fConvTrackPID(pset.get<fhicl::ParameterSet>("ctpHelper")),
fParticleLabel(pset.get<std::string>("particleLabel"))
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

CTPEvaluator::~CTPEvaluator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPEvaluator::beginJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPEvaluator::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPEvaluator::analyze(const art::Event &evt)
{

    // Get all of the PFParticles
    const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleLabel);

    for (const art::Ptr<recob::PFParticle> &particle : particles)
    {
        // Returns a dummy value if not a track or not suitable
        CTPResult thisPID = fConvTrackPID.RunConvolutionalTrackPID(particle,evt);

        if(!thisPID.IsValid()) continue;

        std::cout << "Got a track PID for this particle" << std::endl;
        std::cout << " - Muon score   = " << thisPID.GetMuonScore() << std::endl;       
        std::cout << " - Pion score   = " << thisPID.GetPionScore() << std::endl;       
        std::cout << " - Proton score = " << thisPID.GetProtonScore() << std::endl;       

    }

}

} //namespace ctp

