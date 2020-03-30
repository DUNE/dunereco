/**
 *  @file   larpandora/DUNEAnaAnalysis/CTPTrackDump_module.cc
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
 *  @brief  CTPTrackDump class
 */
class CTPTrackDump : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     CTPTrackDump(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~CTPTrackDump();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);

private:

  void WriteTextFile(const art::Event &evt,const std::vector<std::vector<float>> &inputs, const std::pair<const simb::MCParticle*,float> &trueParticle,
                     const unsigned int &trackNumber, const unsigned int &nHits, const CTPResult &thisPID) const;

  CTPHelper fConvTrackPID;
  std::string fParticleLabel;
};

DEFINE_ART_MODULE(CTPTrackDump)

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

CTPTrackDump::CTPTrackDump(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset),
fConvTrackPID(pset.get<fhicl::ParameterSet>("ctpHelper")),
fParticleLabel(pset.get<std::string>("particleLabel"))
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

CTPTrackDump::~CTPTrackDump()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPTrackDump::beginJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPTrackDump::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPTrackDump::analyze(const art::Event &evt)
{

    // Get all of the PFParticles
    const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleLabel);

    unsigned int nTracks = 0;
    for (const art::Ptr<recob::PFParticle> &particle : particles)
    {
        // Returns a dummy value if not a track or not suitable
        std::vector<std::vector<float>> netInputs = fConvTrackPID.GetNetworkInputs(particle,evt);

        if(netInputs.empty()) continue;

        const std::pair<const simb::MCParticle*,float> trueParticle = fConvTrackPID.GetTrueParticle(particle,evt);
        const unsigned int nHits = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(particle,evt,fParticleLabel).size();

        // DELETE THIS BEFORE COMMITTING!!!
        CTPResult thisPID = fConvTrackPID.RunConvolutionalTrackPID(particle,evt);
        thisPID.Print();
        if(!thisPID.IsValid()) continue;

//        std::cout << "Got valid output from the network, writing output..." << std::endl;

        this->WriteTextFile(evt,netInputs,trueParticle,nTracks,nHits,thisPID);
        ++nTracks;   
    }

}

void CTPTrackDump::WriteTextFile(const art::Event &evt, const std::vector<std::vector<float>> &inputs, const std::pair<const simb::MCParticle*,float> &trueParticle,
                                 const unsigned int &trackNumber, const unsigned int &nHits, const CTPResult &thisPID) const{

       // Open our output file stream
        std::stringstream filename;
        filename << "tracks_" << evt.id().run() << "_" << evt.id().subRun() << "_" <<  evt.id().event() << "_" << trackNumber << "_" << std::floor(inputs.at(1).at(3)*1000) << ".dat";
        std::ofstream output_file(filename.str());

        for(const float dedx : inputs.at(0))
//        for(const float val : chargeVector)
        {
            output_file << dedx << "\n";
        }
        output_file << "End of dE/dx\n";

        for(const float var : inputs.at(1)){
          output_file << var << "\n";
        }

        output_file << trueParticle.first->PdgCode() << "\n";
        output_file << nHits << "\n";
        output_file << trueParticle.first->P() << "\n";
        output_file << trueParticle.second << "\n";
        output_file << thisPID.GetMuonScore() << "\n";
        output_file << thisPID.GetPionScore() << "\n";
        output_file << thisPID.GetProtonScore() << "\n";
        output_file.close();

}

} //namespace ctp

