/**
 *  @file   larpandora/LArPandoraInterface/PFParticleStitcher.cxx
 *
 *  @brief  PFParticle stitcher (ultimately to be replaced by Pandora algorithms)
 */

#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/cpu_timer.h"
#include "cetlib/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "larpandoracontent/LArObjects/LArTrackPfo.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "dune/DUNEPandora/PFParticleStitcher.h"

#include <iostream>
#include <limits>

namespace lar_pandora
{

PFParticleStitcher::PFParticleStitcher(fhicl::ParameterSet const &pset) : art::EDProducer()
{
    produces< std::vector<recob::PFParticle> >();
    produces< std::vector<recob::Seed> >();
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Vertex> >();

    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::Track> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();

    produces< art::Assns<recob::Track, recob::Hit> >();

    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::~PFParticleStitcher()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring = pset.get<bool>("EnableMonitoring", false);
    m_particleLabel = pset.get<std::string>("PFParticleModuleLabel", "pandora");
    m_trackLabel = pset.get<std::string>("TrackModuleLabel", "pandora");

    m_enableStitching = pset.get<bool>("EnableStitching", true);
    m_useXcoordinate = pset.get<bool>("UseXCoordinate", true);
    m_minCosRelativeAngle = pset.get<float>("MinCosRelativeAngle", 0.966);
    m_maxLongitudinalDisplacementX = pset.get<float>("MaxLongitudinalDisplacementX", 15.f);
    m_maxTransverseDisplacement = pset.get<float>("MaxTransverseDisplacement", 5.f);
    m_relaxCosRelativeAngle = pset.get<float>("RelaxCosRelativeAngle", 0.906);
    m_relaxTransverseDisplacement = pset.get<float>("RelaxTransverseDisplacement", 2.5f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::beginJob()
{  
    if (m_enableMonitoring)
        this->InitializeMonitoring();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::endJob()
{   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::produce(art::Event &evt)
{ 
    // Run/Event Numbers
    // =================
    m_run   = evt.run();
    m_event = evt.id().event();

    mf::LogDebug("LArPandora") << " *** PFParticleStitcher::produce(...)  [Run=" << m_run << ", Event=" << m_event << "] *** " << std::endl;


    // Collect reconstructed particles and tracks
    // ==========================================
    TrackVector              inputTracks;
    VertexVector             inputVertices;
    PFParticleVector         inputParticles;
    PFParticleVector         inputParticles2;
    PFParticleVector         parentParticles;

    HitsToPFParticles        hitsToParticles;
    PFParticlesToHits        particlesToHits;
    PFParticlesToSpacePoints particlesToSpacePoints;
    PFParticlesToClusters    particlesToClusters; 
    PFParticlesToTracks      particlesToTracks;
    PFParticlesToVertices    particlesToVertices;
    PFParticleTrackMap       particleTrackMap;
    PFParticleVolumeMap      particleVolumeMap;
    PFParticleMap            particleMap;

    LArPandoraHelper::CollectTracks(evt, m_trackLabel, inputTracks, particlesToTracks);
    LArPandoraHelper::CollectVertices(evt, m_particleLabel, inputVertices, particlesToVertices);
    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, inputParticles, particlesToClusters);
    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, inputParticles2, particlesToSpacePoints);
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_particleLabel, m_particleLabel, particlesToHits, hitsToParticles);

    // Build mapping from particle to particle ID for parent/daughter navigation
    for (PFParticleVector::const_iterator pIter = inputParticles.begin(), pIterEnd = inputParticles.end(); pIter != pIterEnd; ++pIter)
    {
        const art::Ptr<recob::PFParticle> particle = *pIter;
        particleMap[particle->Self()] = particle;
    }
   
    // Select final-state particles
    for (PFParticleVector::const_iterator pIter = inputParticles.begin(), pIterEnd = inputParticles.end(); pIter != pIterEnd; ++pIter)
    {
        const art::Ptr<recob::PFParticle> particle = *pIter;

        if (LArPandoraHelper::IsFinalState(particleMap, particle))
            parentParticles.push_back(particle);
    }

    // Create particle track objects, and associate them with drift volumes
    for (PFParticleVector::const_iterator pIter = parentParticles.begin(), pIterEnd = parentParticles.end(); pIter != pIterEnd; ++pIter)
    {
        const art::Ptr<recob::PFParticle> particle = *pIter;

        if (particlesToTracks.end() == particlesToTracks.find(particle))
        continue;

        const art::Ptr<recob::Track> track = LArPandoraHelper::GetPrimaryTrack(particlesToTracks, particle);

    PFParticlesToHits::const_iterator hIter = particlesToHits.find(particle);
        if (particlesToHits.end() == hIter)
            throw cet::exception("LArPandora") << " PFParticleStitcher::produce --- Found a particle without any associated hits";

        const HitVector &hits = hIter->second;

        particleTrackMap.insert(PFParticleTrackMap::value_type(particle, PFParticleTrack(track)));
        particleVolumeMap.insert(PFParticleVolumeMap::value_type(particle, this->GetVolumeID(hits)));
    }


    // Match PFParticles across drift volumes
    // ======================================
    ParticleAssociationMatrix particleAssociationMatrix;
    PFParticleMergeMap particleMatches, particleMerges;
    PFParticleTrajectoryMap particleTrajectories;

    this->CreateParticleMatches(parentParticles, particleVolumeMap, particleTrackMap, particleAssociationMatrix);
    this->SelectParticleMatches(particleMap, particleAssociationMatrix, particleMatches);

    if (m_enableMonitoring)
        this->WriteParticleMatches(particleMatches, particleTrackMap);


    // Create lists of PFParticles and track trajectory points to be merged
    // ====================================================================
    this->SelectParticleMerges(parentParticles, particleMatches, particleMerges);
    this->OrderParticleMerges(particleMerges, particlesToTracks, particlesToHits, particleTrajectories);


    // Merge PFParticles and output to ART framework
    // =============================================
    this->ProduceArtOutput(evt, particleMap, particleMerges, particleTrajectories, 
        particlesToVertices, particlesToClusters, particlesToSpacePoints, particlesToHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("particle1", &m_particle1, "particle1/I");
    m_pRecoTree->Branch("particle2", &m_particle2, "particle2/I");
    m_pRecoTree->Branch("cosRelativeAngle", &m_cosRelativeAngle, "cosRelativeAngle/F");
    m_pRecoTree->Branch("transverseDisplacement", &m_transverseDisplacement, "transverseDisplacement/F");
    m_pRecoTree->Branch("longitudinalDisplacement", &m_longitudinalDisplacement, "longitudinalDisplacement/F");
    m_pRecoTree->Branch("longitudinalDisplacementCut", &m_longitudinalDisplacementCut, "longitudinalDisplacementCut/F");
    m_pRecoTree->Branch("deltaX", &m_deltaX, "deltaX/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::ProduceArtOutput(art::Event &evt, const PFParticleMap &particleMap, const PFParticleMergeMap &particleMerges, 
    const PFParticleTrajectoryMap &particleTrajectories, const PFParticlesToVertices &particlesToVertices,
    const PFParticlesToClusters &particlesToClusters, const PFParticlesToSpacePoints &particlesToSpacePoints, 
    const PFParticlesToHits &particlesToHits)
{
    mf::LogDebug("LArPandora") << " **** PFParticleStitcher::ProduceArtOutput(...) **** " << std::endl;

    // Set up ART outputs
    // ==================
    std::unique_ptr< std::vector<recob::PFParticle> > outputParticles( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::Seed> >       outputSeeds( new std::vector<recob::Seed> );
    std::unique_ptr< std::vector<recob::Track> >      outputTracks( new std::vector<recob::Track> );
    std::unique_ptr< std::vector<recob::Vertex> >     outputVertices( new std::vector<recob::Vertex> );

    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >    outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> >     outputParticlesToVertices( new art::Assns<recob::PFParticle, recob::Vertex> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> >      outputParticlesToTracks( new art::Assns<recob::PFParticle, recob::Track> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> >       outputParticlesToSeeds( new art::Assns<recob::PFParticle, recob::Seed> );
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> >             outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );


    // Select final-state particles
    // ============================
    PFParticleVector  parentParticles;
    PFParticleVector  daughterParticles;
    
    for (PFParticleMap::const_iterator iter = particleMap.begin(), iterEnd = particleMap.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = iter->second;

        if (LArPandoraHelper::IsFinalState(particleMap, particle))
        {
            parentParticles.push_back(particle);
        }
        else if (LArPandoraHelper::IsNeutrino(particle))
        {
            // TODO: Recover neutrinos
        }
        else
        {
            daughterParticles.push_back(particle);
        }
    }

 
    // Create a new set of particle ID codes
    // =====================================
    typedef std::map<art::Ptr<recob::PFParticle>, size_t> PFParticleIDMap;

    PFParticleIDMap particleIdMap;
    PFParticleIDMap mergeIdMap;

    int vertexCounter(0);
    int trackCounter(0);
    size_t particleCounter(0);
 
    for (PFParticleMergeMap::const_iterator iter1 = particleMerges.begin(), iterEnd1 = particleMerges.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = iter1->first;
        particleIdMap[particle] = particleCounter;

        for (PFParticleList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> mergedParticle = *iter2;

            if (particle->Self() == mergedParticle->Self())
                continue;

            mergeIdMap[mergedParticle] = particleCounter;
        }

        particleCounter++;
    }

    for (PFParticleVector::const_iterator iter1 = daughterParticles.begin(), iterEnd1 = daughterParticles.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = *iter1;
        particleIdMap[particle] = (particleCounter++);
    }


    // Create new vertices
    // ===================
    typedef std::map<art::Ptr<recob::PFParticle>, unsigned int> PFParticleElementMap; 
    typedef std::map<art::Ptr<recob::Vertex>, unsigned int> VertexElementMap;

    PFParticleElementMap particleElementMap;
    VertexElementMap vertexElementMap;

    std::vector<TVector3> vertexVector;

    for (PFParticleIDMap::const_iterator iter1 = particleIdMap.begin(), iterEnd1 = particleIdMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = iter1->first;

        // Get the old particles
        PFParticleList oldParticleList;
        PFParticleMergeMap::const_iterator iter2 = particleMerges.find(particle);

        if (particleMerges.end() == iter2)
        {
            oldParticleList.insert(particle);
        }
        else
        {
            oldParticleList.insert(iter2->second.begin(), iter2->second.end());
        }

        // Get the old vertex positions
        VertexVector oldVertexList;

        for (PFParticleList::const_iterator iter3 = oldParticleList.begin(), iterEnd3 = oldParticleList.end(); iter3 != iterEnd3; ++iter3)
        {
            PFParticlesToVertices::const_iterator iter4 = particlesToVertices.find(*iter3);

            if (particlesToVertices.end() == iter4 || iter4->second.empty())
                continue;

            if (iter4->second.size() != 1)
            throw cet::exception("LArPandora") << " PFParticleStitcher::ProduceArtOutput --- Found a particle with multiple vertices ";

            oldVertexList.push_back(*(iter4->second.begin()));
        }

        VertexVector::const_iterator vIter = (oldVertexList.size() != 1 ? oldVertexList.end() : oldVertexList.begin());

        // Get the new vertex position
        bool foundNewVertex(false);
        TVector3 newVertexPosition(0.0, 0.0, 0.0);

        if (oldVertexList.end() == vIter)
        {
            PFParticleTrajectoryMap::const_iterator iter5 = particleTrajectories.find(particle);

            if (!(particleTrajectories.end() == iter5 || iter5->second.empty()))
            {
                const pandora::CartesianVector &vtxPosition = (iter5->second.front()).GetPosition();
                const pandora::CartesianVector &endPosition = (iter5->second.back()).GetPosition();
                float bestDistanceSquared(0.5f * (endPosition - vtxPosition).GetMagnitudeSquared());

                newVertexPosition.SetXYZ(vtxPosition.GetX(), vtxPosition.GetY(), vtxPosition.GetZ());
                foundNewVertex = true;

                for (VertexVector::const_iterator iter6 = oldVertexList.begin(), iterEnd6 = oldVertexList.end(); iter6 != iterEnd6; ++iter6)
                {
                    const art::Ptr<recob::Vertex> vertex = *iter6;
                    double oldPosition[3] = {0.0, 0.0, 0.0};
                    vertex->XYZ(oldPosition);
                    const TVector3 oldVertexPosition(oldPosition[0], oldPosition[1], oldPosition[2]);
                    const float thisDistanceSquared((oldVertexPosition - newVertexPosition).Mag2());
                    if (thisDistanceSquared < bestDistanceSquared)
                    {
                        bestDistanceSquared = thisDistanceSquared;
                        vIter = iter6;
                    }
                }
            }
        }

        // Choose which vertices will be created
        bool foundNewElement(false);
        unsigned int element(vertexVector.size());
        TVector3 newElementPosition(0.0, 0.0, 0.0);

        if (oldVertexList.end() != vIter) // Recreate an old vertex
        {
            const art::Ptr<recob::Vertex> vertex = *vIter;

            VertexElementMap::const_iterator vIter1 = vertexElementMap.find(vertex);

            if (vertexElementMap.end() != vIter1)
            {
                element = vIter1->second;
            }
            else
            {
                vertexElementMap.insert(VertexElementMap::value_type(vertex, element));
                double position[3] = {0.0, 0.0, 0.0};
                vertex->XYZ(position);
                newElementPosition.SetXYZ(position[0], position[1], position[2]);
                foundNewElement = true;
            }
        }
        else if (foundNewVertex) // Create a new vertex
        {
            newElementPosition = newVertexPosition;
            foundNewElement = true;
        }
        else // There is no vertex (skip event)
        {
        continue;
        }
    
        // Record vertex positions
        particleElementMap.insert(PFParticleElementMap::value_type(particle, element));

        if (foundNewElement)
            vertexVector.push_back(newElementPosition);
    }

    // Create vertices
    for (std::vector<TVector3>::const_iterator vIter = vertexVector.begin(), vIterEnd = vertexVector.end(); vIter != vIterEnd; ++vIter)
    {
        const TVector3 position(*vIter);
        double pos[3] = { position.x(), position.y(), position.z() };
        recob::Vertex newVertex(pos, vertexCounter++);
        outputVertices->push_back(newVertex);
    }

    // Create new particles
    // ====================
    for (PFParticleIDMap::const_iterator iter1 = particleIdMap.begin(), iterEnd1 = particleIdMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> oldParticle = iter1->first;

        // Get the old particles
        PFParticleList oldParticleList;
        PFParticleMergeMap::const_iterator iter2A = particleMerges.find(oldParticle);

        if (particleMerges.end() == iter2A)
        {
            oldParticleList.insert(oldParticle);
        }
        else
        {
            oldParticleList.insert(iter2A->second.begin(), iter2A->second.end());
        }

        // Get the new particle codes
        const int newPdgCode(oldParticle->PdgCode());
        const size_t newSelfCode(iter1->second);

        // Get the new parent particle
        size_t newParentCode(recob::PFParticle::kPFParticlePrimary);

        if (!oldParticle->IsPrimary())
        {    
            PFParticleMap::const_iterator iter3A = particleMap.find(oldParticle->Parent());
            if (particleMap.end() == iter3A)
                throw cet::exception("LArPandora") << " PFParticleStitcher::ProduceArtOutput --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> oldParentParticle = iter3A->second;

            if (particleIdMap.find(oldParentParticle) != particleIdMap.end())
            {
                newParentCode = particleIdMap[oldParentParticle];
            }
            else if (mergeIdMap.find(oldParentParticle) != mergeIdMap.end())
            {
                newParentCode = mergeIdMap[oldParentParticle];
            }
        }

        // Get the new daughter particles
        std::vector<size_t> newDaughterCodes;

        const std::vector<size_t> &oldDaughterCodes = oldParticle->Daughters();

        for (std::vector<size_t>::const_iterator iter2B = oldDaughterCodes.begin(), iterEnd2B = oldDaughterCodes.end(); iter2B != iterEnd2B; ++iter2B)
        {
            int oldDaughterCode = static_cast<int>(*iter2B);

            PFParticleMap::const_iterator iter3B = particleMap.find(oldDaughterCode);
            if (particleMap.end() == iter3B)
                throw cet::exception("LArPandora") << " PFParticleStitcher::ProduceArtOutput --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> oldDaughterParticle = iter3B->second;

            if (particleIdMap.find(oldDaughterParticle) != particleIdMap.end())
            {
                newDaughterCodes.push_back(particleIdMap[oldDaughterParticle]);
            }
            else if (mergeIdMap.find(oldDaughterParticle) != mergeIdMap.end())
            {
                newDaughterCodes.push_back(mergeIdMap[oldDaughterParticle]);
            }
        }

        // Build new particle
        recob::PFParticle newParticle(newPdgCode, newSelfCode, newParentCode, newDaughterCodes);
        outputParticles->push_back(newParticle);

        // Associate new particle with clusters and space points
        HitVector hitVector;
        ClusterVector clusterVector;
        SpacePointVector spacePointVector;

        for (PFParticleList::const_iterator iter4 = oldParticleList.begin(), iterEnd4 = oldParticleList.end(); iter4 != iterEnd4; ++iter4)
        {
            const art::Ptr<recob::PFParticle> particle = *iter4;
          
            // hits
            PFParticlesToHits::const_iterator hIter1 = particlesToHits.find(particle);
            if (particlesToHits.end() == hIter1)
                continue;

            for (HitVector::const_iterator hIter2 = hIter1->second.begin(), hIterEnd2 = hIter1->second.end(); hIter2 != hIterEnd2; ++hIter2)
            {
                const art::Ptr<recob::Hit> hit = *hIter2;
                hitVector.push_back(hit);
            }

            // clusters
            PFParticlesToClusters::const_iterator cIter1 = particlesToClusters.find(particle);
            if (particlesToClusters.end() == cIter1)
                continue;

            for (ClusterVector::const_iterator cIter2 = cIter1->second.begin(), cIterEnd2 = cIter1->second.end(); cIter2 != cIterEnd2; ++cIter2)
            {
                const art::Ptr<recob::Cluster> cluster = *cIter2;
                clusterVector.push_back(cluster);
            }

            // space points
            PFParticlesToSpacePoints::const_iterator sIter1 = particlesToSpacePoints.find(particle);
            if (particlesToSpacePoints.end() == sIter1)
                continue;

            for (SpacePointVector::const_iterator sIter2 = sIter1->second.begin(), sIterEnd2 = sIter1->second.end(); sIter2 != sIterEnd2; ++sIter2)
            {
                const art::Ptr<recob::SpacePoint> spacepoint = *sIter2;
                spacePointVector.push_back(spacepoint);
            }
        }

        util::CreateAssn(*this, evt, *(outputParticles.get()), clusterVector, *(outputParticlesToClusters.get()));
        util::CreateAssn(*this, evt, *(outputParticles.get()), spacePointVector, *(outputParticlesToSpacePoints.get()));

        // Associate new particle with vertex
        PFParticleElementMap::const_iterator iter5 = particleElementMap.find(oldParticle);

        if (particleElementMap.end() != iter5)
        {
            const unsigned int vtxElement(iter5->second);

            util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputVertices.get()), *(outputParticlesToVertices.get()),
            vtxElement, vtxElement + 1);
        }

        // Build new track and new seed, and associate to new particle
        PFParticleTrajectoryMap::const_iterator iter6 = particleTrajectories.find(oldParticle);

        if (particleTrajectories.end() != iter6)
        {
            const lar_content::LArTrackStateVector &trackStateVector = iter6->second;

            if (trackStateVector.empty())
                throw cet::exception("LArPandora") << " PFParticleStitcher::ProduceArtOutput --- Found a track without any trajectory points";

            // Build track
            recob::Track newTrack(this->BuildTrack(trackCounter++, &trackStateVector)); 
            outputTracks->push_back(newTrack);  

            util::CreateAssn(*this, evt, *(outputTracks.get()), hitVector, *(outputTracksToHits.get()));
            util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputTracks.get()), *(outputParticlesToTracks.get()),
            outputTracks->size() - 1, outputTracks->size());

            // Build seed
            const lar_content::LArTrackState &trackState = trackStateVector.front();
            const pandora::CartesianVector &vtxPos = trackState.GetPosition();
            const pandora::CartesianVector &vtxDir = trackState.GetDirection();

            double pos[3]     = { vtxPos.GetX(), vtxPos.GetY(), vtxPos.GetZ() };
            double posErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors
            double dir[3]     = { vtxDir.GetX(), vtxDir.GetY(), vtxDir.GetZ() };
            double dirErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors

            recob::Seed newSeed(pos, dir, posErr, dirErr);
            outputSeeds->push_back(newSeed);

            util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputSeeds.get()), *(outputParticlesToSeeds.get()),
                outputSeeds->size() - 1, outputSeeds->size());
        }
    }

    mf::LogDebug("LArPandora") << "   Number of new particles: " << outputParticles->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new tracks: " << outputTracks->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new seeds: " << outputSeeds->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new vertices: " << outputVertices->size() << std::endl;
    
    evt.put(std::move(outputParticles));
    evt.put(std::move(outputSeeds));
    evt.put(std::move(outputTracks));
    evt.put(std::move(outputVertices));

    evt.put(std::move(outputParticlesToSpacePoints));
    evt.put(std::move(outputParticlesToClusters));
    evt.put(std::move(outputParticlesToVertices));
    evt.put(std::move(outputParticlesToTracks));
    evt.put(std::move(outputParticlesToSeeds));
    evt.put(std::move(outputTracksToHits));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int PFParticleStitcher::GetVolumeID(const HitVector &hitVector) const
{
    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;
        const geo::WireID hit_WireID(hit->WireID());
        return this->GetVolumeID(hit_WireID.Cryostat, hit_WireID.TPC);
    }

    throw cet::exception("LArPandora") << " PFParticleStitcher::GetVolumeID --- No volume ID for this collection of hits";
}

//------------------------------------------------------------------------------------------------------------------------------------------

TVector3 PFParticleStitcher::GetVolumeCenter(const HitVector &hitVector) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;  

    TVector3 sumPosition(0.0, 0.0, 0.0);
    float sumCharge(0.f);

    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;
        const geo::WireID hit_WireID(hit->WireID());
        const float thisCharge(hit->Integral());
        sumPosition += thisCharge * theGeometry->Cryostat(hit_WireID.Cryostat).TPC(hit_WireID.TPC).LocalToWorld(TVector3(0.0, 0.0, 0.0));
        sumCharge += thisCharge;
    }

    if (sumCharge > std::numeric_limits<float>::epsilon())
        return (sumPosition * (1.0 / sumCharge));

    throw cet::exception("LArPandora") << " PFParticleStitcher::GetVolumeCenter --- No volume center for this collection of hits";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::WriteParticleMatches(const PFParticleMergeMap &particleMergeMap, const PFParticleTrackMap &particleTrackMap)
{
    for (PFParticleMergeMap::const_iterator iter1 = particleMergeMap.begin(), iterEnd1 = particleMergeMap.end(); iter1 != iterEnd1; ++iter1)    
    {
        const art::Ptr<recob::PFParticle> particle1 = iter1->first;

        for (PFParticleList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> particle2 = *iter2;

            m_particle1 = particle1->Self();
            m_particle2 = particle2->Self();

            PFParticleTrackMap::const_iterator iter3 = particleTrackMap.find(particle1);
            PFParticleTrackMap::const_iterator iter4 = particleTrackMap.find(particle2);

            if (particleTrackMap.end() == iter3 || particleTrackMap.end() == iter4)
                throw cet::exception("LArPandora") << " PFParticleStitcher::WriteParticleMatches --- No trace of particle seed in seed maps";

            const PFParticleTrack seed1(iter3->second);
            const PFParticleTrack seed2(iter4->second);

            // Get closest pair of vertices
            ParticleAssociation::VertexType vertexType1(ParticleAssociation::UNDEFINED);
            ParticleAssociation::VertexType vertexType2(ParticleAssociation::UNDEFINED);

            this->GetClosestVertices(seed1, seed2, vertexType1, vertexType2);

            // Get vertex and direction at closest pair of vertices
            const pandora::CartesianVector vtx1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerPosition() :  seed1.GetOuterPosition());
            const pandora::CartesianVector dir1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerDirection() : seed1.GetOuterDirection());

            const pandora::CartesianVector vtx2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerPosition() :  seed2.GetOuterPosition());
            const pandora::CartesianVector dir2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerDirection() : seed2.GetOuterDirection());

            if (m_useXcoordinate)
            {
                this->GetImpactParameters3D(vtx1, dir1, vtx2, m_longitudinalDisplacement, m_transverseDisplacement);
                m_longitudinalDisplacementCut = this->GetLongitudinalDisplacementCut3D(dir1);
            }
            else
            {
                this->GetImpactParameters2D(vtx1, dir1, vtx2, m_longitudinalDisplacement, m_transverseDisplacement);
                m_longitudinalDisplacementCut = this->GetLongitudinalDisplacementCut2D(dir1);
            }
            
            m_cosRelativeAngle = -dir1.GetDotProduct(dir2);
            m_deltaX = this->GetDeltaX(vtx1, dir1, vtx2);

            m_pRecoTree->Fill();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::CreateParticleMatches(const PFParticleVector &particleVector, const PFParticleVolumeMap &particleVolumeMap,
    const PFParticleTrackMap &particleTrackMap, ParticleAssociationMatrix &particleAssociationMatrix) const
{
    for (PFParticleVector::const_iterator pIter1 = particleVector.begin(), pIterEnd1 = particleVector.end(); pIter1 != pIterEnd1; ++pIter1)
    {
        const art::Ptr<recob::PFParticle> particle1 = *pIter1;

        for (PFParticleVector::const_iterator pIter2 = pIter1, pIterEnd2 = particleVector.end(); pIter2 != pIterEnd2; ++pIter2)
        {
            const art::Ptr<recob::PFParticle> particle2 = *pIter2;

            if (particle1->Self() == particle2->Self())
                continue;

            if (particle1->PdgCode() != particle2->PdgCode())
                continue;

            PFParticleVolumeMap::const_iterator vIter1 = particleVolumeMap.find(particle1);
            PFParticleVolumeMap::const_iterator vIter2 = particleVolumeMap.find(particle2);

            if (particleVolumeMap.end() == vIter1 || particleVolumeMap.end() == vIter2)
                continue;

            const unsigned int volume1(vIter1->second);
            const unsigned int volume2(vIter2->second);

            // TODO: Require NEIGHBOURING volumes
            if (volume1 == volume2)
                continue;

            this->CreateParticleMatches(particle1, particle2, particleTrackMap, particleAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::CreateParticleMatches(const art::Ptr<recob::PFParticle> particle1, const art::Ptr<recob::PFParticle> particle2,
    const PFParticleTrackMap &particleTrackMap, ParticleAssociationMatrix &particleAssociationMatrix) const
{
    if (!m_enableStitching)
        return;

    PFParticleTrackMap::const_iterator iter1 = particleTrackMap.find(particle1);
    PFParticleTrackMap::const_iterator iter2 = particleTrackMap.find(particle2);

    if (particleTrackMap.end() == iter1 || particleTrackMap.end() == iter2)
        return;

    const PFParticleTrack seed1(iter1->second);
    const PFParticleTrack seed2(iter2->second);

    // Get closest pair of vertices
    ParticleAssociation::VertexType vertexType1(ParticleAssociation::UNDEFINED);
    ParticleAssociation::VertexType vertexType2(ParticleAssociation::UNDEFINED);

    try
    {
        this->GetClosestVertices(seed1, seed2, vertexType1, vertexType2);
    }
    catch (pandora::StatusCodeException& )
    {
        return;
    }
    
    // Get vertex and direction at closest pair of vertices
    const pandora::CartesianVector vtx1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerPosition() :  seed1.GetOuterPosition());
    const pandora::CartesianVector dir1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerDirection() : seed1.GetOuterDirection());

    const pandora::CartesianVector vtx2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerPosition() :  seed2.GetOuterPosition());
    const pandora::CartesianVector dir2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerDirection() : seed2.GetOuterDirection());

    // Relative Angles
    const float cosRelativeAngle(-dir1.GetDotProduct(dir2));

    if (cosRelativeAngle < m_relaxCosRelativeAngle)
        return;

    // Impact Parameters
    float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);

    if (m_useXcoordinate)
    {
        this->GetImpactParameters3D(vtx1, dir1, vtx2, rL1, rT1);
        this->GetImpactParameters3D(vtx2, dir2, vtx1, rL2, rT2);
    }
    else
    {
        try
        {
            this->GetImpactParameters2D(vtx1, dir1, vtx2, rL1, rT1);
            this->GetImpactParameters2D(vtx2, dir2, vtx1, rL2, rT2);
        }
        catch (pandora::StatusCodeException& )
        {
            return;
        }
    }

    // Selection Cuts on longitudinal Impact Parameters
    float rCutL1(-std::numeric_limits<float>::max()), rCutL2(-std::numeric_limits<float>::max());

    try
    {
        if (m_useXcoordinate)
        {
            rCutL1 = this->GetLongitudinalDisplacementCut3D(dir1);
            rCutL2 = this->GetLongitudinalDisplacementCut3D(dir2);
        }
        else
        {
            rCutL1 = this->GetLongitudinalDisplacementCut2D(dir1);
            rCutL2 = this->GetLongitudinalDisplacementCut2D(dir2);
        }
    }
    catch (pandora::StatusCodeException& )
    {
    }

    // Cut on longitudinal displacement
    if (rL1 < -1.f || rL1 > rCutL1 || rL2 < -1.f || rL2 > rCutL2)
        return;

    // Cut on transverse displacement and relative angle
    const bool minPass(std::min(rT1, rT2) < m_relaxTransverseDisplacement && cosRelativeAngle > m_relaxCosRelativeAngle);
    const bool maxPass(std::max(rT1, rT2) < m_maxTransverseDisplacement && cosRelativeAngle > m_minCosRelativeAngle);

    if (!minPass && !maxPass)
        return;

    // Store Association 
    const float particleLength1((seed1.GetInnerPosition() - seed1.GetOuterPosition()).GetMagnitudeSquared());
    const float particleLength2((seed2.GetInnerPosition() - seed2.GetOuterPosition()).GetMagnitudeSquared());

    mf::LogDebug("LArPandora") << " *** ParticleStitcher::MatchParticles(...) *** " << std::endl
                               << "     id1=" << particle1->Self() << ", id2=" << particle2->Self() << std::endl
                               << "     cosTheta=" << -dir1.GetDotProduct(dir2) << std::endl
                               << "     rL1=" << rL1 << ", rT1=" << rT1 << ", rL2=" << rL2 << ", rT2=" << rT2 << std::endl
                               << "     Length1=" << std::sqrt(particleLength1) << " Length2=" << std::sqrt(particleLength2) << std::endl;

    (void) particleAssociationMatrix[particle1].insert(ParticleAssociationMap::value_type(particle2, 
            ParticleAssociation(vertexType1, vertexType2, particleLength2)));
    (void) particleAssociationMatrix[particle2].insert(ParticleAssociationMap::value_type(particle1, 
            ParticleAssociation(vertexType2, vertexType1, particleLength1)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::SelectParticleMatches(const PFParticleMap &particleMap, const ParticleAssociationMatrix &particleAssociationMatrix,
    PFParticleMergeMap &particleMergeMap) const
{
    // First step: find the best associations A -> X and B -> Y
    // ========================================================
    ParticleAssociationMatrix intermediateAssociationMatrix;

    for (ParticleAssociationMatrix::const_iterator iter1 = particleAssociationMatrix.begin(), iterEnd1 = particleAssociationMatrix.end(); 
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> parentParticle(iter1->first);
        const ParticleAssociationMap &particleAssociationMap(iter1->second);

        int bestParticleInner(-1);
        ParticleAssociation bestAssociationInner(ParticleAssociation::UNDEFINED, ParticleAssociation::UNDEFINED, 0.f);

        int bestParticleOuter(-1);
        ParticleAssociation bestAssociationOuter(ParticleAssociation::UNDEFINED, ParticleAssociation::UNDEFINED, 0.f);

        for (ParticleAssociationMap::const_iterator iter2 = particleAssociationMap.begin(), iterEnd2 = particleAssociationMap.end(); 
            iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> daughterParticle(iter2->first);
            const ParticleAssociation &particleAssociation(iter2->second);

            // Inner associations
            if (particleAssociation.GetParent() == ParticleAssociation::INNER)
            {
                if (particleAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = particleAssociation;
                    bestParticleInner = daughterParticle->Self();
                }
            }

            // Outer associations
            if (particleAssociation.GetParent() == ParticleAssociation::OUTER)
            {
                if (particleAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = particleAssociation;
                    bestParticleOuter = daughterParticle->Self();
                }
            }
        }

        if (bestAssociationInner.GetFigureOfMerit() > std::numeric_limits<float>::epsilon())
        {
            PFParticleMap::const_iterator pIter = particleMap.find(bestParticleInner);
            if (particleMap.end() == pIter)
                throw cet::exception("LArPandora") << " PFParticleStitcher::SelectParticleMatches --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> daughterParticle = pIter->second;
            (void) intermediateAssociationMatrix[parentParticle].insert(ParticleAssociationMap::value_type(daughterParticle, bestAssociationInner));
        }

        if (bestAssociationOuter.GetFigureOfMerit() > std::numeric_limits<float>::epsilon())
        {
            PFParticleMap::const_iterator pIter = particleMap.find(bestParticleOuter);
            if (particleMap.end() == pIter)
                throw cet::exception("LArPandora") << " PFParticleStitcher::SelectParticleMatches --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> daughterParticle = pIter->second;
            (void) intermediateAssociationMatrix[parentParticle].insert(ParticleAssociationMap::value_type(daughterParticle, bestAssociationOuter));
        }
    }

    // Second step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    // =============================================================================
    for (ParticleAssociationMatrix::const_iterator iter3 = intermediateAssociationMatrix.begin(), iterEnd3 = intermediateAssociationMatrix.end(); iter3 != iterEnd3; ++iter3)
    {
        const art::Ptr<recob::PFParticle> parentParticle(iter3->first);
        const ParticleAssociationMap &parentAssociationMap(iter3->second);

        for (ParticleAssociationMap::const_iterator iter4 = parentAssociationMap.begin(), iterEnd4 = parentAssociationMap.end(); iter4 != iterEnd4; ++iter4)
        {
            const art::Ptr<recob::PFParticle> daughterParticle(iter4->first);
            const ParticleAssociation &parentToDaughterAssociation(iter4->second);

            ParticleAssociationMatrix::const_iterator iter5 = intermediateAssociationMatrix.find(daughterParticle);

            if (intermediateAssociationMatrix.end() == iter5)
                continue;

            const ParticleAssociationMap &daughterAssociationMap(iter5->second);

            ParticleAssociationMap::const_iterator iter6 = daughterAssociationMap.find(parentParticle);

            if (daughterAssociationMap.end() == iter6)
                continue;

            const ParticleAssociation &daughterToParentAssociation(iter6->second);

            if (parentToDaughterAssociation.GetParent() == daughterToParentAssociation.GetDaughter() &&
                parentToDaughterAssociation.GetDaughter() == daughterToParentAssociation.GetParent())
            {
                particleMergeMap[parentParticle].insert(daughterParticle);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::SelectParticleMerges(const PFParticleVector &particleVector, const PFParticleMergeMap &inputMergeMap,
    PFParticleMergeMap &outputMergeMap) const
{
    PFParticleList vetoList;

    for (PFParticleVector::const_iterator iter1 = particleVector.begin(), iterEnd1 = particleVector.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> seedParticle = *iter1;

        if (vetoList.count(seedParticle))
            continue;

        PFParticleList mergeList;
        this->CollectAssociatedParticles(seedParticle, seedParticle, inputMergeMap, vetoList, mergeList);

        vetoList.insert(seedParticle);
        outputMergeMap[seedParticle].insert(seedParticle);

        for (PFParticleList::const_iterator iter2 = mergeList.begin(), iterEnd2 = mergeList.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> associatedParticle = *iter2;

            if (vetoList.count(associatedParticle))
                throw cet::exception("LArPandora") << " PFParticleStitcher::SelectParticleMerges --- This particle has been counted twice!";

            vetoList.insert(associatedParticle);
            outputMergeMap[seedParticle].insert(associatedParticle);

            mf::LogDebug("LArPandora") << " *** ParticleStitcher::SelectMerges(...) *** " << std::endl
                                       << "     Seed=" << seedParticle->Self() << ", Associated=" << associatedParticle->Self() << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::OrderParticleMerges(const PFParticleMergeMap &particleMerges, const PFParticlesToTracks &particlesToTracks, 
    const PFParticlesToHits &particlesToHits, PFParticleTrajectoryMap &particleTrajectories) const
{
    // Create typedef for track trajectory
    typedef std::map< const float, lar_content::LArTrackStateVector > TrackTrajectoryMap;

    // Loop over top-level particles 
    for (PFParticleMergeMap::const_iterator iter = particleMerges.begin(), iterEnd = particleMerges.end(); iter != iterEnd; ++iter)    
    {
        const art::Ptr<recob::PFParticle> primaryParticle = iter->first;

        if (particlesToTracks.end() == particlesToTracks.find(primaryParticle))
            continue;

        // First loop: find direction of longest track
        TVector3 trackDirection(0.0, 0.0, 0.0);
        float trackLengthSquared(0.f);

        for (PFParticleList::const_iterator pIter = iter->second.begin(), pIterEnd = iter->second.end(); pIter != pIterEnd; ++pIter)
        {
            const art::Ptr<recob::PFParticle> particle = *pIter;  
            const art::Ptr<recob::Track> track = LArPandoraHelper::GetPrimaryTrack(particlesToTracks, particle);
            const float thisLengthSquared((track->End() - track->Vertex()).Mag2());

            if (thisLengthSquared > trackLengthSquared)
            {
                trackDirection = track->VertexDirection();
                trackLengthSquared = thisLengthSquared;
            }
        }

        // Second loop: build ordered list of trajectory points for each track
        TrackTrajectoryMap trajectoryMap;

        for (PFParticleList::const_iterator pIter = iter->second.begin(), pIterEnd = iter->second.end(); pIter != pIterEnd; ++pIter)
        {
            const art::Ptr<recob::PFParticle> particle = *pIter;  
            const art::Ptr<recob::Track> track = LArPandoraHelper::GetPrimaryTrack(particlesToTracks, particle);

            PFParticlesToHits::const_iterator hIter = particlesToHits.find(particle);
            if (particlesToHits.end() == hIter)
            continue;

            const HitVector &hits = hIter->second;

            if(hits.empty())
                continue;

            const TVector3 detectorPosition(this->GetVolumeCenter(hits));
            const bool isForward(trackDirection.Dot(track->End() - track->Vertex()) > 0.0);
            const float displacement(trackDirection.Dot(detectorPosition));
            const float propagation(isForward ? +1.0 : -1.0);

            // Fill vector of TrackState objects
            lar_content::LArTrackStateVector trackStateVector;

            for (unsigned int p = 0; p < track->NumberTrajectoryPoints(); ++p)
            {
                const unsigned int pEntry(isForward ? p : track->NumberTrajectoryPoints() - 1 - p);
                const TVector3 position(track->LocationAtPoint(pEntry));
                const TVector3 direction(propagation * track->DirectionAtPoint(pEntry));

                trackStateVector.push_back(lar_content::LArTrackState(
                pandora::CartesianVector(position.x(), position.y(), position.z()), 
                pandora::CartesianVector(direction.x(), direction.y(), direction.z())));
            }

            trajectoryMap.insert(TrackTrajectoryMap::value_type(displacement, trackStateVector));
        }

        // Third loop: concatenate trajectory points
        for (TrackTrajectoryMap::const_iterator tIter1 = trajectoryMap.begin(), tIterEnd1 = trajectoryMap.end(); tIter1 != tIterEnd1; ++tIter1)
        {
            for (lar_content::LArTrackStateVector::const_iterator tIter2 = tIter1->second.begin(), tIterEnd2 = tIter1->second.end(); tIter2 != tIterEnd2; ++tIter2)
            {
                    particleTrajectories[primaryParticle].push_back(*tIter2);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::CollectAssociatedParticles(art::Ptr<recob::PFParticle> seedParticle, art::Ptr<recob::PFParticle> currentParticle, 
    const PFParticleMergeMap &particleMergeMap, const PFParticleList &vetoList, PFParticleList &associatedParticleList) const
{
    if (vetoList.count(currentParticle))
        return;

    PFParticleMergeMap::const_iterator iter1 = particleMergeMap.find(currentParticle);

    if (particleMergeMap.end() == iter1)
        return;

    for (PFParticleList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        const art::Ptr<recob::PFParticle> associatedParticle = *iter2;

        if (associatedParticle->Self() == seedParticle->Self())
            continue;

        if (!associatedParticleList.insert(associatedParticle).second)
            continue;

        this->CollectAssociatedParticles(seedParticle, associatedParticle, particleMergeMap, vetoList, associatedParticleList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::GetClosestVertices(const PFParticleTrack& seed1, const PFParticleTrack &seed2, 
    ParticleAssociation::VertexType &vertexType1, ParticleAssociation::VertexType &vertexType2) const
{
    for (unsigned int inner1 = 0; inner1 < 2; ++inner1)
    {
        vertexType1 = ((0 == inner1) ? ParticleAssociation::INNER : ParticleAssociation::OUTER);
        const pandora::CartesianVector vtx1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerPosition()  : seed1.GetOuterPosition());
        const pandora::CartesianVector dir1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerDirection() : seed1.GetOuterDirection());
        const pandora::CartesianVector end1((ParticleAssociation::INNER == vertexType1) ? seed1.GetOuterPosition()  : seed1.GetInnerPosition());

        for (unsigned int inner2 = 0; inner2 < 2; ++inner2)
        {
            vertexType2 = ((0 == inner2) ? ParticleAssociation::INNER : ParticleAssociation::OUTER);
            const pandora::CartesianVector vtx2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerPosition()  : seed2.GetOuterPosition());
            const pandora::CartesianVector dir2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerDirection() : seed2.GetOuterDirection());
            const pandora::CartesianVector end2((ParticleAssociation::INNER == vertexType2) ? seed2.GetOuterPosition()  : seed2.GetInnerPosition());

            float vtxT(0.f), vtxL(0.f), endT(0.f), endL(0.f);

            // Project seed1 onto seed2
            if (m_useXcoordinate)
            {
                this->GetImpactParameters3D(vtx1, dir1, vtx2, vtxL, vtxT);
                this->GetImpactParameters3D(vtx1, dir1, end2, endL, endT);
            }
            else
                {
                try
                {
                    this->GetImpactParameters2D(vtx1, dir1, vtx2, vtxL, vtxT);
                    this->GetImpactParameters2D(vtx1, dir1, end2, endL, endT);
                }
                catch (pandora::StatusCodeException& )
                {
                    continue;
                }
            }

            if (vtxL > endL || endL < 0.f)
                continue;
            
            // Project seed2 onto seed1
            if (m_useXcoordinate)
            {
                this->GetImpactParameters3D(vtx2, dir2, vtx1, vtxL, vtxT);
                this->GetImpactParameters3D(vtx2, dir2, end1, endL, endT);
            }
            else
            {
                try
                {
                    this->GetImpactParameters2D(vtx2, dir2, vtx1, vtxL, vtxT);
                    this->GetImpactParameters2D(vtx2, dir2, end1, endL, endT);
                }
                catch (pandora::StatusCodeException& )
                {
                    continue;
                }
            }

            if (vtxL > endL || endL < 0.f)
                continue;

            // These are the closest vertices [return]
            return;
        }
    }

    // Can't find closest vertices [bail out] (ATTN: throw Pandora exception here, as code is based on pandora objects) 
    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::GetImpactParameters3D(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse) const
{
    // sign convention for longitudinal distance:
    // -positive value means initial position is downstream of target position
    transverse = initialDirection.GetCrossProduct(targetPosition-initialPosition).GetMagnitude();
    longitudinal = -initialDirection.GetDotProduct(targetPosition-initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::GetImpactParameters2D(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse) const
{
    if (std::fabs(initialDirection.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    const pandora::CartesianVector initialPosition2D(0.f, initialPosition.GetY(), initialPosition.GetZ());
    const pandora::CartesianVector initialDirection2D(0.f, initialDirection.GetY(), initialDirection.GetZ());
    const pandora::CartesianVector targetPosition2D(0.f, targetPosition.GetY(), targetPosition.GetZ());
   
    this->GetImpactParameters3D(initialPosition2D, initialDirection2D.GetUnitVector(), targetPosition2D, longitudinal, transverse);  
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::GetLongitudinalDisplacementCut2D(const pandora::CartesianVector &direction) const
{
    if (std::fabs(direction.GetX()) < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    if (std::fabs(direction.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        return m_maxLongitudinalDisplacementX;

    const pandora::CartesianVector directionX(direction.GetX(), 0.f, 0.f);
    const pandora::CartesianVector directionYZ(0.f, direction.GetY(), direction.GetZ());

    return (m_maxLongitudinalDisplacementX * directionYZ.GetMagnitude() / directionX.GetMagnitude());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::GetLongitudinalDisplacementCut3D(const pandora::CartesianVector &direction) const
{
    if (std::fabs(direction.GetX()) < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    if (std::fabs(direction.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        return m_maxLongitudinalDisplacementX;

    const pandora::CartesianVector directionX(direction.GetX(), 0.f, 0.f);
    const pandora::CartesianVector directionYZ(0.f, direction.GetY(), direction.GetZ());

    return (m_maxLongitudinalDisplacementX / directionX.GetMagnitude());    
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::GetDeltaX(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition) const
{
    if (std::fabs(initialDirection.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    const pandora::CartesianVector targetPositionYZ(0.f, targetPosition.GetY(), targetPosition.GetZ());
    const pandora::CartesianVector initialPositionYZ(0.f, initialPosition.GetY(), initialPosition.GetZ());
    const pandora::CartesianVector initialDirectionYZ(0.f, initialDirection.GetY(), initialDirection.GetZ());
    const pandora::CartesianVector initialDirectionX(initialDirection.GetX(), 0.f, 0.f);
    
    const float R(initialDirectionYZ.GetUnitVector().GetDotProduct(initialPositionYZ - targetPositionYZ) / initialDirectionYZ.GetMagnitude());

    return (-1.f * initialDirectionX.GetUnitVector().GetDotProduct(targetPosition - (initialPosition - initialDirection * R)));
}

//------------------------------------------------------------------------------------------------------------------------------------------
         
recob::Track PFParticleStitcher::BuildTrack(const int id, const lar_content::LArTrackStateVector *const pTrackStateVector)
{
    mf::LogDebug("LArPandora") << "   Building Track [" << id << "], Number of trajectory points = " << pTrackStateVector->size() << std::endl;
            
    if (pTrackStateVector->empty())
        throw cet::exception("LArPandora") << " PFParticleStitcher::BuildTrack --- No input trajectory points were provided ";
            
    // Fill list of track properties
    std::vector<TVector3>               xyz;
    std::vector<TVector3>               pxpypz;
    std::vector< std::vector<double> >  dQdx(3);
    std::vector<double>                 momentum = std::vector<double>(2, util::kBogusD);
    
    // Loop over trajectory points
    for (const lar_content::LArTrackState &nextPoint : *pTrackStateVector)
    {      
        const TVector3 trackPosition(nextPoint.GetPosition().GetX(), nextPoint.GetPosition().GetY(), nextPoint.GetPosition().GetZ());
        const TVector3 trackDirection(nextPoint.GetDirection().GetX(), nextPoint.GetDirection().GetY(), nextPoint.GetDirection().GetZ());
      
        xyz.push_back(trackPosition);
        pxpypz.push_back(trackDirection);      
    }       

    // Return a new recob::Track object (of the Bezier variety)
    return recob::Track(xyz, pxpypz, dQdx, momentum, id);
}   

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::PFParticleTrack::PFParticleTrack(const art::Ptr<recob::Track> track) : 
    m_innerPosition(0.f, 0.f, 0.f),
    m_innerDirection(0.f, 0.f, 0.f),
    m_outerPosition(0.f, 0.f, 0.f),
    m_outerDirection(0.f, 0.f, 0.f)
{
    const bool isForward(track->Vertex().z() <= track->End().z());
    const TVector3 innerPosition(isForward ? track->Vertex() : track->End()); 
    const TVector3 outerPosition(isForward ? track->End() : track->Vertex());    
    const TVector3 innerDirection(isForward ? track->VertexDirection() : -1.0 * track->EndDirection());
    const TVector3 outerDirection(isForward ? -1.0 * track->EndDirection() : track->VertexDirection());
                
    m_innerPosition  = pandora::CartesianVector(innerPosition.x(), innerPosition.y(), innerPosition.z());
    m_outerPosition  = pandora::CartesianVector(outerPosition.x(), outerPosition.y(), outerPosition.z());
    m_innerDirection = pandora::CartesianVector(innerDirection.x(), innerDirection.y(), innerDirection.z());
    m_outerDirection = pandora::CartesianVector(outerDirection.x(), outerDirection.y(), outerDirection.z());
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::PFParticleTrack::~PFParticleTrack()
{
}
  
//------------------------------------------------------------------------------------------------------------------------------------------
    
pandora::CartesianVector PFParticleStitcher::PFParticleTrack::GetInnerPosition() const 
{ 
    return m_innerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PFParticleStitcher::PFParticleTrack::GetInnerDirection() const 
{ 
    return m_innerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PFParticleStitcher::PFParticleTrack::GetOuterPosition() const 
{ 
    return m_outerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PFParticleStitcher::PFParticleTrack::GetOuterDirection() const 
{ 
    return m_outerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::ParticleAssociation::ParticleAssociation(const VertexType parent, const VertexType daughter, const float fom) :
    m_parent(parent),
    m_daughter(daughter),
    m_fom(fom)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::ParticleAssociation::VertexType PFParticleStitcher::ParticleAssociation::GetParent() const
{
    return m_parent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::ParticleAssociation::VertexType PFParticleStitcher::ParticleAssociation::GetDaughter() const
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::ParticleAssociation::GetFigureOfMerit() const
{
    return m_fom;
}

} // namespace lar_pandora
