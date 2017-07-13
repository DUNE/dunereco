/**
 *  @file   ParallelHitsSeedFinderAlg.cxx
 * 
 *  @brief  Implementation of the Seed Finder Algorithm
 *          The intent of this algorithm is to take an input list of 3D space points and from those
 *          to find candidate track start points and directions
 */

// The main include
#include "dune/RecoAlgDUNE/Cluster3DAlgs/ParallelHitsSeedFinderAlg.h"
// Framework Includes

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// ROOT includes
#include "TTree.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

ParallelHitsSeedFinderAlg::ParallelHitsSeedFinderAlg(fhicl::ParameterSet const &pset) :
     m_maxNumEdgeHits(1000),
     m_gapDistance(20.),
     m_numSeed2DHits(80),
     m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->reconfigure(pset);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    //    auto const* detectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    m_geometry = &*geometry;
    //    m_detector = detectorProperties->provider();
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParallelHitsSeedFinderAlg::~ParallelHitsSeedFinderAlg()
{
}
    
void ParallelHitsSeedFinderAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_maxNumEdgeHits = pset.get<size_t>("MaxNumEdgeHits", 1000);
    m_gapDistance    = pset.get<double>("GapDistance",     20.);
    m_numSeed2DHits  = pset.get<size_t>("NumSeed2DHits",    80);
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
}
    
bool ParallelHitsSeedFinderAlg::findTrackSeeds(reco::HitPairListPtr&      inputHitPairListPtr,
                                               reco::PrincipalComponents& inputPCA,
                                               SeedHitPairListPairVec&    seedHitPairVec) const
{
    // This routine can fail...
    bool foundGoodSeeds(false);
    
    // Make sure we are using the right pca
    reco::PrincipalComponents pca = inputPCA;
    
    if (pca.getSvdOK())
    {
        // Presume CR muons will be "downward going"...
        if (pca.getEigenVectors()[0][1] > 0.) pca.flipAxis(0);
        
        // This routine is typically called when there are LOTS of hits... so we are going to try
        // to reduce the number of operations on the full list as much as possible. However, we
        // really need the input hit list to be sorted by th input PCA so do that now
        m_pcaAlg.PCAAnalysis_calc3DDocas(inputHitPairListPtr, pca);
        
        // Use this info to sort the hits along the principle axis
        // Note that this will sort hits from the "middle" to the "outside"
        inputHitPairListPtr.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
        
        // We need to determine the number of hits to drop... and this is really driven by the number of unique
        // 2D hits we want to keep at one end of the cluster. So, make some containers and do some looping
        reco::HitPairListPtr                seedHit3DList;
        std::set<const reco::ClusterHit2D*> hit2DSet;
        
        // Look for numSeed2DHits which are continuous
        double lastArcLen = inputHitPairListPtr.front()->getArclenToPoca();
        
        for(const auto& hit3D : inputHitPairListPtr)
        {
            double arcLen = hit3D->getArclenToPoca();
            
            if (fabs(arcLen-lastArcLen) > m_gapDistance)
            {
                seedHit3DList.clear();
                hit2DSet.clear();
            }
            
            seedHit3DList.push_back(hit3D);
            
            for(const auto& hit2D : hit3D->getHits())
            {
                hit2DSet.insert(hit2D);
            }
            
            if (hit2DSet.size() > m_numSeed2DHits) break;
            
            lastArcLen = arcLen;
        }
        
        // We require a minimum number of seed hits in order to proceed
        if (hit2DSet.size() > m_numSeed2DHits)
        {
            // The above derivation determines the number of hits to keep each side of the cloud
            size_t num3DHitsToKeep = std::min(2*seedHit3DList.size(), inputHitPairListPtr.size());
            size_t numEdgeHits     = std::min(size_t(num3DHitsToKeep/2), m_maxNumEdgeHits);
        
            // Get an iterator for the hits
            reco::HitPairListPtr::iterator edgeHitItr = inputHitPairListPtr.begin();
        
            std::advance(edgeHitItr, numEdgeHits);
        
            // Make a container for copying in the edge hits and size it to hold all the hits
            reco::HitPairListPtr hit3DList;
        
            hit3DList.resize(2*numEdgeHits);
        
            // Copy the low edge hit range into the list
            reco::HitPairListPtr::iterator nextHit3DItr = std::copy(inputHitPairListPtr.begin(), edgeHitItr, hit3DList.begin());
            
            // reset the "seed hits"
            seedHit3DList.clear();
            seedHit3DList.resize(numEdgeHits);
            
            std::copy(inputHitPairListPtr.begin(), edgeHitItr, seedHit3DList.begin());
        
            // Now advance the iterator into the main container and copy the rest of the elements
            std::advance(edgeHitItr, inputHitPairListPtr.size() - 2 * numEdgeHits);
        
            std::copy(edgeHitItr, inputHitPairListPtr.end(), nextHit3DItr);
        
            // At this point we should now have trimmed out the bulk of the 3D hits from the middle
            // of the input cloud of hits. Next step is to rerun the PCA on our reduced set of hits
            reco::PrincipalComponents seedPCA;
        
            m_pcaAlg.PCAAnalysis_3D(hit3DList, seedPCA, false);
        
            if (seedPCA.getSvdOK())
            {
                // Still looking to point "down"
                if (seedPCA.getEigenVectors()[0][1] > 0.) seedPCA.flipAxis(0);
            
                // Recompute the 3D docas and arc lengths
                m_pcaAlg.PCAAnalysis_calc3DDocas(hit3DList, seedPCA);
            
                // Now sort hits in regular order
                //hit3DList.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
                seedHit3DList.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
                
                // Now translate the seedCenter by the arc len to the first hit
                double seedDir[3]   = {seedPCA.getEigenVectors()[0][0], seedPCA.getEigenVectors()[0][1], seedPCA.getEigenVectors()[0][2]};
                double seedStart[3] = {seedHit3DList.front()->getX(), seedHit3DList.front()->getY(), seedHit3DList.front()->getZ()};
                
                // Move from the first hit to the halfway point but from the first hit, not from the axis position
                double halfArcLen = 0.5 * fabs(seedHit3DList.back()->getArclenToPoca() - seedHit3DList.front()->getArclenToPoca());
                
                seedStart[0] += halfArcLen * seedDir[0];
                seedStart[1] += halfArcLen * seedDir[1];
                seedStart[2] += halfArcLen * seedDir[2];
                
                for(const auto& hit3D : seedHit3DList) hit3D->setStatusBit(0x40000000);
            
                seedHitPairVec.emplace_back(std::pair<recob::Seed, reco::HitPairListPtr>(recob::Seed(seedStart, seedDir), seedHit3DList));
                
                inputPCA = seedPCA;
                
                foundGoodSeeds = true;
            }
        }
    }
    
    return foundGoodSeeds;
}

} // namespace lar_cluster3d
