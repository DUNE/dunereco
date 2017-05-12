/**
 *  @file   SkeletonAlg.cxx
 * 
 *  @brief  This algorithm accepts a list of 3D space points and attempts to find the "medial skeleton" 
 *          describing the "best" path through them. Here, "medial" is defined by the best time match
 *          between associated 2D hits along the "best" directions.
 * 
 */

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

#include "dune/RecoAlgDUNE/Cluster3DAlgs/SkeletonAlg.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{

SkeletonAlg::SkeletonAlg(fhicl::ParameterSet const &pset)
{
    reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

SkeletonAlg::~SkeletonAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SkeletonAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_minimumDeltaTicks = pset.get<double>("MinimumDeltaTicks",  0.05);
    m_maximumDeltaTicks = pset.get<double>("MaximumDeltaTicks", 10.0 );
}
    
double SkeletonAlg::FindFirstAndLastWires(std::vector<const reco::ClusterHit3D*>& hitVec,
                                          int                                     viewToCheck,
                                          int                                     referenceWire,
                                          double                                  referenceTicks,
                                          int&                                    firstWire,
                                          int&                                    lastWire) const
{
    // In the simple case the first and last wires are simply the front and back of the input vector
    firstWire = hitVec.front()->getHits()[viewToCheck]->getHit().WireID().Wire;
    lastWire  = hitVec.back()->getHits()[viewToCheck]->getHit().WireID().Wire;
    
    double maxDeltaTicks = referenceTicks - hitVec.front()->getHits()[viewToCheck]->getTimeTicks();
    double minDeltaTicks = referenceTicks - hitVec.back()->getHits()[viewToCheck]->getTimeTicks();
    
    if (minDeltaTicks > maxDeltaTicks) std::swap(maxDeltaTicks, minDeltaTicks);

    double bestDeltaTicks = 1000000.;
    
    // Can't have a gap if only one element
    if (hitVec.size() > 1)
    {
        // The issue is that there may be a gap in wires which we need to find.
        // Reset the last wire
        lastWire = firstWire;
        
        // Keep track of all the deltas
        double nextBestDeltaTicks = bestDeltaTicks;
        
        for(const auto& hitPair : hitVec)
        {
            int    curWire    = hitPair->getHits()[viewToCheck]->getHit().WireID().Wire;
            double deltaTicks = referenceTicks - hitPair->getHits()[viewToCheck]->getTimeTicks();
            
            maxDeltaTicks = std::max(maxDeltaTicks, deltaTicks);
            minDeltaTicks = std::min(minDeltaTicks, deltaTicks);
            
            if (bestDeltaTicks > fabs(deltaTicks))
            {
                // By definition, the "next best" is now the current best
                nextBestDeltaTicks = bestDeltaTicks;
                bestDeltaTicks     = fabs(deltaTicks);
            }
            else if (nextBestDeltaTicks > fabs(deltaTicks))
            {
                nextBestDeltaTicks = fabs(deltaTicks);
            }
        
            // if gap detected, take action depending on where the input wire is
            // But remember we need to be willing to accept a 1 wire gap... for efficiency
            if (fabs(curWire - lastWire) > 2000)
            {
                if (referenceWire <= lastWire) break;
            
                // Stepped over gap, reset everything
                firstWire          = curWire;
                maxDeltaTicks      = deltaTicks;
                minDeltaTicks      = deltaTicks;
                bestDeltaTicks     = fabs(deltaTicks);
                nextBestDeltaTicks = bestDeltaTicks;
            }
            
            lastWire = curWire;
        }
        
        bestDeltaTicks = nextBestDeltaTicks;
    }
    
    bestDeltaTicks = std::max(std::min(bestDeltaTicks,m_maximumDeltaTicks), m_minimumDeltaTicks);
    
    if (minDeltaTicks * maxDeltaTicks > 0. && bestDeltaTicks > 30.) bestDeltaTicks = 0.;
    
    return bestDeltaTicks;
}
    
class OrderHitsAlongWire
{
public:
    OrderHitsAlongWire(int view = 0) : m_view(view) {}
    
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        int viewToCheck = (m_view + 1) % 3;
        
        return left->getHits()[viewToCheck]->getHit().WireID().Wire < right->getHits()[viewToCheck]->getHit().WireID().Wire;
    }
private:
    int m_view;
};
    
struct OrderBestViews
{
    bool operator()(const std::pair<size_t,size_t>& left, const std::pair<size_t,size_t>& right)
    {
        return left.second < right.second;
    }
    
};
    
int SkeletonAlg::FindMedialSkeleton(reco::HitPairListPtr& hitPairList) const
{
    // Our mission is to try to find the medial skeletion of the input list of hits
    // We define that as the set of hit pairs where the pairs share the same hit in a given direction
    // and the selected medial hit is equal distance from the edges.
    // The first step in trying to do this is to build a map which relates a give 2D hit to an ordered list of hitpairs using it
    typedef std::map<const reco::ClusterHit2D*, std::vector<const reco::ClusterHit3D*> > Hit2DtoHit3DMapU;
    
    Hit2DtoHit3DMapU hit2DToHit3DMap[3];
    
    for(const auto& hitPair : hitPairList)
    {
        // Don't consider points "rejected" earlier
        if (hitPair->bitsAreSet(reco::ClusterHit3D::REJECTEDHIT)) continue;
        
        for(const auto& hit2D : hitPair->getHits())
        {
            size_t view = hit2D->getHit().View();
            
            hit2DToHit3DMap[view][hit2D].push_back(hitPair);
        }
    }
    
    // Do an explicit loop through to sort?
    for(size_t idx = 0; idx < 3; idx++)
    {
        for(auto& mapPair : hit2DToHit3DMap[idx])
        {
            size_t numHitPairs = mapPair.second.size();
            
            if (numHitPairs > 1) std::sort(mapPair.second.begin(), mapPair.second.end(), OrderHitsAlongWire(idx));
        }
    }
    
    // Keep a count of the number of skeleton points to be returned
    int nSkeletonPoints(0);
    
    // The idea is go through all the hits again and determine if they could be "skeleton" elements
    for(const auto& hitPair : hitPairList)
    {
        // Don't consider points "rejected" earlier
        if (hitPair->bitsAreSet(reco::ClusterHit3D::REJECTEDHIT)) continue;
        
        // If a hit pair we skip for now
        if (hitPair->getHits().size() < 3) continue;
        
        // Hopefully I am not confusing myself here.
        // The goal is to know, for a given 3D hit, how many other 3D hits share the 2D hits that it is made of
        // We want this to be a bit more precise in that we don't count other 3D hits if there is a gap between
        // the current hit and the one we are comparing.
        // How do we do this?
        // We consider each 2D hit comprising the 3D hit in turn... for each 2D hit we have a list of the 3D
        // hits which share this hit. Note that for the 3D hits to be unique, there must be an equal number of 2D
        // hits from the other two views. So, to count "wires" to the boundary, we can simply pick one of the other views
        
        // so the idea is to through each of the views to determine
        // a) the difference in distance to the last hit on this wire (in each direction)
        // b) the number of 3D hits using the 2D hit in the given 3D hit
        size_t numHitPairs[3]    = {0,  0,  0};    // This keeps track of the number of 3D hits a given 2D hit is associated with
        int    deltaWires[3]     = {0,  0,  0};
        double viewDeltaT[3]     = {0., 0., 0.};
        double bestDeltaTicks[3] = {0., 0., 0.};
        int    wireNumByView[3]  = {int(hitPair->getHits()[0]->getHit().WireID().Wire),
                                    int(hitPair->getHits()[1]->getHit().WireID().Wire),
                                    int(hitPair->getHits()[2]->getHit().WireID().Wire)};
        
        size_t bestViewIdx(0);
        
        // Initiate the loop over views
        for(size_t viewIdx = 0; viewIdx < 3; viewIdx++)
        {
            const reco::ClusterHit2D* hit2D          = hitPair->getHits()[viewIdx];
            double                    hit2DTimeTicks = hit2D->getTimeTicks();
            
            std::vector<const reco::ClusterHit3D*>& hitVec(hit2DToHit3DMap[viewIdx][hit2D]);
            
            numHitPairs[viewIdx] = hitVec.size();
            
            if (numHitPairs[viewIdx] > 1)
            {
                int viewToCheck = (viewIdx + 1) % 3;
                int firstWire;
                int lastWire;
                
                bestDeltaTicks[viewIdx] = FindFirstAndLastWires(hitVec,
                                                                viewToCheck,
                                                                wireNumByView[viewToCheck],
                                                                hit2DTimeTicks,
                                                                firstWire,
                                                                lastWire);
                
                int deltaFirst = wireNumByView[viewToCheck] - firstWire;
                int deltaLast  = wireNumByView[viewToCheck] - lastWire;
                
                // If either distance to the first or last wire is zero then we are an edge point
                if (deltaFirst == 0 || deltaLast == 0) hitPair->setStatusBit(reco::ClusterHit3D::EDGEHIT);
                
                deltaWires[viewIdx]  = deltaFirst + deltaLast;
                numHitPairs[viewIdx] = lastWire - firstWire + 1;
                viewDeltaT[viewIdx]  = fabs(hit2DTimeTicks - hitPair->getHits()[viewToCheck]->getTimeTicks());
            }
            // Otherwise, by definition, it is both a skeleton point and an edge point
            else hitPair->setStatusBit(reco::ClusterHit3D::SKELETONHIT | reco::ClusterHit3D::EDGEHIT);
            
            if (numHitPairs[viewIdx] < numHitPairs[bestViewIdx])
            {
                bestViewIdx = viewIdx;
            }
        }
        
        // Make sure a real index was found (which should always happen)
        if (bestViewIdx < 3 && numHitPairs[bestViewIdx] > 0)
        {
            // Make a sliding value for the width of the core region
            int maxBestWires = 0.05 * numHitPairs[bestViewIdx] + 1;
            
            maxBestWires = 5000;
            
            // Check condition for a "skeleton" point
            if (fabs(deltaWires[bestViewIdx]) < maxBestWires + 1)
            {
                
                // If exactly in the middle then we are done
                // Set a bit in the ClusterHit3D word to signify
//                if (fabs(deltaWires[bestViewIdx]) < maxBestWires || numHitPairs[bestViewIdx] < 3)
//                    hitPair->setStatusBit(reco::ClusterHit3D::SKELETONHIT);

                // Otherwise we can try to look at the other view
//                else
                {
                    // Find the next best view
                    int nextBestIdx = (bestViewIdx + 1) % 3;
                    
                    for(size_t idx = 0; idx < 3; idx++)
                    {
                        if (idx == bestViewIdx) continue;
                        
                        if (numHitPairs[idx] < numHitPairs[nextBestIdx]) nextBestIdx = idx;
                    }
                    
                    if (nextBestIdx > 2)
                    {
                        std::cout << "***** invalid next best view: " << nextBestIdx << " *******" << std::endl;
                        continue;
                    }
                    
                    if (viewDeltaT[bestViewIdx] < 1.01*bestDeltaTicks[bestViewIdx] && viewDeltaT[nextBestIdx] < 6.01*bestDeltaTicks[nextBestIdx])
                        hitPair->setStatusBit(reco::ClusterHit3D::SKELETONHIT);
                }
            }
        }
        
        // We want to keep count of "pure" skeleton points only
        if (hitPair->bitsAreSet(reco::ClusterHit3D::SKELETONHIT) && !hitPair->bitsAreSet(reco::ClusterHit3D::EDGEHIT)) nSkeletonPoints++;
    }
    
    return nSkeletonPoints;
}
    
void SkeletonAlg::GetSkeletonHits(const reco::HitPairListPtr& inputHitList, reco::HitPairListPtr& skeletonHitList) const
{
    for(const auto& hit3D : inputHitList) if (hit3D->bitsAreSet(reco::ClusterHit3D::SKELETONHIT)) skeletonHitList.emplace_back(hit3D);
    
    return;
}

void SkeletonAlg::AverageSkeletonPositions(reco::HitPairListPtr& skeletonHitList) const
{
    // NOTE: This method assumes the list being given to it is comprised of skeleton hits
    //       YMMV if you send in a complete hit collection!

    // Define a mapping between 2D hits and the 3D hits they make
    typedef std::map<const reco::ClusterHit2D*, std::vector<const reco::ClusterHit3D*> > Hit2DtoHit3DMap;
    
    // We want to keep track of this mapping by view
    Hit2DtoHit3DMap hit2DToHit3DMap[3];
    
    // Keep count of the number of skeleton hits selected
    unsigned int nSkeletonHits(0);
    
    // Execute the first loop through the hits to build the map
    for(const auto& hitPair : skeletonHitList)
    {
        // Don't consider points "rejected" earlier
        if (hitPair->bitsAreSet(reco::ClusterHit3D::REJECTEDHIT)) continue;
        
        // Count only those skeleton hits which have not been averaged
        if (!hitPair->bitsAreSet(reco::ClusterHit3D::SKELETONPOSAVE)) nSkeletonHits++;
        
        for(const auto& hit2D : hitPair->getHits())
        {
            size_t view = hit2D->getHit().View();
            
            hit2DToHit3DMap[view][hit2D].push_back(hitPair);
        }
    }
    
    // Exit early if no skeleton hits to process
    if (!nSkeletonHits) return;
    
    // The list of 3D hits associated to any given 2D hit is most useful to us if it is ordered
    // So this will loop through entries in the map and do the ordering
    for(size_t idx = 0; idx < 3; idx++)
    {
        for(auto& mapPair : hit2DToHit3DMap[idx])
        {
            size_t numHitPairs = mapPair.second.size();
            
            if (numHitPairs > 1) std::sort(mapPair.second.begin(), mapPair.second.end(), OrderHitsAlongWire(idx));
        }
    }
    
    // Ok, so the basic strategy is not entirely different from that used to build the skeleton hits in the first
    // place. The idea is to loop through all skeleton hits and then determine the average position for the hit
    // based on the wire direction with the fewest hits to the edge.
    // Currently this is a pretty simple minded approach and it would appear that some effort could be spent
    // here to improve this.
    // NOTE: the averaging being performed is ONLY in the Y-Z plane since this is where the hit ambiguity
    //       issue arises.
    reco::HitPairListPtr::iterator hitPairItr = skeletonHitList.begin();
    
    for(int bestViewVecIdx = 0; bestViewVecIdx < 2; bestViewVecIdx++)
    {
        std::list<reco::ClusterHit3D> tempHitPairList;
        reco::HitPairListPtr          tempHitPairListPtr;
        
        std::map<const reco::ClusterHit3D*, const reco::ClusterHit3D*> hit3DToHit3DMap;
        
        while(hitPairItr != skeletonHitList.end())
        {
            const reco::ClusterHit3D* hit3D = *hitPairItr++;
            
            std::vector<std::pair<size_t,size_t> > bestViewVec;
            
            for(const auto& hit2D : hit3D->getHits())
            {
                bestViewVec.push_back(std::pair<size_t,size_t>(hit2D->getHit().View(),hit2DToHit3DMap[hit2D->getHit().View()][hit2D].size()));
            }
            
            std::sort(bestViewVec.begin(), bestViewVec.end(), OrderBestViews());
            
            size_t bestViewIdx = bestViewVec[bestViewVecIdx].first;
            size_t bestViewCnt = bestViewVec[bestViewVecIdx].second;
            
            if (bestViewCnt > 5) continue;
            
            std::vector<const reco::ClusterHit3D*>& hit3DVec = hit2DToHit3DMap[bestViewIdx][hit3D->getHits()[bestViewIdx]];
            
            double avePosition[3] = {hit3D->getPosition()[0],0.,0.};
            
            for(const auto& tempHit3D : hit3DVec)
            {
                avePosition[1] += tempHit3D->getPosition()[1];
                avePosition[2] += tempHit3D->getPosition()[2];
            }
            
            avePosition[1] *= 1./double(hit3DVec.size());
            avePosition[2] *= 1./double(hit3DVec.size());
            
            tempHitPairList.emplace_back(reco::ClusterHit3D(hit3D->getID(),
                                                            hit3D->getStatusBits(),
                                                            avePosition,
                                                            hit3D->getTotalCharge(),
                                                            hit3D->getAvePeakTime(),
                                                            hit3D->getDeltaPeakTime(),
                                                            hit3D->getSigmaPeakTime(),
                                                            hit3D->getDocaToAxis(),
                                                            hit3D->getArclenToPoca(),
                                                            hit3D->getOverlapFraction(),
                                                            hit3D->getHits()));
            
            tempHitPairListPtr.push_back(&tempHitPairList.back());
            
            hit3DToHit3DMap[tempHitPairListPtr.back()] = hit3D;
        }
        
        for(const auto& pair : hit3DToHit3DMap)
        {
            pair.second->setPosition(pair.first->getPosition());
            pair.second->setStatusBit(reco::ClusterHit3D::SKELETONPOSAVE);
        }
    }
    
    return;
}
    

} // namespace lar_cluster3d
