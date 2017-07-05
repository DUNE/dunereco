/**
 *  @file   Cluster3D_module.cc
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "dune/RecoAlgDUNE/Cluster3DAlgs/DBScanAlg_DUNE35t.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

DBScanAlg_DUNE35t::DBScanAlg_DUNE35t(fhicl::ParameterSet const &pset)
{
  this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DBScanAlg_DUNE35t::~DBScanAlg_DUNE35t()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DBScanAlg_DUNE35t::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring       = pset.get<bool>("EnableMonitoring",  true);
    m_minPairPts             = pset.get<size_t>("MinPairPts",        2 );
    m_timeAdvanceGap         = pset.get<double>("TimeAdvanceGap",   50.);
    m_numSigmaPeakTime       = pset.get<double>("NumSigmaPeakTime",  5.);
    m_EpsMaxDist             = pset.get<double>("EpsilonDistanceDBScan", 5.);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    //    auto const* detectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    m_geometry = &*geometry;
    //    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    m_timeVector.resize(NUMTIMEVALUES, 0.);

    
}
    
void DBScanAlg_DUNE35t::expandCluster(EpsPairNeighborhoodMapVec&          epsNeighborhoodMapVec,
                              EpsPairNeighborhoodMapVec::iterator epsVecItr,
                              reco::HitPairListPtr&               curCluster,
                              size_t                              minPts) const
{
    // This is the main inside loop for the DBScan based clustering algorithm
    //
    // Add the current hit to the current cluster
    epsVecItr->second.first.setInCluster();
    curCluster.push_back(epsVecItr->first);
    
    // Get the list of points in this hit's epsilon neighborhood
    // Note this is a copy so we can modify locally
    EpsPairNeighborhoodList epsNeighborhoodList = epsVecItr->second.second;
    
    while(!epsNeighborhoodList.empty())
    {
        // Dereference the point so we can see in the debugger...
        const reco::ClusterHit3D* neighborPtr = epsNeighborhoodList.front();
        
        // Use that to look up the iterator in the general neighborhood map
        EpsPairNeighborhoodMapVec::iterator curPtEpsVecItr = epsNeighborhoodMapVec.begin();
        
        std::advance(curPtEpsVecItr, neighborPtr->getID());
        
        // If we've not been here before then take action...
        if (!curPtEpsVecItr->second.first.visited())
        {
            curPtEpsVecItr->second.first.setVisited();
                
            // If this epsilon neighborhood of this point is large enough then add its points to our list
            if (curPtEpsVecItr->second.first.getCount() >= minPts)
            {
                // Plan is to loop through the hits in this point's neighborhood and add them to our list provided
                // they have not already been added, or are part of a cluster, etc.
                // So, get the list of points in the neighborhood
                EpsPairNeighborhoodList& locList = curPtEpsVecItr->second.second;
                    
                // And loop through them...
                for(EpsPairNeighborhoodList::iterator hitItr = locList.begin(); hitItr != locList.end(); hitItr++)
                {
                    epsNeighborhoodList.push_back(*hitItr);
                }
            }
        }
            
        // If the point is not yet in a cluster then we now add
        if (!curPtEpsVecItr->second.first.inCluster())
        {
            curPtEpsVecItr->second.first.setInCluster();
            curCluster.push_back(curPtEpsVecItr->first);
        }
        
        epsNeighborhoodList.pop_front();
    }
    
    return;
}

bool SetPositionOrder(const std::unique_ptr<reco::ClusterHit3D>& left, const std::unique_ptr<reco::ClusterHit3D>& right)
{
    // The positions in the Y-Z plane are quantized so take advantage of that for ordering
    // First check that we are in the same "bin" in the z direction
    if (left->getHits().back()->getHit().WireID().Wire == right->getHits().back()->getHit().WireID().Wire) // These hits are "on the same w wire"
    {
        // We can use the U wires as a proxy for ordering hits in increasing Y
        // where we remember that as Y increases the U wire number decreases
        if (left->getHits().front()->getHit().WireID().Wire == right->getHits().front()->getHit().WireID().Wire)
        {
            // In the time direction we are not quantized at a large enough level to check
            return left->getX() < right->getX();
        }
        
        return left->getHits().front()->getHit().WireID().Wire > right->getHits().front()->getHit().WireID().Wire;
    }

    return left->getHits().back()->getHit().WireID().Wire < right->getHits().back()->getHit().WireID().Wire;
}


//This function sorts by Z first. It is the default, and is also used if we have particles passing through the EW
//scintillation counters on the exterior of the detector
bool SetPositionOrderZ( const std::unique_ptr<reco::ClusterHit3D>& left, const std::unique_ptr<reco::ClusterHit3D>& right ) 
{
  //Sort by Z first, then by Y, then by X.
  if( left->getZ() == right->getZ() ){
    if( left->getY() == right->getY() ){
      return left->getX() < right->getX();
    }
    return left->getY() < right->getY();
  }
  return left->getZ() < right->getZ();
}


//This function sorts by X first. It is used if we have particles passing through the NS scintillation counters on the
//exterior of the detector.
bool SetPositionOrderX( const std::unique_ptr<reco::ClusterHit3D>& left, const std::unique_ptr<reco::ClusterHit3D>& right ) 
{
  //Sort by X first, then by Z, then by Y.
  if( left->getX() == right->getX() ){
    if( left->getZ() == right->getZ() ){
      return left->getY() < right->getY();
    }
    return left->getZ() < right->getZ();
  }
  return left->getX() < right->getX();
}

//This function sorts by Y first. It is used if we have particles passing through the top scintillation counters above
//the detector.
bool SetPositionOrderY( const std::unique_ptr<reco::ClusterHit3D>& left, const std::unique_ptr<reco::ClusterHit3D>& right ) 
{
  //Sort by Y first, then by Z, then by X.
  if( left->getY() == right->getY() ){
    if( left->getZ() == right->getZ() ){
      return left->getX() < right->getX();
    }
    return left->getZ() < right->getZ();
  }
  return left->getY() < right->getY();
}


    
bool SetPeakHitPairIteratorOrder(const HitPairList::iterator& left, const HitPairList::iterator& right)
{
    return (*left)->getAvePeakTime() < (*right)->getAvePeakTime();
}
    
struct HitPairClusterOrder
{
    bool operator()(const reco::HitPairClusterMap::iterator& left, const reco::HitPairClusterMap::iterator& right)
    {
        // Watch out for the case where two clusters can have the same number of hits!
        if (left->second.size() == right->second.size())
            return left->first < right->first;
        
        return left->second.size() > right->second.size();
    }
};
    
void DBScanAlg_DUNE35t::ClusterHitsDBScan(ViewToHitVectorMap&       viewToHitVectorMap,
                                  ViewToWireToHitSetMap&    viewToWiretoHitSetMap,
                                  HitPairList&              hitPairList,
                                  reco::HitPairClusterMap&  hitPairClusterMap)
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    cet::cpu_timer theClockMakeHits;
    cet::cpu_timer theClockBuildNeighborhood;
    cet::cpu_timer theClockDBScan;
    
    m_timeVector.resize(NUMTIMEVALUES, 0.);
    
    if (m_enableMonitoring) theClockMakeHits.start();


    // The first task is to take the lists of input 2D hits (a map of view to sorted lists of 2D hits)
    // and then to build a list of 3D hits to be used in downstream processing
    size_t numHitPairs = BuildHitPairMap(viewToHitVectorMap, viewToWiretoHitSetMap, hitPairList);

    //Monitoring activity
    if (m_enableMonitoring)
    {
        theClockMakeHits.stop();
        theClockBuildNeighborhood.start();
    }

    // The container of pairs and those in each pair's epsilon neighborhood
    EpsPairNeighborhoodMapVec epsPairNeighborhoodMapVec;
    epsPairNeighborhoodMapVec.resize(numHitPairs, EpsPairNeighborhoodMapPair(0,EpsPairNeighborhoodPair()));  //<-- initialize too!
    
    // DBScan is driven of its "epsilon neighborhood". Computing adjacency within DBScan can be time
    // consuming so the idea is the prebuild the adjaceny map and then run DBScan.
    // The following call does this work
    BuildNeighborhoodMap(hitPairList, epsPairNeighborhoodMapVec);

    //Monitoring activity
    if (m_enableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        theClockDBScan.start();
    }

    // With the neighborhood built we can now form "clusters" with DBScan
    int pairClusterIdx(0);
    
    // Clear the cluster list just for something to do here...
    hitPairClusterMap.clear();
    
    // Ok, here we go!
    // We can simply iterate over the map we have just built to loop through the hits "simply"
    for(EpsPairNeighborhoodMapVec::iterator epsPairVecItr  = epsPairNeighborhoodMapVec.begin();
                                            epsPairVecItr != epsPairNeighborhoodMapVec.end();
                                            epsPairVecItr++)
    {
        // Skip the null entries (they were filtered out)
        if (!epsPairVecItr->first) continue;
        
        // If this hit has been "visited" already then skip
        if (epsPairVecItr->second.first.visited()) continue;
        
        // We are now visiting it so mark it as so
        epsPairVecItr->second.first.setVisited();
        
        // Check that density is sufficient
        if (epsPairVecItr->second.first.getCount() < m_minPairPts)
        {
            epsPairVecItr->second.first.setNoise();
        }
        else
        {
            // "Create" a new cluster and get a reference to it
            reco::HitPairListPtr& curCluster = hitPairClusterMap[pairClusterIdx++];
            
            // expand the cluster
            expandCluster(epsPairNeighborhoodMapVec, epsPairVecItr, curCluster, m_minPairPts);
        }
    }

    //Monitoring activity    
    if (m_enableMonitoring)
    {
        theClockDBScan.stop();
        
        m_timeVector[BUILDTHREEDHITS]  = theClockMakeHits.accumulated_real_time();
        m_timeVector[BUILDHITTOHITMAP] = theClockBuildNeighborhood.accumulated_real_time();
        m_timeVector[RUNDBSCAN]        = theClockDBScan.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> DBScan done, found " << hitPairClusterMap.size() << " clusters" << std::endl;
    
    return;
}


//This is second-least restrictive: U+W+V required, now with a corrected version of nearestwireid
size_t DBScanAlg_DUNE35t::BuildHitPairMap(ViewToHitVectorMap& viewToHitVectorMap, ViewToWireToHitSetMap& viewToWireToHitSetMap, HitPairList& hitPairList) const
{
  /**
   *  @brief Given input 2D hits, build out the lists of possible 3D hits
   *
   *         The current strategy: ideally all 3D hits would be comprised of a triplet of 2D hits, one from each view
     *         However, we have concern that, in particular, the v-plane may have some inefficiency which we have to be
     *         be prepared to deal with. The idea, then, is to first make the association of hits in the U and W planes
     *         and then look for the match in the V plane. In the event we don't find the match in the V plane then we 
     *         will evaluate the situation and in some instances keep the U-W pairs in order to keep efficiency high.
     */
  

    // Should we set a minimum total charge for a hit?
    size_t totalNumHits(0);
    size_t hitPairCntr(0);
    double minCharge[] = {0., 0., 0.};
    
    //-********************************************************************************
    // First job is to find all pairs of hits in different views which are "consistent"
    // Start loops over rows in HitVectorMap
    ViewToHitVectorMap::iterator mapItrU = viewToHitVectorMap.find(geo::kU);

    if (mapItrU != viewToHitVectorMap.end())
    {
        HitVector& hitVectorU = mapItrU->second;
        
        totalNumHits += hitVectorU.size();

        ViewToHitVectorMap::iterator mapItrW = viewToHitVectorMap.find(geo::kW);
        
        if (mapItrW != viewToHitVectorMap.end())
        {
            HitVector& hitVectorW = mapItrW->second;

            totalNumHits += hitVectorW.size();
            
            if (viewToHitVectorMap.find(geo::kV) != viewToHitVectorMap.end())
                totalNumHits += viewToHitVectorMap[geo::kV].size();
            
            // Take advantage that hits are sorted in "start time order"
            // Set the inner loop iterator before starting loop over outer hits
            HitVector::iterator hitVectorWStartItr = hitVectorW.begin();
            
            // Now we loop over the hits in these two layers
            for (HitVector::iterator hitItrU = hitVectorU.begin(); hitItrU != hitVectorU.end(); hitItrU++)
            {

                const reco::ClusterHit2D* hitPtrU = *hitItrU;
                
                if (hitPtrU->getHit().Integral() < minCharge[hitPtrU->getHit().View()]) continue;
                

                // This will be used in each loop so dereference the peak time here
                double hitUPeakTime = hitPtrU->getTimeTicks();
                
                // Inner loop is over hits within the time range of the outer hit
                for (HitVector::iterator hitItrW = hitVectorWStartItr; hitItrW != hitVectorW.end(); hitItrW++)
                {
                    const reco::ClusterHit2D* hitPtrW = *hitItrW;
                    
		    //Check to make sure that this hit is in the same TPC as the previous hit
		    if( hitPtrW->getHit().WireID().TPC != hitPtrU->getHit().WireID().TPC ) continue;

                    if (hitPtrW->getHit().Integral() < minCharge[hitPtrW->getHit().View()]) continue;
                    
                    // Hits are sorted in "peak time order" which we can take advantage of to try to speed
                    // the loops. Basically, we can compare the peak time for the outer loop hit against
                    // the current hit's peak time and if outside the range of interest we can take action.
                    // The range of interest is eyeballed but is meant to account for large pulse widths
                    // Start by dereferencing the inner hit peak time
                    double hitWPeakTime = hitPtrW->getTimeTicks();
                    
                    // If the outer loop's peak time is well past the current hit's then
                    // we should advance the inner loop's start iterator and keep going
                    if (hitUPeakTime >  hitWPeakTime + m_timeAdvanceGap)
                    {

		      
                        hitVectorWStartItr++;
                        continue;
                    }

                    // If the inner loop hit start time is past the end of the outer's end time, then break out loop
                    if (hitUPeakTime + m_timeAdvanceGap < hitWPeakTime) break;
                    
                    // We have a candidate hit pair combination, try to make a hit
                    reco::ClusterHit3D pair = makeHitPair(hitPtrU, hitPtrW);

                    
                    // The sign of success here is that the average peak time of the combined hits is > 0
                    // (note that when hits are combined the first window offset is accounted for)
                    if (pair.getAvePeakTime() > 0.)
                    {
		      //bool hitInVNotFound(true);


			//I need to know the plane ID to find the nearest wire in this TPC. In order to do that, I need to
			//loop through the viewToHitVector map and find the plane ID corresponding to a hit in this plane and this TPC
			size_t thePlane = 9999;
			for( HitVector::iterator hitVect_iter = viewToHitVectorMap[geo::kV].begin(); hitVect_iter != viewToHitVectorMap[geo::kV].end(); ++hitVect_iter ){
			  if( (*(hitVect_iter))->getHit().WireID().TPC == hitPtrW->getHit().WireID().TPC ){
			    thePlane = (*(hitVect_iter))->getHit().WireID().Plane;
			    break;
			  }
			}
			geo::PlaneID thePlaneID(0,hitPtrW->getHit().WireID().TPC,thePlane);
           
                        // Recover the WireID nearest in the V plane to the position of the pair
			//REL Use plane ids here to avoid TPC ambiguity introduced by just using view type
                        const geo::WireID wireIDV = NearestWireID_mod(pair.getPosition(), thePlaneID);

			//Some debug printing
			if( wireIDV.TPC != hitPtrW->getHit().WireID().TPC )
			  std::cout << "TPC of third (nearest) wire is not the same as the TPC of the first two." << std::endl;
			if( hitPtrU->getHit().WireID().TPC != hitPtrW->getHit().WireID().TPC )
			  std::cout << "TPC of wire 1 is not TPC of wire 2." << std::endl;
                            
                        // We believe the code that returns the ID is offset by one
			//Get the hits associated with the nearest wire
                        WireToHitSetMap::iterator wireToHitSetMapVItr = viewToWireToHitSetMap[geo::kV].find(wireIDV.Wire);
			
                        if (wireToHitSetMapVItr != viewToWireToHitSetMap[geo::kV].end())
                        {
			  const reco::ClusterHit2D* hit2DV = FindBestMatchingHit(wireToHitSetMapVItr->second, pair, m_numSigmaPeakTime*pair.getSigmaPeakTime());


                            // If a V hit found then it should be straightforward to make the triplet

			  if (hit2DV){
			    if( hitPtrW->getHit().WireID().TPC == hit2DV->getHit().WireID().TPC ){ 
			      
			      reco::ClusterHit3D triplet = makeHitTriplet(pair, hit2DV);
			      
			      if (triplet.getAvePeakTime() > 0.)
				{
				  
				  triplet.setID(hitPairCntr++);
				  hitPairList.emplace_back(std::unique_ptr<reco::ClusterHit3D>(new reco::ClusterHit3D(triplet)));
				  //				    goodBool = true;
				}
			      
			      //hitInVNotFound = false;
			    }
			  }
			}
			
		    }
		}
	    }
	}
    }
    hitPairList.sort(SetPositionOrderZ);
    return hitPairList.size();
    
}

 
    


//REL modified
size_t DBScanAlg_DUNE35t::BuildNeighborhoodMap(HitPairList& hitPairList,
					       EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec )const
					      
{
  size_t consistentPairsCnt = 0;
  size_t pairsChecked = 0;
  
  for (HitPairList::const_iterator pairItrO = hitPairList.begin(); pairItrO != hitPairList.end(); pairItrO++){
     
    const reco::ClusterHit3D* hitPairO   = (*pairItrO).get();
    const size_t              hitPairOID = hitPairO->getID();
    
    // Need to initialize the "first" part of the vector pseudo map
    epsPairNeighborhoodMapVec[hitPairOID].first = hitPairO;
    
    // Get reference to the list for this hit so we don't look it up inside the loop
    EpsPairNeighborhoodPair& hitPairOPair(epsPairNeighborhoodMapVec[hitPairOID].second);
    
    HitPairList::const_iterator pairItrI = pairItrO;
    
    std::map<int, std::pair<double, const reco::ClusterHit3D*> > bestTripletMap;
    
    // Get the X,Y, and Z for this triplet
    double pairO_X = hitPairO->getPosition()[0];
    double pairO_Y = hitPairO->getPosition()[1];
    double pairO_Z = hitPairO->getPosition()[2];

    // Set maximums
    double maxDist = m_EpsMaxDist; //REL 4.0; //Was 2, then 4 REL

    //Loop over second iterator
    while (++pairItrI != hitPairList.end())
      {
	const reco::ClusterHit3D* hitPairI = (*pairItrI).get();

	//OLD WAY
	// Note that the hits have been sorted by z, then y and then in x
	// Translation: hits are sorted by W wire, then in Y (increasing u, decreasing v)
	// and then in time.
	
	//NEW WAY
	// Hits have been sorted by Z position, then by Y, then by X (each in ascending order).
	// There is no reference to wire number here, since we need to consider clustering
	// across TPCs.
	
	// Get the X,Y, and Z for this triplet
	double pairI_X = hitPairI->getPosition()[0];
	double pairI_Y = hitPairI->getPosition()[1];
	double pairI_Z = hitPairI->getPosition()[2];
	
	//Sorting allows us to break this loop after the maxDist is surpassed in z
	if( pairI_Z - pairO_Z > maxDist ) break;

	//Form the distance so we can check ranges
	double distance = pow(
			      pow(pairO_X-pairI_X,2) +
			      pow(pairO_Y-pairI_Y,2) +
			      pow(pairO_Z-pairI_Z,2),0.5);
	
	// If we have passed the 3d distance then we are done with the loop
	if( distance > maxDist )continue;
	
	
	// This is the tight constraint on the hits
	if (consistentPairs(hitPairO, hitPairI))
	  {
	    double    bestBin = distance;
	    double bestDist = 10000.;
            
	    if (bestTripletMap.find(bestBin) != bestTripletMap.end())  //Look at a given case of error (a given bin) and pull out the best (shortest) distance
	      {
		bestDist = bestTripletMap[bestBin].first;
	      }
	    
	    double newDist = fabs(hitPairI->getX() - hitPairO->getX());
            
	    // This is an attempt to "prefer" triplets over pairs
	    if (hitPairI->getHits().size() < 3) newDist += 25.;
            
	    if (newDist < bestDist) bestTripletMap[bestBin] = std::pair<double, const reco::ClusterHit3D*>(newDist, hitPairI); //Select the pair of hits that is closest in time
            
	  }
	
	bestTripletMap[distance] = std::pair<double,const reco::ClusterHit3D*>(distance,hitPairI);
	
      }
    
    for(const auto& bestMapItr : bestTripletMap)
      {
	const reco::ClusterHit3D* hitPairI(bestMapItr.second.second);
        
	hitPairOPair.first.incrementCount();
	hitPairOPair.second.emplace_back(hitPairI);
	
	epsPairNeighborhoodMapVec[hitPairI->getID()].second.first.incrementCount();
	epsPairNeighborhoodMapVec[hitPairI->getID()].second.second.emplace_back(hitPairO);
        
	consistentPairsCnt++;
      }
  }
  
  mf::LogDebug("Cluster3D") << "Consistent pairs: " << consistentPairsCnt << " of " << pairsChecked << " checked." << std::endl;
  
  return consistentPairsCnt;
  
}



    
bool DBScanAlg_DUNE35t::consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2) const
{
  // Most likely this will fail
  bool consistent(false);
  
  // Strategy: We assume that our mission is to check only in the time direction
  //           The idea is to check each 2D hit individually but to account for
  //           skipped wires if necessary. The criterion applied here is that 2
  //           of 3 times in the views must match the "tight" tolerance with the
    //           largest deltaT must still in "in range"
    //           Note that the input values can be either pairs of hits in U-W plane
    //           or a triplet of U-V-W 2D hits
    //           Start with U plane
    std::vector<bool> matchedHits = {false, false, false};
    double            maxDeltaT(0.);
    double            maxAllowedDeltaT(0.);
    
    const recob::Hit& pair1UHit(pair1->getHits().front()->getHit());
    const recob::Hit& pair2UHit(pair2->getHits().front()->getHit());
    
    int    maxDeltaWires   = pair1->getHits().size() < pair2->getHits().size()
                           ? int(pair1->getHits().size())
                           : int(pair2->getHits().size());
    
    int    pair1UWire      = pair1UHit.WireID().Wire;
    double pair1UPeakTime  = pair1UHit.PeakTime();
    double pair1ULimit     = pair1UHit.RMS() * 2.;
    
    int    pair2UWire      = pair2UHit.WireID().Wire;
    double pair2UPeakTime  = pair2UHit.PeakTime();
    double pair2ULimit     = pair2UHit.RMS() * 2.;
    
    int    deltaUWire      = fabs(pair1UWire - pair2UWire);
    double limitTotalU     = pair1ULimit + pair2ULimit;
    
    // Do we need to constrain this?
    limitTotalU = std::min(limitTotalU, 100.);
    
    // Special conditions
    if      (deltaUWire == 0) limitTotalU *= 1.5; //3.;
    else if (deltaUWire >  1) limitTotalU *= 3.;  //2.;
    
    // If we are happy with U wires then go to the next step
    if (deltaUWire < maxDeltaWires && fabs(pair1UPeakTime - pair2UPeakTime) < limitTotalU) matchedHits[0] = true;

    maxDeltaT        = std::max(maxDeltaT, fabs(pair1UPeakTime - pair2UPeakTime));
    maxAllowedDeltaT = std::max(maxAllowedDeltaT, limitTotalU);
    
    if (deltaUWire < maxDeltaWires)
    {
        // Now do the W plane which will always be the "back" element
        const recob::Hit& pair1WHit(pair1->getHits().back()->getHit());
        const recob::Hit& pair2WHit(pair2->getHits().back()->getHit());
        
        int    pair1WWire      = pair1WHit.WireID().Wire;
        double pair1WPeakTime  = pair1WHit.PeakTime();
        double pair1WLimit     = pair1WHit.RMS() * 2.;
        
        int    pair2WWire      = pair2WHit.WireID().Wire;
        double pair2WPeakTime  = pair2WHit.PeakTime();
        double pair2WLimit     = pair2WHit.RMS() * 2.;
        
        int    deltaWWire      = fabs(pair1WWire - pair2WWire);
        double limitTotalW     = pair1WLimit + pair2WLimit;
        
        // Do we need to constrain this?
        limitTotalW = std::min(limitTotalW, 100.);
        
        // Special conditions
        if      (deltaWWire == 0) limitTotalW *= 1.5; //3.;
        else if (deltaWWire >  1) limitTotalW *= 3.0; //2.;
        
        // If we are happy with U wires then go to the next step
        if (deltaWWire < maxDeltaWires && fabs(pair1WPeakTime - pair2WPeakTime) < limitTotalW) matchedHits[1] = true;

        maxDeltaT        = std::max(maxDeltaT, fabs(pair1WPeakTime - pair2WPeakTime));
        maxAllowedDeltaT = std::max(maxAllowedDeltaT, limitTotalW);
        
        if (deltaWWire < maxDeltaWires)
        {
            // If we are triplet then keep going
            if (pair1->getHits().size() > 2 && pair2->getHits().size() > 2)
            {
                // Now do the W plane which will always be the "back" element
                const recob::Hit& pair1VHit(pair1->getHits()[1]->getHit());
                const recob::Hit& pair2VHit(pair2->getHits()[1]->getHit());
                
                int    pair1VWire      = pair1VHit.WireID().Wire;
                double pair1VPeakTime  = pair1VHit.PeakTime();
                double pair1VLimit     = pair1VHit.RMS() * 2.;
                
                int    pair2VWire      = pair2VHit.WireID().Wire;
                double pair2VPeakTime  = pair2VHit.PeakTime();
                double pair2VLimit     = pair2VHit.RMS() * 2.;
                
                int    deltaVWire      = fabs(pair1VWire - pair2VWire);
                double limitTotalV     = pair1VLimit + pair2VLimit;
                
                // Do we need to constrain this?
                limitTotalV = std::min(limitTotalV, 100.);
                
                // Special conditions
                if      (deltaVWire == 0) limitTotalV *= 2.; //4.;
                else if (deltaVWire >  1) limitTotalV *= 4.; //2.;
                
                // If we are happy with U wires then go to the next step
                if (deltaVWire < maxDeltaWires)
                {
                    if (fabs(pair1VPeakTime - pair2VPeakTime) < limitTotalV) matchedHits[2] = true;
                    
                    maxDeltaT        = std::max(maxDeltaT, fabs(pair1VPeakTime - pair2VPeakTime));
                    maxAllowedDeltaT = std::max(maxAllowedDeltaT, limitTotalV);
                }
                else
                {
                    // If wires are not in range then we zap the total vector
                    matchedHits[0] = false;
                    matchedHits[1] = false;
                }
            }
        }
    }
    
    // Now count the number of matches
    int numMatches = std::count(matchedHits.begin(), matchedHits.end(), true);
    
    if      (numMatches > 2) consistent = true;
    else if (numMatches > 1)
    {
        maxAllowedDeltaT = std::min(5.*maxAllowedDeltaT, 50.);

        if (maxDeltaT < maxAllowedDeltaT) consistent = true;
    }
    
    return consistent;
}


reco::ClusterHit3D DBScanAlg_DUNE35t::makeHitPair(const reco::ClusterHit2D* hit1,
						  const reco::ClusterHit2D* hit2,
						  double                    hitWidthSclFctr,
						  size_t                    hitPairCntr) const
{
  // No matter what happens, we will return a HitPair object
    unsigned statusBits(0);
    double   position[] = {0.,0.,0.};
    double   totalCharge(0.);
    double   avePeakTime(0.);
    double   deltaPeakTime(0.);
    double   sigmaPeakTime(0.);
    double   overlapFraction(0.);
    double   hitDocaToAxis(9999.);
    double   hitArclenToPoca(0.);
    
    std::vector<const reco::ClusterHit2D*> hitVector;
    
    hitVector.push_back(hit1);
    hitVector.push_back(hit2);
    
    // We assume in this routine that we are looking at hits in different views
    // The first mission is to check that the wires intersect
    const geo::WireID& hit1WireID = hit1->getHit().WireID();
    const geo::WireID& hit2WireID = hit2->getHit().WireID();
    
    geo::WireIDIntersection widIntersect;
            
    if (m_geometry->WireIDsIntersect(hit1WireID, hit2WireID, widIntersect))
    {
        // Wires intersect so now we can check the timing
        // We'll need the plane offsets to proceed, these can be inferred from the "time ticks" stored in the 2D hit
        double hit1Offset = hit1->getHit().PeakTime() - hit1->getTimeTicks();
        double hit2Offset = hit2->getHit().PeakTime() - hit2->getTimeTicks();
        
        // Basically, require that the hit times "overlap"
        // Check here is that they are inconsistent
        double hit1Peak  = hit1->getHit().PeakTime()  - hit1Offset;
        double hit1Sigma = hit1->getHit().RMS() * 2. / 2.355;
        
        double hit2Peak  = hit2->getHit().PeakTime()  - hit2Offset;
        double hit2Sigma = hit2->getHit().RMS() * 2. / 2.355;
        
//        double hit1Width = 2.*0.5*hit1->getHit().RMS() * 2.;
//        double hit2Width = 2.*0.5*hit2->getHit().RMS() * 2.;
        double hit1Width = hitWidthSclFctr * hit1->getHit().RMS() * 2.;
        double hit2Width = hitWidthSclFctr * hit2->getHit().RMS() * 2.;
        
        // Check hit times are consistent
        if (fabs(hit1Peak - hit2Peak) <= (hit1Width + hit2Width))
        {
            double maxUpper        = std::min(hit1Peak+hit1Width,hit2Peak+hit2Width);
            double minLower        = std::max(hit1Peak-hit1Width,hit2Peak-hit2Width);
            double overlap         = maxUpper - minLower;
            
            overlapFraction = 0.5 * overlap / std::min(hit1Width,hit2Width);
            
            if (overlapFraction > 0.)
            {
                avePeakTime     = 0.5 * (hit1Peak + hit2Peak);
                deltaPeakTime   = fabs(hit1Peak - hit2Peak);
                sigmaPeakTime   = sqrt(hit1Sigma*hit1Sigma + hit2Sigma*hit2Sigma);
                totalCharge     = hit1->getHit().Integral() + hit2->getHit().Integral();
            
                sigmaPeakTime   = std::max(sigmaPeakTime, 0.1);
            
                double xPositionHit1(hit1->getXPosition());
                double xPositionHit2(hit2->getXPosition());
                double xPosition = 0.5 * (xPositionHit1 + xPositionHit2);
            
                position[0] = xPosition;
                position[1] = widIntersect.y;
                position[2] = widIntersect.z;
            
                // If to here then we need to sort out the hit pair code telling what views are used
                statusBits = (unsigned(hit1->getHit().View()) + 1) * unsigned(hit2->getHit().View());
            }
        }
    }
    
    // Create the 3D cluster hit
    reco::ClusterHit3D hitPair(hitPairCntr,
                               statusBits,
                               position,
                               totalCharge,
                               avePeakTime,
                               deltaPeakTime,
                               sigmaPeakTime,
                               overlapFraction,
                               hitDocaToAxis,
                               hitArclenToPoca,
                               hitVector);
    // Send it back
    return hitPair;
}

    
reco::ClusterHit3D DBScanAlg_DUNE35t::makeHitTriplet(const reco::ClusterHit3D& pairUW,
                                             const reco::ClusterHit2D* hitV) const
{
    // No matter what happens, we will return a HitPair object
    unsigned statusBits(0x8);              // This one is easy... 0x8 means triplet
    double   position[] = {0.,0.,0.};
    double   totalCharge(0.);
    double   avePeakTime(0.);
    double   deltaPeakTime(1000.);
    double   sigmaPeakTime(0.);
    double   overlapFraction(0.);
    double   hitDocaToAxis(9999.);
    double   hitArclenToPoca(0.);
        
    std::vector<const reco::ClusterHit2D*> hitVector;
    
    // Recover all the hits involved
    const reco::ClusterHit2D* hitU(pairUW.getHits()[0]);
    const reco::ClusterHit2D* hitW(pairUW.getHits()[1]);
    
    // Let's do a quick consistency check on the input hits to make sure we are in range...
    //double uvDeltaT   = std::max(pairUV.getDeltaPeakTime(), 2.);  // Set a lower limit
    //double uvSigma    = pairUV.getSigmaPeakTime();
    //double wSigma     = (hitW->getHit().RMS() * 2) / 2.355;
    double uwDeltaT   = 1.;
    double uwSigma    = 2.355 * pairUW.getSigmaPeakTime();
    double vSigma     = 2. * hitV->getHit().RMS();
    
    // Require the W hit to be "in range" with the UV Pair
    if (fabs(hitV->getTimeTicks() - pairUW.getAvePeakTime()) < uwDeltaT * (uwSigma + vSigma))
    {
        // Use the existing code to see the U and W hits are willing to pair with the V hit
        reco::ClusterHit3D pairUV = makeHitPair(hitU, hitV, 1.); //3.);
        reco::ClusterHit3D pairVW = makeHitPair(hitW, hitV, 1.); //3.);
    
        // If pairs made here then we can make a triplet
        // This condition means one of the two must be good
        if (pairUV.getAvePeakTime() > 0. || pairVW.getAvePeakTime() > 0.)
        {
            double numGoodPairs(1.);
            
            // Position is simply the average
            position[0]     = pairUW.getPosition()[0];
            position[1]     = pairUW.getPosition()[1];
            position[2]     = pairUW.getPosition()[2];
        
            // Ave, delta and sigmas
            avePeakTime     = pairUW.getAvePeakTime();
            deltaPeakTime   = pairUW.getDeltaPeakTime();
            sigmaPeakTime   = pairUW.getSigmaPeakTime() * pairUW.getSigmaPeakTime();
        
            // Overlap fraction... hmmm....
            overlapFraction = pairUW.getOverlapFraction();
        
            // Finally, total charge
            totalCharge     = pairUW.getTotalCharge();
            
            // Is the first W pair good?
            if (pairUV.getAvePeakTime() > 0.)
            {
                position[0]     += pairUV.getPosition()[0];
                position[1]     += pairUV.getPosition()[1];
                position[2]     += pairUV.getPosition()[2];
                avePeakTime     += pairUV.getAvePeakTime();
                deltaPeakTime    = std::max(deltaPeakTime, pairUV.getDeltaPeakTime());
                sigmaPeakTime   += pairUV.getSigmaPeakTime() * pairUV.getSigmaPeakTime();
                overlapFraction  = std::max(overlapFraction, pairUV.getOverlapFraction());
                totalCharge     += pairUV.getTotalCharge();
                numGoodPairs++;
            }
            
            // Is the second W pair good?
            if (pairVW.getAvePeakTime() > 0.)
            {
                position[0]     += pairVW.getPosition()[0];
                position[1]     += pairVW.getPosition()[1];
                position[2]     += pairVW.getPosition()[2];
                avePeakTime     += pairVW.getAvePeakTime();
                deltaPeakTime    = std::max(deltaPeakTime, pairVW.getDeltaPeakTime());
                sigmaPeakTime   += pairUV.getSigmaPeakTime() * pairVW.getSigmaPeakTime();
                overlapFraction  = std::max(overlapFraction, pairVW.getOverlapFraction());
                totalCharge     += pairVW.getTotalCharge();
                numGoodPairs++;
            }
            
            // Get averages
            position[0]   /= numGoodPairs;
            position[1]   /= numGoodPairs;
            position[2]   /= numGoodPairs;
            avePeakTime   /= numGoodPairs;
            sigmaPeakTime  = sqrt(sigmaPeakTime);

            // Remember our hits!
            hitVector.push_back(hitU);
            hitVector.push_back(hitV);
            hitVector.push_back(hitW);
        }
    }
        
    // Create the 3D cluster hit
    reco::ClusterHit3D hitTriplet(0,
                                  statusBits,
                                  position,
                                  totalCharge,
                                  avePeakTime,
                                  deltaPeakTime,
                                  sigmaPeakTime,
                                  overlapFraction,
                                  hitDocaToAxis,
                                  hitArclenToPoca,
                                  hitVector);
    // Send it back
    return hitTriplet;
}

    
const reco::ClusterHit2D* DBScanAlg_DUNE35t::FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, double pairDeltaTimeLimits) const
{
    static const double minCharge(0.);
    
    const reco::ClusterHit2D* bestVHit(0);
    
    double pairAvePeakTime(pair.getAvePeakTime());
    
    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit().Integral() < minCharge) continue;
        
        double hitVPeakTime(hit2D->getTimeTicks());
        double deltaPeakTime(pairAvePeakTime-hitVPeakTime);
        
        if (deltaPeakTime >  pairDeltaTimeLimits) continue;
        
	//        if (deltaPeakTime < -pairDeltaTimeLimits) break;
	//REL
        if (deltaPeakTime < -pairDeltaTimeLimits) continue;

        pairDeltaTimeLimits = fabs(deltaPeakTime);
        bestVHit            = hit2D;
    }

    return bestVHit;
}
    
int DBScanAlg_DUNE35t::FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, double range) const
{
    static const double minCharge(0.);
    
    int    numberInRange(0);
    double pairAvePeakTime(pair.getAvePeakTime());
    
    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit().Integral() < minCharge) continue;
        
        double hitVPeakTime(hit2D->getTimeTicks());
        double deltaPeakTime(pairAvePeakTime-hitVPeakTime);
        
        if (deltaPeakTime >  range) continue;
        
        if (deltaPeakTime < -range) break;
        
        numberInRange++;
    }
    
    return numberInRange;
}


geo::WireID DBScanAlg_DUNE35t::NearestWireID(const double* position, const geo::View_t& view) const
{
  geo::WireID wireID(0,view,0,0);
  
    // Embed the call to the geometry's services nearest wire id method in a try-catch block
    try {
        wireID =  m_geometry->NearestWireID(position, view);
    }
    catch(std::exception& exc)
    {
        // This can happen, almost always because the coordinates are **just** out of range
        mf::LogWarning("Cluster3D") << "Exception caught finding nearest wire, position - " << exc.what() << std::endl;
        
        // Assume extremum for wire number depending on z coordinate
        if (position[2] < 0.5 * m_geometry->DetLength()) wireID.Wire = 0;
        else                                             wireID.Wire = m_geometry->Nwires(view) - 1;
    }
    
    return wireID;
}


geo::WireID DBScanAlg_DUNE35t::NearestWireID_mod(const double* position, const geo::PlaneID & thePlaneID ) const
{
    geo::WireID wireID(0,0,0,0);
    
    // Embed the call to the geometry's services nearest wire id method in a try-catch block
    try {
      wireID =  m_geometry->NearestWireID(position, thePlaneID.Plane, thePlaneID.TPC, thePlaneID.Cryostat);
    }
    catch(std::exception& exc)
    {
        // This can happen, almost always because the coordinates are **just** out of range
        mf::LogWarning("Cluster3D") << "Exception caught finding nearest wire, position - " << exc.what() << std::endl;
	std::cout << "Exception caught finding nearest wire." << std::endl;


        // Assume extremum for wire number depending on z coordinate
        if (position[2] < 0.5 * m_geometry->DetLength()) wireID.Wire = 0;
        else                                             wireID.Wire = m_geometry->Nwires(thePlaneID.Plane) - 1; //Could be a problem here.Testing for now...

	std::cout << "WireID: " << wireID.Wire << std::endl;

    }
    
    return wireID;
}

}





  //size_t DBScanAlg_DUNE35t::BuildHitPairMap(ViewToHitVectorMap& viewToHitVectorMap, ViewToWireToHitSetMap& viewToWireToHitSetMap, HitPairList& hitPairList) const
  //{
    /**
     *  @brief Given input 2D hits, build out the lists of possible 3D hits
     *
     *         The current strategy: ideally all 3D hits would be comprised of a triplet of 2D hits, one from each view
     *         However, we have concern that, in particular, the v-plane may have some inefficiency which we have to be
     *         be prepared to deal with. The idea, then, is to first make the association of hits in the U and W planes
     *         and then look for the match in the V plane. In the event we don't find the match in the V plane then we 
     *         will evaluate the situation and in some instances keep the U-W pairs in order to keep efficiency high.
     */
  /*    
    // Should we set a minimum total charge for a hit?
    size_t totalNumHits(0);
    size_t hitPairCntr(0);
    double minCharge[] = {0., 0., 0.};
    
    //o*********************************************************************************
    // First job is to find all pairs of hits in different views which are "consistent"
    // Start loops over rows in HitVectorMap
    ViewToHitVectorMap::iterator mapItrU = viewToHitVectorMap.find(geo::kU);
    
    if (mapItrU != viewToHitVectorMap.end())
    {
        HitVector& hitVectorU = mapItrU->second;
        
        totalNumHits += hitVectorU.size();
        
        ViewToHitVectorMap::iterator mapItrW = viewToHitVectorMap.find(geo::kW);
        
        if (mapItrW != viewToHitVectorMap.end())
        {
            HitVector& hitVectorW = mapItrW->second;
            
            totalNumHits += hitVectorW.size();
            
            if (viewToHitVectorMap.find(geo::kV) != viewToHitVectorMap.end())
                totalNumHits += viewToHitVectorMap[geo::kV].size();
            
            // Take advantage that hits are sorted in "start time order"
            // Set the inner loop iterator before starting loop over outer hits
            HitVector::iterator hitVectorWStartItr = hitVectorW.begin();
            
            // Now we loop over the hits in these two layers
            for (HitVector::iterator hitItrU = hitVectorU.begin(); hitItrU != hitVectorU.end(); hitItrU++)
            {
                const reco::ClusterHit2D* hitPtrU = *hitItrU;
                
                if (hitPtrU->getHit().Integral() < minCharge[hitPtrU->getHit().View()]) continue;
                
                // This will be used in each loop so dereference the peak time here
                double hitUPeakTime = hitPtrU->getTimeTicks();
                
                // Inner loop is over hits within the time range of the outer hit
                for (HitVector::iterator hitItrW = hitVectorWStartItr; hitItrW != hitVectorW.end(); hitItrW++)
                {
                    const reco::ClusterHit2D* hitPtrW = *hitItrW;
                    
                    if (hitPtrW->getHit().Integral() < minCharge[hitPtrW->getHit().View()]) continue;
                    
                    // Hits are sorted in "peak time order" which we can take advantage of to try to speed
                    // the loops. Basically, we can compare the peak time for the outer loop hit against
                    // the current hit's peak time and if outside the range of interest we can take action.
                    // The range of interest is eyeballed but is meant to account for large pulse widths
                    // Start by dereferencing the inner hit peak time
                    double hitWPeakTime = hitPtrW->getTimeTicks();
                    
                    // If the outer loop's peak time is well past the current hit's then
                    // we should advance the inner loop's start iterator and keep going
                    if (hitUPeakTime >  hitWPeakTime + m_timeAdvanceGap)
                    {
                        hitVectorWStartItr++;
                        continue;
                    }
                    
                    // If the inner loop hit start time is past the end of the outer's end time, then break out loop
                    if (hitUPeakTime + m_timeAdvanceGap < hitWPeakTime) break;
                    
                    // We have a candidate hit pair combination, try to make a hit
                    reco::ClusterHit3D pair = makeHitPair(hitPtrU, hitPtrW);
                    
                    // The sign of success here is that the average peak time of the combined hits is > 0
                    // (note that when hits are combined the first window offset is accounted for)
                    if (pair.getAvePeakTime() > 0.)
                    {
                        bool hitInVNotFound(true);
                        
                        // Recover the WireID nearest in the V plane to the position of the pair
                        const geo::WireID wireIDV = NearestWireID(pair.getPosition(), geo::kV);
                            
                        // We believe the code that returns the ID is offset by one
                        WireToHitSetMap::iterator wireToHitSetMapVItr = viewToWireToHitSetMap[geo::kV].find(wireIDV.Wire);
                        
                        if (wireToHitSetMapVItr != viewToWireToHitSetMap[geo::kV].end())
                        {
                            const reco::ClusterHit2D* hit2DV = FindBestMatchingHit(wireToHitSetMapVItr->second, pair, m_numSigmaPeakTime*pair.getSigmaPeakTime());

                            // If a V hit found then it should be straightforward to make the triplet
                            if (hit2DV)
                            {
                                reco::ClusterHit3D triplet = makeHitTriplet(pair, hit2DV);
                            
                                if (triplet.getAvePeakTime() > 0.)
                                {
                                    triplet.setID(hitPairCntr++);
                                    hitPairList.emplace_back(std::unique_ptr<reco::ClusterHit3D>(new reco::ClusterHit3D(triplet)));
                                }
                                
                                hitInVNotFound = false;
                            }
                        }
                        
                        // NO V hit found, let's consider allowing for the pair to be kept
                        // WARNING: ugly code coming since I've not thought about how to do this more elegantly
                        if (hitInVNotFound)
                        {
                            // Idea will be to search either side of the wire in question looking to see how big a "gap"
                            // to any hits on adjacent wires that are "in time"
                            int lowSideGap(0);
                            int nLowSideHits(0);
                            int hiSideGap(0);
                            int nHiSideHits(0);
                            
                            // Idea is to check up to 3 wires either side of the target wire
                            // Do this in two pieces since my brain is not seeing the easy way at the moment
                            for(int loopIdx = 0; loopIdx < 3; loopIdx++)
                            {
                                // Set and check wire index
                                int wireIdx = wireIDV.Wire - loopIdx - 1;
                                
                                // Make sure we are in range
                                if (wireIdx < 0) break;
                                
                                // Find the set of hits along the current wire
                                WireToHitSetMap::iterator wireToHitSetMapVItr = viewToWireToHitSetMap[geo::kV].find(wireIdx);
                                
                                // Check there is something there and if there is find the closest hit
                                if (wireToHitSetMapVItr != viewToWireToHitSetMap[geo::kV].end())
                                {
                                    // Set a range that increases as we move further away in wires
                                    double range      = 3. * double(loopIdx + 1) * pair.getSigmaPeakTime();
                                    int    numInRange = FindNumberInRange(wireToHitSetMapVItr->second, pair, range);
                                    
                                    // Increment count of wires with hits near target wire
//                                    nLowSideHits++;
                                    
                                    // If hits were found in range then we are done with this loop
                                    if (numInRange > 0)
                                    {
                                        nLowSideHits++;
                                        break;
                                    }
                                    // Otherwise it is a gap and we do the accounting
                                    else lowSideGap++;
                                }
                            }
                            
                            // Second piece tries to check the high side
                            for(int loopIdx = 0; loopIdx < 3; loopIdx++)
                            {
                                int wireIdx = wireIDV.Wire + loopIdx + 1;
                                
                                if (wireIdx >= int(m_geometry->Nwires(geo::kV))) break;
                                
                                // Find the set of hits along the current wire
                                WireToHitSetMap::iterator wireToHitSetMapVItr = viewToWireToHitSetMap[geo::kV].find(wireIdx);
                                
                                // Check there is something there and if there is find the closest hit
                                if (wireToHitSetMapVItr != viewToWireToHitSetMap[geo::kV].end())
                                {
                                    double range      = 3. * double(loopIdx + 1) * pair.getSigmaPeakTime();
                                    int    numInRange = FindNumberInRange(wireToHitSetMapVItr->second, pair, range);
                                    
                                    // Again, increment counter
//                                    nHiSideHits++;
                                    
                                    // If a V hit found then check how close to wire
                                    if (numInRange)
                                    {
                                        nHiSideHits++;
                                        break;
                                    }
                                    else hiSideGap++;
                                }
                            }
                        
                            // Condition to keep the hit pair instead of requiring a triplet is that there is a partner hit
                            // on one side and a gap of no more than 2 wires on the other
                            if ((nLowSideHits > 0 && nHiSideHits > 0) && (lowSideGap == 0 || hiSideGap == 0))
                            {
                                pair.setID(hitPairCntr++);
                                hitPairList.emplace_back(std::unique_ptr<reco::ClusterHit3D>(new reco::ClusterHit3D(pair)));
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Return the hit pair list but sorted by z and y positions (faster traversal in next steps)
    hitPairList.sort(SetPositionOrder);
    
    // Where are we?
    mf::LogDebug("Cluster3D") << "Total number hits: " << totalNumHits << std::endl;
    mf::LogDebug("Cluster3D") << "Created a total of " << hitPairList.size() << " hit pairs, counted: " << hitPairCntr << std::endl;
    
    return hitPairList.size();


  */



//size_t DBScanAlg_DUNE35t::BuildNeighborhoodMap(HitPairList& hitPairList, EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec) const
//{
    /**
     *  @brief build out the epsilon neighborhood map to be used by DBScan
     */
  /*    
    size_t consistentPairsCnt(0);
    size_t pairsChecked(0);
    
    //o**********************************************************************************
    // Given the list of pairs of hits which are consistent with each other, build out the
    // epsilon neighbor maps
    // Now we go through the pair list and basically repeat the same exercise as above
    // The following assumes that the HitPairList is ordered
    // a) in increasing Z for hits which are not on the "same W wire",
    // b) in increasing U (Y) for hits on the same W wire
    
    for (HitPairList::const_iterator pairItrO = hitPairList.begin(); pairItrO != hitPairList.end(); pairItrO++)
    {
        const reco::ClusterHit3D* hitPairO   = (*pairItrO).get();
        const size_t              hitPairOID = hitPairO->getID();
        
        // Need to initialize the "first" part of the vector pseudo map
        epsPairNeighborhoodMapVec[hitPairOID].first = hitPairO;
        
        // Get reference to the list for this hit so we don't look it up inside the loop
        EpsPairNeighborhoodPair& hitPairOPair(epsPairNeighborhoodMapVec[hitPairOID].second);
        
        HitPairList::const_iterator pairItrI = pairItrO;
        
        std::map<int, std::pair<double, const reco::ClusterHit3D*> > bestTripletMap;
        
        // Get the wire numbers for this triplet
        int pairOWireU(hitPairO->getHits().front()->getHit().WireID().Wire);
        int pairOWireV(0);
        int pairOWireW(hitPairO->getHits().back()->getHit().WireID().Wire);
        
        if (hitPairO->getHits().size() > 2)
            pairOWireV = hitPairO->getHits()[1]->getHit().WireID().Wire;
        else
        {
            const geo::WireID wireIDV = NearestWireID(hitPairO->getPosition(), geo::kV);
            
            pairOWireV = wireIDV.Wire;
        }
        
        // Set maximums
        int maxDeltaW(2);
        int maxDeltaU(2);
        int maxDeltaV(2);
        int maxSumAbsUV(4);
        
        while (++pairItrI != hitPairList.end())
        {
            const reco::ClusterHit3D* hitPairI = (*pairItrI).get();
            
            // Note that the hits have been sorted by z, then y and then in x
            // Translation: hits are sorted by W wire, then in Y (increasing u, decreasing v)
            // and then in time.
            int pairIWireU(hitPairI->getHits().front()->getHit().WireID().Wire);
            int pairIWireV(0);
            int pairIWireW(hitPairI->getHits().back()->getHit().WireID().Wire);
            
            if (hitPairI->getHits().size() > 2)
                pairIWireV = hitPairI->getHits()[1]->getHit().WireID().Wire;
            else
            {
                const geo::WireID wireIDV = NearestWireID(hitPairI->getPosition(), geo::kV);
                
                pairIWireV = wireIDV.Wire;
            }
            
            // Form the wire number differences so we can check ranges
            int deltaU(pairIWireU - pairOWireU);
            int deltaV(pairIWireV - pairOWireV);
            int deltaW(pairIWireW - pairOWireW);
            
            // If we have passed the max delta W then we are done with the loop
            if (deltaW >  maxDeltaW) break;
            
            // Check limits on U,V differences and if past then continue to next
            if (fabs(deltaU) > maxDeltaU || fabs(deltaV) > maxDeltaV) continue;
            
            // Special case
            if (fabs(deltaU) + fabs(deltaV) > maxSumAbsUV) continue;
            
            // Keep count...
            pairsChecked++;
            
            // This is the tight constraint on the hits
            if (consistentPairs(hitPairO, hitPairI))
            {
                int    bestBin  = 100 * deltaW + 10 * deltaU + deltaV;
                double bestDist = 10000.;
                
                if (bestTripletMap.find(bestBin) != bestTripletMap.end())
                {
                    bestDist = bestTripletMap[bestBin].first;
                }
                
                double newDist = fabs(hitPairI->getX() - hitPairO->getX());
                
                // This is an attempt to "prefer" triplets over pairs
                if (hitPairI->getHits().size() < 3) newDist += 25.;
                
                if (newDist < bestDist) bestTripletMap[bestBin] = std::pair<double, const reco::ClusterHit3D*>(newDist, hitPairI);
                
                // Check limits
                if      (deltaW == 0 && deltaU == -1 && deltaV == 1) maxSumAbsUV = 2;
                else if (deltaW == 1 && ((deltaU == 0 && deltaV == 1) || (deltaU == 1 && deltaV == 0))) maxDeltaW = 1;
            }
        }
        
        for(const auto& bestMapItr : bestTripletMap)
        {
            const reco::ClusterHit3D* hitPairI(bestMapItr.second.second);
            
            hitPairOPair.first.incrementCount();
            hitPairOPair.second.emplace_back(hitPairI);
            
            epsPairNeighborhoodMapVec[hitPairI->getID()].second.first.incrementCount();
            epsPairNeighborhoodMapVec[hitPairI->getID()].second.second.emplace_back(hitPairO);
            
            consistentPairsCnt++;
        }
    }
    
    mf::LogDebug("Cluster3D") << "Consistent pairs: " << consistentPairsCnt << " of " << pairsChecked << " checked." << std::endl;
    
    return consistentPairsCnt;
}
*/



