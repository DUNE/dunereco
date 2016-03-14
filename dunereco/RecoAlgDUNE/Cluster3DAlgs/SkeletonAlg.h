/**
 *  @file   SkeletonAlg.h
 * 
 *  @brief  Header file to define the interface to the SkeletonAlg
 * 
 */
#ifndef SkeletonAlg_h
#define SkeletonAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoObjects/Cluster3D.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{

/**
 *  @brief  Cluster3D class
 */
class SkeletonAlg
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    SkeletonAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~SkeletonAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    void reconfigure(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief This is intended to find the medial skeleton given a list of input hit pairs
     *
     *  @param hitPairList - input list of pointers to internal Cluster3D 3D hits
     */
    int FindMedialSkeleton(reco::HitPairListPtr& hitPairList) const;
    
    /**
     *  @brief Return the skeleton hits from the input list
     *         - note that this presumes the skeleton hits have been found already
     *
     *  @param inputHitList - input list of pointers to internal Cluster3D 3D hits
     *  @param skeletonHitList - output list of skeleton hits
     */
    void GetSkeletonHits(const reco::HitPairListPtr& inputHitList, reco::HitPairListPtr& skeletonHitList) const;
    
    /**
     *  @brief Modifies the position of input skeleton hits by averaging along the "best" 
     *         wire direction.
     *
     *  @param skeletonHitList - input list of skeleton hits
     */
    void AverageSkeletonPositions(reco::HitPairListPtr& skeletonHitList) const;
    
private:
    
    /**
     *  @brief A function to find the bounding wires in a given view 
     *
     */
    double FindFirstAndLastWires(std::vector<const reco::ClusterHit3D*>& hitVec,
                                 int                                     viewToCheck,
                                 int                                     referenceWire,
                                 double                                  referenceTicks,
                                 int&                                    firstWire,
                                 int&                                    lastWire) const;
    
    double                    m_minimumDeltaTicks;
    double                    m_maximumDeltaTicks;

    fhicl::ParameterSet       m_pset;
};

} // namespace lar_cluster3d
#endif
