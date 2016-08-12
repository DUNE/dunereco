/**
 *  @file   SeedFinderAlgBase.h
 * 
 *  @brief  This is intended to define an interface to all Seed finder algorithms employed
 *          by the 3D clustering
 * 
 */
#ifndef SeedFinderAlgBase_h
#define SeedFinderAlgBase_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardata/RecoObjects/Cluster3D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
typedef std::pair<recob::Seed, reco::HitPairListPtr> SeedHitPairListPair;
typedef std::vector<SeedHitPairListPair >            SeedHitPairListPairVec;

/**
 *  @brief  SeedFinderAlgBase class
 */
class SeedFinderAlgBase
{
public:
    /**
     *  @brief Require that a handler is definied in case the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const &pset) = 0;

    /**
     *  @brief Define the interface to take an input list of 3D hits and return seed candidates
     *         so hits are ordered along the axis
     */
    virtual bool findTrackSeeds(reco::HitPairListPtr&       hitPairListPtr,
                                reco::PrincipalComponents&  inputPCA,
                                SeedHitPairListPairVec&     seedHitPairVec) const = 0;

protected:

    /**
     *  @brief Define a comparator which will sort hits by arc length along a PCA axis
     */
    struct Sort3DHitsByArcLen3D
    {
        bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
        {
            return left->getArclenToPoca() < right->getArclenToPoca();
        }
        
    };
    
    /**
     *  @brief Define a comparator which will sort hits by the absolute value of arc length
     *         so hits are ordered closed to PCA origin to furthest
     */
    struct Sort3DHitsByAbsArcLen3D
    {
        bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
        {
            return fabs(left->getArclenToPoca()) < fabs(right->getArclenToPoca());
        }
        
    };
    
private:
};

} // namespace lar_cluster3d
#endif
