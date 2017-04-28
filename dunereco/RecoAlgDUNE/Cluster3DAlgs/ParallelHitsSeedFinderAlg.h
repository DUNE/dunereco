/**
 *  @file   ParallelHitsSeedFinderAlg.h
 * 
 *  @brief  This is an algorithm for finding recob::Seed objects in 3D clusters
 * 
 */
#ifndef ParallelHitsSeedFinderAlg_h
#define ParallelHitsSeedFinderAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

#include "dune/RecoAlgDUNE/Cluster3DAlgs/SeedFinderAlgBase.h"
#include "dune/RecoAlgDUNE/Cluster3DAlgs/PrincipalComponentsAlg.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardata/RecoObjects/Cluster3D.h"

// ROOT includes
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
    
//typedef std::pair<recob::Seed, reco::HitPairListPtr> SeedHitListPair;
//typedef std::list<SeedHitListPair >                  SeedHitPairList;


/**
 *  @brief  ParallelHitsSeedFinderAlg class
 */
class ParallelHitsSeedFinderAlg : virtual public SeedFinderAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    ParallelHitsSeedFinderAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~ParallelHitsSeedFinderAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const &pset);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    virtual bool findTrackSeeds(reco::HitPairListPtr&      hitPairListPtr,
                                reco::PrincipalComponents& inputPCA,
                                SeedHitPairListPairVec&    seedHitMap) const;

private:

    geo::Geometry*                         m_geometry;        // pointer to the Geometry service
    //    const detinfo::DetectorProperties*    m_detector;        // Pointer to the detector properties

    size_t                                 m_maxNumEdgeHits;  ///< Maximum number hits each end of PCA axis
    double                                 m_gapDistance;     ///< Maximum allowed distance between hits
    size_t                                 m_numSeed2DHits;   ///< Number 2D seed hits desired
    
    PrincipalComponentsAlg                 m_pcaAlg;          // For running Principal Components Analysis
};

} // namespace lar_cluster3d
#endif
