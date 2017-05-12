/**
 *  @file   PCASeedFinderAlg.h
 * 
 *  @brief  This is an algorithm for finding recob::Seed objects in 3D clusters
 * 
 */
#ifndef PCASeedFinderAlg_h
#define PCASeedFinderAlg_h

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

/**
 *  @brief  PCASeedFinderAlg class
 */
class PCASeedFinderAlg : virtual public SeedFinderAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    PCASeedFinderAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~PCASeedFinderAlg();
    
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
    
    /**
     *  @brief Separate function to find hits at the ends of the input hits
     */
    bool getHitsAtEnd(reco::HitPairListPtr& hit3DList, reco::PrincipalComponents& seedPca) const;
    
    void LineFit2DHits(const reco::HitPairListPtr& hitList, double XOrigin, TVector3& Pos, TVector3& Dir, double& ChiDOF) const;

    geo::Geometry*                         m_geometry;         // pointer to the Geometry service
    //    const detinfo::DetectorProperties*    m_detector;         // Pointer to the detector properties

    double                                 m_gapDistance;      ///<
    size_t                                 m_numSeed2DHits;    ///<
    double                                 m_minAllowedCosAng; ///< The minimum cos(ang) between input and seed axes
    
    PrincipalComponentsAlg                 m_pcaAlg;           // For running Principal Components Analysis
};

} // namespace lar_cluster3d
#endif
