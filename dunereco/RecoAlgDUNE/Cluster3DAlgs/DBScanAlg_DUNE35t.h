/**
 *  @file   DBScanAlg_DUNE35t.h
 * 
 *  @brief  This algorithm will create and then cluster 3D hits using DBScan
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef DBScanAlg_DUNE35t_h
#define DBScanAlg_DUNE35t_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"

// std includes
#include <vector>
#include <list>
#include <set>
#include <map>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief a utility class for keeping track of the state of a hit for DBScan
 */
class DBScanParams
{
public:
    DBScanParams() : m_visited(false), m_noise(false), m_inCluster(false), m_count(0) {}
        
    void setVisited()         {m_visited   = true;}
    void setNoise()           {m_noise     = true;}
    void setInCluster()       {m_inCluster = true;}
    void setCount(int count)  {m_count     = count;}
    
    void clearVisited()                   const {m_visited   = false;}
    void incrementCount(size_t count = 1) const {m_count    += count;}
        
    bool   visited()   const {return m_visited;}
    bool   isNoise()   const {return m_noise;}
    bool   inCluster() const {return m_inCluster;}
    size_t getCount()  const {return m_count;}
        
private:
    mutable bool   m_visited;
    bool           m_noise;
    bool           m_inCluster;
    mutable size_t m_count;
};
    
/**
 *   @brief What follows are several highly useful typedefs which we 
 *          want to expose to the outside world
 */

typedef std::vector<reco::ClusterHit2D*>       HitVector;
typedef std::vector<const reco::ClusterHit2D*> HitVectorConst;
typedef std::map<int, HitVector >              HitClusterMap;
typedef std::map<geo::View_t, HitVector>       ViewToHitVectorMap;

// forward declaration to define an ordering function for our hit set
struct Hit2DSetCompare
{
    bool operator() (const reco::ClusterHit2D*, const reco::ClusterHit2D*) const;
};
    
typedef std::vector<reco::ClusterHit2D>                      Hit2DVector;
typedef std::set<const reco::ClusterHit2D*, Hit2DSetCompare> Hit2DSet;
typedef std::map<unsigned int, Hit2DSet >                    WireToHitSetMap;
typedef std::map<geo::View_t, WireToHitSetMap >              ViewToWireToHitSetMap;
typedef std::map<geo::View_t, HitVector >                    HitVectorMap;
    
typedef std::vector<std::unique_ptr<reco::ClusterHit3D> >    HitPairVector;
typedef std::list<std::unique_ptr<reco::ClusterHit3D> >      HitPairList;

/**
 *  @brief  DBScanAlg_DUNE35t class definiton
 */
class DBScanAlg_DUNE35t
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    DBScanAlg_DUNE35t(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~DBScanAlg_DUNE35t();

    void reconfigure(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     */
    void ClusterHitsDBScan(ViewToHitVectorMap&      viewToHitVectorMap,
                           ViewToWireToHitSetMap&   viewToWiretoHitSetMap,
                           HitPairList&             hitPairList,
                           reco::HitPairClusterMap& hitPairClusterMap);

    /**
     *  @brief enumerate the possible values for time checking if monitoring timing
     */
    enum TimeValues {BUILDTHREEDHITS  = 0,
                     BUILDHITTOHITMAP = 1,
                     RUNDBSCAN        = 2,
                     NUMTIMEVALUES
    };
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    double getTimeToExecute(TimeValues index) const {return m_timeVector.at(index);}
    
private:
    
    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    size_t BuildHitPairMap(ViewToHitVectorMap& viewToHitVectorMap, ViewToWireToHitSetMap& viewToWiretoHitSetMap, HitPairList& hitPairList) const;
    
    /**
     *  @brief Make a HitPair object by checking two hits
     */
    reco::ClusterHit3D makeHitPair(const reco::ClusterHit2D* hit1,
                                   const reco::ClusterHit2D* hit2,
                                   double                    hitWidthSclFctr = 1.,
                                   size_t                    hitPairCntr = 0) const;
    
    /**
     *  @brief Make a HitPair object by checking two hits
     */
    reco::ClusterHit3D makeHitTriplet(const reco::ClusterHit3D& pair,
                                      const reco::ClusterHit2D* hit2) const;
    
    /**
     *  @brief A utility routine for finding a 2D hit closest in time to the given pair
     */
    const reco::ClusterHit2D* FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, double pairDeltaTimeLimits) const;
    
    /**
     *  @brief A utility routine for returning the number of 2D hits from the list in a given range
     */
    int FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, double range) const;
    
    /**
     *  @brief The bigger question: are two pairs of hits consistent?
     */
    bool consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2) const;
    
    typedef std::list<const reco::ClusterHit3D*>                           EpsPairNeighborhoodList;
    typedef std::pair<DBScanParams, EpsPairNeighborhoodList >              EpsPairNeighborhoodPair;
    typedef std::map<const reco::ClusterHit3D*, EpsPairNeighborhoodPair >  EpsPairNeighborhoodMap;
    typedef std::pair<const reco::ClusterHit3D*, EpsPairNeighborhoodPair > EpsPairNeighborhoodMapPair;
    typedef std::vector<EpsPairNeighborhoodMapPair >                       EpsPairNeighborhoodMapVec;
    
    /**
     *  @brief the main routine for DBScan
     */
    void expandCluster(EpsPairNeighborhoodMapVec&          nMap,
                       EpsPairNeighborhoodMapVec::iterator vecItr,
                       reco::HitPairListPtr&               cluster,
                       size_t                              minPts) const;
    
    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    size_t BuildNeighborhoodMap(HitPairList& hitPairList,
				EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec )const;

    
    /** 
     *  @brief Jacket the calls to finding the nearest wire in order to intercept the exceptions if out of range
     */
    geo::WireID NearestWireID(const double* position, const geo::View_t& view) const;
    geo::WireID NearestWireID_mod(const double* position, const geo::PlaneID & thePlaneID ) const;

    
    /**
     *  @brief Data members to follow
     */
    
    size_t                    m_minPairPts;
    double                    m_timeAdvanceGap;
    double                    m_numSigmaPeakTime;
    double                    m_EpsMaxDist;

    bool                      m_enableMonitoring;      ///<
    int                       m_hits;                  ///<
    std::vector<float>        m_timeVector;            ///<

    geo::Geometry*            m_geometry;  // pointer to the Geometry service


    //    const detinfo::DetectorProperties* m_detector;  // Pointer to the detector properties
};

} // namespace lar_cluster3d
#endif
