/**
 *  @file   HoughSeedFinderAlg.cxx
 * 
 *  @brief  Implementation of the Seed Finder Algorithm based on a Hough Transform
 *
 *          The intent of this algorithm is to take an input list of 3D space points and from those
 *          to find candidate track start points and directions. It does so by first performing a 
 *          Principal Components Analysis of the input 3D hits and then projects them to the plane
 *          of largest spread. A standard Hough Transform method is then applied to attempt to identify
 *          straight line segments which can be used as seeds to the kalman filter tracker.
 */

// The main include
#include "dune/RecoAlgDUNE/Cluster3DAlgs/HoughSeedFinderAlg.h"
// Framework Includes

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// ROOT includes
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2D.h"
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

HoughSeedFinderAlg::HoughSeedFinderAlg(fhicl::ParameterSet const &pset) :
     m_minimum3DHits(5),
     m_thetaBins(360),
     m_rhoBins(75),
     m_hiThresholdMin(5),
     m_hiThresholdFrac(.05),
     m_loThresholdFrac(0.85),
     m_numSeed2DHits(80),
     m_numAveDocas(6.),
     m_numSkippedHits(10),
     m_maxLoopsPerCluster(3),
     m_maximumGap(5.),
     m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg")),
     m_displayHist(false)
{
    this->reconfigure(pset);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    //    auto const* detectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    m_geometry = &*geometry;
    //    m_detector = detectorProperties->provider();
}

//------------------------------------------------------------------------------------------------------------------------------------------

HoughSeedFinderAlg::~HoughSeedFinderAlg()
{
}
    
void HoughSeedFinderAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_minimum3DHits      = pset.get<size_t>("Minimum3DHits",           5);
    m_thetaBins          = pset.get<int>   ("ThetaBins",             360);
    m_rhoBins            = pset.get<int>   ("HalfRhoBins",            75);
    m_hiThresholdMin     = pset.get<size_t>("HiThresholdMin",          5);
    m_hiThresholdFrac    = pset.get<double>("HiThresholdFrac",      0.05);
    m_loThresholdFrac    = pset.get<double>("LoThresholdFrac",      0.85);
    m_numSeed2DHits      = pset.get<size_t>("NumSeed2DHits",          80);
    m_numAveDocas        = pset.get<double>("NumAveDocas",            6.);
    m_numSkippedHits     = pset.get<int>   ("NumSkippedHits",         10);
    m_maxLoopsPerCluster = pset.get<int>("MaxLoopsPerCluster",         3);
    m_maximumGap         = pset.get<double>("MaximumGap",             5.);
    m_displayHist        = pset.get<bool>  ("DisplayHoughHist",    false);
    
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
}
    
class AccumulatorValues
{
    /**
     *  @brief A utility class to contain the values of a given "bin" in Hough Space
     *
     *         Specifically, this is keeping track of the projected x,y coordinates of a given
     *         3D hit projected to the plane of largest spread in PCA and an interator to 
     *         that hit in the input container
     */
public:
    AccumulatorValues() : m_position(0.,0.,0.) {}
    AccumulatorValues(const TVector3& position, const reco::HitPairListPtr::const_iterator& itr) :
         m_position(position), m_hit3DIterator(itr) {}
     
    const TVector3&                      getPosition()    const {return m_position;}
    reco::HitPairListPtr::const_iterator getHitIterator() const {return m_hit3DIterator;}
     
private:
    TVector3                             m_position;       ///< We really only need the x,y coordinates here but keep all three for now
    reco::HitPairListPtr::const_iterator m_hit3DIterator;  ///< This will be used to take us back to our 3D hit
};
     
typedef std::vector<AccumulatorValues> AccumulatorValuesVec;
     
class HoughSeedFinderAlg::AccumulatorBin
{
    /**
     *  @brief A utility class used to accumulate the above values 
     *
     *         One of these objects will exist for each "bin" in rho-theta space and this will
     *         be used to accumulate the 3D hits which contribute to this bin
     */
public:
    AccumulatorBin() : m_visited(false), m_noise(false), m_inCluster(false) {}
    AccumulatorBin(AccumulatorValues& values) : m_visited(false), m_noise(false), m_inCluster(false)
    {
        m_accumulatorValuesVec.push_back(values);
    }
     
    void setVisited()   {m_visited   = true;}
    void setNoise()     {m_noise     = true;}
    void setInCluster() {m_inCluster = true;}
     
    void addAccumulatorValue(AccumulatorValues& value) {m_accumulatorValuesVec.push_back(value);}
     
    bool isVisited()   const {return m_visited;}
    bool isNoise()     const {return m_noise;}
    bool isInCluster() const {return m_inCluster;}
     
    const AccumulatorValuesVec& getAccumulatorValues() const {return m_accumulatorValuesVec;}
private:
    bool                 m_visited;
    bool                 m_noise;
    bool                 m_inCluster;
    AccumulatorValuesVec m_accumulatorValuesVec;
};

class HoughSeedFinderAlg::SortHoughClusterList
{
    /**
     * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
     */
public:
    SortHoughClusterList(HoughSeedFinderAlg::RhoThetaAccumulatorBinMap& accMap) : m_accMap(accMap) {}
     
    bool operator()(const HoughSeedFinderAlg::HoughCluster& left, const HoughSeedFinderAlg::HoughCluster& right)
    {
        size_t peakCountLeft(0);
        size_t peakCountRight(0);
     
        for(const auto& binIndex : left)
            peakCountLeft = std::max(peakCountLeft, m_accMap[binIndex].getAccumulatorValues().size());
        for(const auto& binIndex : right)
            peakCountRight = std::max(peakCountRight, m_accMap[binIndex].getAccumulatorValues().size());
     
        return peakCountLeft > peakCountRight;
    }
private:
    HoughSeedFinderAlg::RhoThetaAccumulatorBinMap& m_accMap;
};
    
struct HoughSeedFinderAlg::SortBinIndexList
{
    /**
     *  @brief This is used to sort bins in Hough Space
     */
    bool operator()(const HoughSeedFinderAlg::RhoThetaAccumulatorBinMap::iterator& left, const HoughSeedFinderAlg::RhoThetaAccumulatorBinMap::iterator& right)
    {
        size_t leftSize  = left->second.getAccumulatorValues().size();
        size_t rightSize = right->second.getAccumulatorValues().size();
     
        return  leftSize > rightSize;
    }
};
     
    
void HoughSeedFinderAlg::HoughRegionQuery(BinIndex&                  curBin,
                                          RhoThetaAccumulatorBinMap& rhoThetaAccumulatorBinMap,
                                          HoughCluster&              neighborPts,
                                          size_t                     threshold) const
{
    /**
     *   @brief Does a query of nearest neighbors to look for matching bins
     */
    
    // We simply loop over the nearest indices and see if we have any friends over threshold
    for(int rhoIdx = curBin.first - 1; rhoIdx <= curBin.first + 1; rhoIdx++)
    {
        for(int jdx = curBin.second - 1; jdx <= curBin.second + 1; jdx++)
        {
            // Skip the self reference
            if (rhoIdx == curBin.first && jdx == curBin.second) continue;
            
            // Theta bin needs to handle the wrap.
            int thetaIdx(jdx);
            
            if      (thetaIdx < 0)             thetaIdx = m_thetaBins - 1;
            else if (thetaIdx > m_thetaBins -1) thetaIdx = 0;
            
            BinIndex                            binIndex(rhoIdx,thetaIdx);
            RhoThetaAccumulatorBinMap::iterator mapItr = rhoThetaAccumulatorBinMap.find(binIndex);
            
            if (mapItr != rhoThetaAccumulatorBinMap.end())
            {
                if (mapItr->second.getAccumulatorValues().size() >= threshold) neighborPts.push_back(binIndex);
            }
        }
    }
    
    return;
}
    
void HoughSeedFinderAlg::expandHoughCluster(BinIndex&                  curBin,
                                            HoughCluster&              neighborPts,
                                            HoughCluster&              houghCluster,
                                            RhoThetaAccumulatorBinMap& rhoThetaAccumulatorBinMap,
                                            size_t                     threshold) const
{
    /** 
     *  @brief The workhorse routine for a DBScan like clustering routine to identify peak bins in Hough Space
     */
    
    // Start by adding the input point to our Hough Cluster
    houghCluster.push_back(curBin);
    
    for(auto& binIndex : neighborPts)
    {
        AccumulatorBin& accBin = rhoThetaAccumulatorBinMap[binIndex];
        
        if (!accBin.isVisited())
        {
            accBin.setVisited();
            
            HoughCluster nextNeighborPts;
            
            HoughRegionQuery(binIndex, rhoThetaAccumulatorBinMap, nextNeighborPts, threshold);
            
            if(!nextNeighborPts.empty())
            {
                for(auto& neighborIdx : nextNeighborPts)
                {
                    neighborPts.push_back(neighborIdx);
                }
            }
        }
        
        if (!accBin.isInCluster())
        {
            houghCluster.push_back(binIndex);
            accBin.setInCluster();
        }
    }
    
    return;
}
    
bool Hit3DCompare(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
{
    return *left < *right;
}
    
struct Hit3DSetCompare
{
    bool operator() (const reco::ClusterHit3D* left, const reco::ClusterHit3D* right) const {return Hit3DCompare(left, right);}
};
    
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
    
void HoughSeedFinderAlg::findHitGaps(reco::HitPairListPtr&      inputHitList,
                                     reco::HitPairListPtr&      outputList) const
{
    // Intention is to try to find the largest "contiguous" block of hits in the input list
    // In a nutshell, the idea is to order input hits according to the pca, then
    // loop through the hits and store them in a new hit list. If a gap is detected then
    // we terminate the previous list, create a new one and continue. After the loop over
    // hits is complete then simply pick the biggest list.
    // Step I is to run the pca and order the hits
    reco::PrincipalComponents pca;
    
    m_pcaAlg.PCAAnalysis_3D(inputHitList, pca, true);
    
    // It would seem that the analysis can fail!
    if (!pca.getSvdOK())
    {
        outputList = inputHitList;
        return;
    }
    
    m_pcaAlg.PCAAnalysis_calc3DDocas(inputHitList, pca);
    
    inputHitList.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
    outputList.clear();
    
    // Create containers to hold our hit lists
    reco::HitPairListPtr continuousHitList;
    
    std::map<size_t, reco::HitPairListPtr > gapHitListMap;
    
    // Remember the distance arc length of the last hit
    double arcLenLastHit = inputHitList.front()->getArclenToPoca();
    
    // Loop through the input hits
    for(const auto& hit3D : inputHitList)
    {
        // Hits in order, delta arclength should always be positive
        double arcLen      = hit3D->getArclenToPoca();
        double deltaArcLen = arcLen - arcLenLastHit;
        
        if (deltaArcLen > m_maximumGap)
        {
            gapHitListMap[continuousHitList.size()] = continuousHitList;
            continuousHitList.clear();
        }
        
        continuousHitList.emplace_back(hit3D);
        
        arcLenLastHit = arcLen;
    }
    
    if (!continuousHitList.empty()) gapHitListMap[continuousHitList.size()] = continuousHitList;
    
    // Get the largest list of hits
    std::map<size_t, reco::HitPairListPtr >::reverse_iterator longestListItr = gapHitListMap.rbegin();
    
    if (longestListItr != gapHitListMap.rend())
    {
        size_t               nContinuousHits = longestListItr->first;
        reco::HitPairListPtr longestList     = longestListItr->second;
    
        outputList.resize(nContinuousHits);
        std::copy(longestList.begin(), longestList.end(), outputList.begin());
    }

    return;
}
    
void HoughSeedFinderAlg::findHoughClusters(const reco::HitPairListPtr& hitPairListPtr,
                                           reco::PrincipalComponents&  pca,
                                           int&                        nLoops,
                                           RhoThetaAccumulatorBinMap&  rhoThetaAccumulatorBinMap,
                                           HoughClusterList&           houghClusters) const
{
    // The goal of this function is to do a basic Hough Transform on the input list of 3D hits.
    // In order to transform this to a 2D problem, the 3D hits are projected to the plane of the two
    // largest eigen values from the input principal components analysis axis.
    // There are two basic steps to the job here:
    // 1) Build and accumlate a rho-theta map
    // 2) Go through that rho-theta map and find candidate Hough "clusters"
    // Unfortunately, the following may not be suitable viewing for those who may be feint of heart
    //
    // Define some constants
    static int   histCount(0);
    const double maxTheta(M_PI);                               // Obviously, 180 degrees
    const double thetaBinSize(maxTheta/double(m_thetaBins));    // around 4 degrees (45 bins)
    const double rhoBinSizeMin(m_geometry->WirePitch());       // Wire spacing gives a natural bin size?
    
    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    TVector3 pcaCenter(pca.getAvePosition()[0],pca.getAvePosition()[1],pca.getAvePosition()[2]);
    TVector3 planeVec0(pca.getEigenVectors()[0][0],pca.getEigenVectors()[0][1],pca.getEigenVectors()[0][2]);
    TVector3 planeVec1(pca.getEigenVectors()[1][0],pca.getEigenVectors()[1][1],pca.getEigenVectors()[1][2]);
    TVector3 pcaPlaneNrml(pca.getEigenVectors()[2][0],pca.getEigenVectors()[2][1],pca.getEigenVectors()[2][2]);
    double   eigenVal0  = 3. * sqrt(pca.getEigenValues()[0]);
    double   eigenVal1  = 3. * sqrt(pca.getEigenValues()[1]);
    double   maxRho     = std::sqrt(eigenVal0*eigenVal0 + eigenVal1*eigenVal1) * 2. / 3.;
    double   rhoBinSize = maxRho / double(m_rhoBins);
    
    if (rhoBinSize < rhoBinSizeMin) rhoBinSize = rhoBinSizeMin;
    
    // **********************************************************************
    // Part I: Accumulate values in the rho-theta map
    // **********************************************************************
    
    size_t maxBinCount(0);
    int    nAccepted3DHits(0);
    
    // Commence looping over the input 3D hits and fill our accumulator bins
    for(reco::HitPairListPtr::const_iterator hit3DItr  = hitPairListPtr.begin();
                                             hit3DItr != hitPairListPtr.end();
                                             hit3DItr++)
    {
        // Recover the hit
        const reco::ClusterHit3D* hit3D(*hit3DItr);
        
        // Skip hits which are not skeleton points
        if (!(hit3D->getStatusBits() & 0x10000000)) continue;
        
        nAccepted3DHits++;
        
        TVector3 hit3DPosition(hit3D->getPosition()[0], hit3D->getPosition()[1], hit3D->getPosition()[2]);
        TVector3 pcaToHitVec = hit3DPosition - pcaCenter;
        TVector3 pcaToHitPlaneVec(pcaToHitVec.Dot(planeVec0), pcaToHitVec.Dot(planeVec1), 0.);
        
        double   xPcaToHit = pcaToHitPlaneVec[0];
        double   yPcaToHit = pcaToHitPlaneVec[1];
        
        // Create an accumulator value
        AccumulatorValues accValue(pcaToHitPlaneVec, hit3DItr);
        
        // Commence loop over theta to fill accumulator bins
        // Note that with theta in the range 0-pi then we can have negative values for rho
        for(int thetaIdx  = 0; thetaIdx < m_thetaBins; thetaIdx++)
        {
            // We need to convert our theta index to an angle
            double theta  = thetaBinSize * double(thetaIdx);
            
            // calculate rho for this angle
            double rho    = xPcaToHit * std::cos(theta) + yPcaToHit * std::sin(theta);
            
            // Get the rho index
            int    rhoIdx = std::round(rho / rhoBinSize);
            
            // Accumulate
            BinIndex binIndex(rhoIdx, thetaIdx);
            
            rhoThetaAccumulatorBinMap[binIndex].addAccumulatorValue(accValue);
            
            if (rhoThetaAccumulatorBinMap[binIndex].getAccumulatorValues().size() > maxBinCount)
                maxBinCount = rhoThetaAccumulatorBinMap[binIndex].getAccumulatorValues().size();
        }
    }
    
    // Accumulation done, if asked now display the hist
    if (m_displayHist)
    {
        std::ostringstream ostr;
        ostr << "Hough Histogram " << histCount++;
        m_Canvases.emplace_back(new TCanvas(ostr.str().c_str(), ostr.str().c_str(), 1000, 1000));
        
        std::ostringstream ostr2;
        ostr2 << "Plane";
        
        m_Canvases.back()->GetFrame()->SetFillColor(46);
        m_Canvases.back()->SetFillColor(19);
        m_Canvases.back()->SetBorderMode(19);
        m_Canvases.back()->cd(1);
        
        double zmin = 0.06;
        double zmax = 0.94;
        double xmin = 0.04;
        double xmax = 0.95;
        TPad* p = new TPad(ostr2.str().c_str(), ostr2.str().c_str(), zmin, xmin, zmax, xmax);
        p->SetBit(kCanDelete);   // Give away ownership.
        p->Range(zmin, xmin, zmax, xmax);
        p->SetFillStyle(4000);   // Transparent.
        p->Draw();
        m_Pads.push_back(p);
        
        TH2D* houghHist = new TH2D("HoughHist", "Hough Space", 2*m_rhoBins, -m_rhoBins+0.5, m_rhoBins+0.5, m_thetaBins, 0., m_thetaBins);
        
        for(const auto& rhoThetaMap : rhoThetaAccumulatorBinMap)
        {
            houghHist->Fill(rhoThetaMap.first.first, rhoThetaMap.first.second+0.5, rhoThetaMap.second.getAccumulatorValues().size());
        }
        
        houghHist->SetBit(kCanDelete);
        houghHist->Draw();
        m_Canvases.back()->Update();
    }
    
    // **********************************************************************
    // Part II: Use DBScan (or a slight variation) to find clusters of bins
    // **********************************************************************
    
    size_t thresholdLo = std::max(size_t(m_hiThresholdFrac*nAccepted3DHits), m_hiThresholdMin);
    size_t thresholdHi = m_loThresholdFrac * maxBinCount;
    
    std::list<RhoThetaAccumulatorBinMap::iterator> binIndexList;
    
    for(RhoThetaAccumulatorBinMap::iterator mapItr  = rhoThetaAccumulatorBinMap.begin();
        mapItr != rhoThetaAccumulatorBinMap.end();
        mapItr++)
        binIndexList.push_back(mapItr);
    
    binIndexList.sort(SortBinIndexList());
    
    for(auto& mapItr : binIndexList)
    {
        // If we have been here before we skip
        //if (mapItr.second.isVisited()) continue;
        if (mapItr->second.isInCluster()) continue;
        
        // Mark this bin as visited
        // Actually, don't mark it since we are double thresholding and don't want it missed
        //mapItr.second.setVisited();
        
        // Make sure over threshold
        if (mapItr->second.getAccumulatorValues().size() < thresholdLo)
        {
            mapItr->second.setNoise();
            continue;
        }
        
        // Set the low threshold to make sure we merge bins that might be either side of a boundary trajectory
        thresholdHi = std::max(size_t(m_loThresholdFrac * mapItr->second.getAccumulatorValues().size()), m_hiThresholdMin);
        
        // Recover our neighborhood
        HoughCluster neighborhood;
        BinIndex     curBin(mapItr->first);
        
        HoughRegionQuery(curBin, rhoThetaAccumulatorBinMap, neighborhood, thresholdHi);
        
        houghClusters.push_back(HoughCluster());
        
        HoughCluster& houghCluster = houghClusters.back();
        
        expandHoughCluster(curBin, neighborhood, houghCluster, rhoThetaAccumulatorBinMap, thresholdHi);
    }
    
    // Sort the clusters using the SortHoughClusterList metric
    if (!houghClusters.empty()) houghClusters.sort(SortHoughClusterList(rhoThetaAccumulatorBinMap));
    
    return;
}
    
bool HoughSeedFinderAlg::buildSeed(reco::HitPairListPtr& seed3DHits, SeedHitPairListPair& seedHitPair) const
{
    if (seed3DHits.size() < m_minimum3DHits) return false;
    
    reco::PrincipalComponents seedFullPca;
    
    m_pcaAlg.PCAAnalysis_3D(seed3DHits, seedFullPca, true);
    
    if (!seedFullPca.getSvdOK()) return false;

    // Use the following to set the 3D doca and arclength for each hit
    m_pcaAlg.PCAAnalysis_calc3DDocas(seed3DHits, seedFullPca);
    
    // Use this info to sort the hits along the principle axis
    //seed3DHits.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
    seed3DHits.sort(SeedFinderAlgBase::Sort3DHitsByAbsArcLen3D());
    
    // The idea here is to search for the first hit that lies "close" to the principle axis
    // At that point we count out n hits to use as the seed
    reco::HitPairListPtr                seedHit3DList;
    std::set<const reco::ClusterHit2D*> seedHitSet;
    double                              aveDocaToAxis = seedFullPca.getAveHitDoca();
    int                                 gapCount(0);
    
    // Now loop through hits to search for a "continuous" block of at least m_numSeed2DHits
    // We'll arrive at that number by collecting 2D hits in an stl set which will keep track of unique occurances
    for(reco::HitPairListPtr::iterator peakBinItr = seed3DHits.begin(); peakBinItr != seed3DHits.end(); peakBinItr++)
    {
        const reco::ClusterHit3D* hit3D = *peakBinItr;
        
        if (hit3D->getDocaToAxis() < m_numAveDocas*aveDocaToAxis)
        {
            // Check if we need to reset because of gap count
            if (gapCount > m_numSkippedHits)
            {
                seedHit3DList.clear();
                seedHitSet.clear();
            }
            
            seedHit3DList.push_back(hit3D);
            
            for(const auto& hit : hit3D->getHits()) seedHitSet.insert(hit);
            
            gapCount = 0;
        }
        else gapCount++;
        
        if (seedHitSet.size() > m_numSeed2DHits) break;
    }
    
    // If not enough hits then we are done
    if (seedHit3DList.size() < m_minimum3DHits) return false;
    
    reco::PrincipalComponents seedPca;

    // Use only the "seed" 3D hits to get new PCA axes
    m_pcaAlg.PCAAnalysis_3D(seedHit3DList, seedPca, true);
    
    if (!seedPca.getSvdOK()) return false;
    
    m_pcaAlg.PCAAnalysis_calc3DDocas(seedHit3DList, seedPca);
    //seedHit3DList.sort(SeedFinderAlgBase::Sort3DHitsByAbsArcLen3D());
    seedHit3DList.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
    
    // Now translate the seedCenter by the arc len to the first hit
    double seedCenter[3] = {seedPca.getAvePosition()[0],     seedPca.getAvePosition()[1],     seedPca.getAvePosition()[2]};
    double seedDir[3]    = {seedPca.getEigenVectors()[0][0], seedPca.getEigenVectors()[0][1], seedPca.getEigenVectors()[0][2]};
    
    double arcLen       = seedHit3DList.front()->getArclenToPoca();
    double seedStart[3] = {seedCenter[0]+arcLen*seedDir[0], seedCenter[1]+arcLen*seedDir[1], seedCenter[2]+arcLen*seedDir[2]};
    
    //seedStart[0] = seedHit3DList.front()->getX();
    //seedStart[1] = seedHit3DList.front()->getY();
    //seedStart[2] = seedHit3DList.front()->getZ();
    
    if (seedHitSet.size() >= 10)
    {
        TVector3 newSeedPos;
        TVector3 newSeedDir;
        double   chiDOF;
    
        LineFit2DHits(seedHitSet, seedStart[0], newSeedPos, newSeedDir, chiDOF);
    
        if (chiDOF > 0.)
        {
            // check angles between new/old directions
            double cosAng = seedDir[0]*newSeedDir[0]+seedDir[1]*newSeedDir[1]+seedDir[2]*newSeedDir[2];
            
            if (cosAng < 0.) newSeedDir *= -1.;
            
            seedStart[0] = newSeedPos[0];
            seedStart[1] = newSeedPos[1];
            seedStart[2] = newSeedPos[2];
            seedDir[0]   = newSeedDir[0];
            seedDir[1]   = newSeedDir[1];
            seedDir[2]   = newSeedDir[2];
        }
    }

    // Keep track of this seed and the 3D hits that make it up
    seedHitPair = SeedHitPairListPair(recob::Seed(seedStart, seedDir), seedHit3DList);
    
    // We are going to drop a few hits off the ends in the hope this facilitates finding more track like objects, provided there are enough hits
    if (seed3DHits.size() > 100)
    {
        // Need to reset the doca/arclen first
        m_pcaAlg.PCAAnalysis_calc3DDocas(seed3DHits, seedFullPca);
        
        // Now resort the hits
        seed3DHits.sort(SeedFinderAlgBase::Sort3DHitsByAbsArcLen3D());
        
        size_t numToKeep = seed3DHits.size() - 10;
        
        reco::HitPairListPtr::iterator endItr = seed3DHits.begin();
        
        std::advance(endItr, numToKeep);
        
        seed3DHits.erase(endItr, seed3DHits.end());
    }
    
    return true;
}
    
    
bool HoughSeedFinderAlg::findTrackSeeds(reco::HitPairListPtr&      inputHitPairListPtr,
                                        reco::PrincipalComponents& inputPCA,
                                        SeedHitPairListPairVec&    seedHitPairVec) const
{
    // This will be a busy routine... the basic tasks are:
    // 1) loop through hits and project to the plane defined by the two largest eigen values, accumulate in Hough space
    // 2) "Cluster" the Hough space to associate hits which are common to a line
    // 3) Process these clusters (still to be defined exactly)
    
    // Create an interim data structure which will allow us to sort our seeds by "best"
    // before we return them in the seedHitPairVec
    typedef std::map<size_t, SeedHitPairListPairVec > SizeToSeedToHitMap;

    SizeToSeedToHitMap  seedHitPairMap;
    SeedHitPairListPair seedHitPair;
    
    // Make sure we are using the right pca
    reco::HitPairListPtr hitPairListPtr = inputHitPairListPtr;
    
    int nLoops(0);
    
    // Make a local copy of the input PCA
    reco::PrincipalComponents pca = inputPCA;
    
    // We loop over hits in our list until there are no more
    while(!hitPairListPtr.empty())
    {
        // We also require that there be some spread in the data, otherwise not worth running?
        double eigenVal0 = 3. * sqrt(pca.getEigenValues()[0]);
        double eigenVal1 = 3. * sqrt(pca.getEigenValues()[1]);
        
        if (eigenVal0 > 5. && eigenVal1 > 0.001)
        {
            // **********************************************************************
            // Part I: Build Hough space and find Hough clusters
            // **********************************************************************
            RhoThetaAccumulatorBinMap rhoThetaAccumulatorBinMap;
            HoughClusterList          houghClusters;
            
            findHoughClusters(hitPairListPtr, pca, nLoops, rhoThetaAccumulatorBinMap, houghClusters);
            
            // If no clusters then done
            if (houghClusters.empty()) break;
            
            // **********************************************************************
            // Part II: Go through the clusters to find the peak bins
            // **********************************************************************
            
            // We need to use a set so we can be sure to have unique hits
            reco::HitPairListPtr                clusterHitsList;
            std::set<const reco::ClusterHit3D*> masterHitPtrList;
            std::set<const reco::ClusterHit3D*> peakBinPtrList;
            
            size_t firstPeakCount(0);
        
            // Loop through the list of all clusters found above
            for(auto& houghCluster : houghClusters)
            {
                BinIndex peakBin   = houghCluster.front();
                size_t   peakCount = 0;
                size_t   totalHits = 0;
                
                // Make a local (to this cluster) set of of hits
                std::set<const reco::ClusterHit3D*> localHitPtrList;
            
                // Now loop through the bins that were attached to this cluster
                for(auto& binIndex : houghCluster)
                {
                    // An even more local list so we can keep track of peak values
                    std::set<const reco::ClusterHit3D*> tempHitPtrList;
                    
                    // Recover the hits associated to this cluster
                    for(auto& hitItr : rhoThetaAccumulatorBinMap[binIndex].getAccumulatorValues())
                    {
                        reco::HitPairListPtr::const_iterator hit3DItr = hitItr.getHitIterator();
                        
                        tempHitPtrList.insert(*hit3DItr);
                    }
                    
                    // count hits before we remove any
                    totalHits += tempHitPtrList.size();
                    
                    // Trim out any hits already used by a bigger/better cluster
                    std::set<const reco::ClusterHit3D*> tempHit3DSet;
                    
                    std::set_difference(tempHitPtrList.begin(),   tempHitPtrList.end(),
                                        masterHitPtrList.begin(), masterHitPtrList.end(),
                                        std::inserter(tempHit3DSet, tempHit3DSet.end()) );
                    
                    tempHitPtrList = tempHit3DSet;
                    
                    size_t binCount = tempHitPtrList.size();
                
                    if (peakCount < binCount)
                    {
                        peakCount      = binCount;
                        peakBin        = binIndex;
                        peakBinPtrList = tempHitPtrList;
                    }
                    
                    // Add this to our local list
                    localHitPtrList.insert(tempHitPtrList.begin(),tempHitPtrList.end());
                }
                
                if (localHitPtrList.size() < m_minimum3DHits) continue;
                
                if (!firstPeakCount) firstPeakCount = peakCount;
                
                // If the peak counts are significantly less than the first cluster's peak then skip
                if (peakCount < firstPeakCount / 10) continue;
            
                // **********************************************************************
                // Part III: Make a Seed from the peak bin hits
                // **********************************************************************
            
                reco::HitPairListPtr allPeakBinHits;
            
                for(const auto& hit3D : localHitPtrList) allPeakBinHits.push_back(hit3D);
                
                reco::HitPairListPtr peakBinHits;
                
                // Find longest "continuous" set of hits and use these for the seed
                findHitGaps(allPeakBinHits, peakBinHits);
                
                if (peakBinHits.size() < m_minimum3DHits) continue;
                
                // We now build the actual seed.
                if (buildSeed(peakBinHits, seedHitPair))
                {
                    // Keep track of this in our map (which will do ordering for us)
                    seedHitPairMap[peakBinHits.size()].push_back(seedHitPair);
                
                    // For visual testing in event display, mark all the hits in the first seed so we can see them
                    if (seedHitPairMap.size() == 1)
                    {
                        for(const auto& hit3D : peakBinHits) hit3D->setStatusBit(0x40000000);
                    }
                    
//                    for(const auto& hit3D : seedHitPair.second) hit3D->setStatusBit(0x40000000);
                }
                
                // Our peakBinHits collection will most likely be a subset of the localHitPtrList collection
                // We want to remove only the "pure" hits which are those in the peakBinHits collection
                // So, sort them and then add to our master list
                peakBinHits.sort();
                
                masterHitPtrList.insert(peakBinHits.begin(),peakBinHits.end());
                
                if (hitPairListPtr.size() - masterHitPtrList.size() < m_minimum3DHits) break;
            } // end loop over hough clusters
            
            // If the masterHitPtrList is empty then nothing happened and we're done
            if (masterHitPtrList.empty()) break;
            
            // **********************************************************************
            // Part IV: Remove remaining peak bin hits from HitPairPtrList
            // **********************************************************************

            hitPairListPtr.sort();
            
            reco::HitPairListPtr::iterator newListEnd =
                std::set_difference(hitPairListPtr.begin(),   hitPairListPtr.end(),
                                    masterHitPtrList.begin(), masterHitPtrList.end(),
                                    hitPairListPtr.begin() );

            hitPairListPtr.erase(newListEnd, hitPairListPtr.end());
            
            if (hitPairListPtr.size() < m_minimum3DHits) break;
            
            if (nLoops++ > m_maxLoopsPerCluster) break;
            
            // ********************************************************
        }
        else break; // eigen values not in range
        
        // At this point run the PCA on the remaining hits
        m_pcaAlg.PCAAnalysis_3D(hitPairListPtr, pca, true);
        
        if (!pca.getSvdOK()) break;
    }
    
    // The final task before returning is to transfer the stored seeds into the output seed vector
    // What we want to do is make sure the first seeds are the "best" seeds which is defined as the
    // seeds which were associated to the most hits by the Hough Transform. Our seed map will have
    // the reverse of this ordering so we simply iterate through it "backwards"
    for(SizeToSeedToHitMap::reverse_iterator seedMapItr = seedHitPairMap.rbegin(); seedMapItr != seedHitPairMap.rend(); seedMapItr++)
    {
        for(const auto& seedHitPair : seedMapItr->second)
        {
            seedHitPairVec.emplace_back(seedHitPair);
        }
    }
    
    return true;
}
    
bool HoughSeedFinderAlg::findTrackHits(reco::HitPairListPtr&      inputHitPairListPtr,
                                       reco::PrincipalComponents& inputPCA,
                                       reco::HitPairListPtrList&  hitPairListPtrList) const
{
    // The goal of this routine is run the Hough Transform on the input set of hits
    // and then to return a list of lists of hits which are associated to a given line
    
    // Make sure we are using the right pca
    reco::HitPairListPtr hitPairListPtr = inputHitPairListPtr;
    
    int nLoops(0);
    
    // Make a local copy of the input PCA
    reco::PrincipalComponents pca = inputPCA;
    
    // We also require that there be some spread in the data, otherwise not worth running?
    double eigenVal0 = 3. * sqrt(pca.getEigenValues()[0]);
    double eigenVal1 = 3. * sqrt(pca.getEigenValues()[1]);
    
    if (eigenVal0 > 5. && eigenVal1 > 0.001)
    {
        // **********************************************************************
        // Part I: Build Hough space and find Hough clusters
        // **********************************************************************
        RhoThetaAccumulatorBinMap rhoThetaAccumulatorBinMap;
        HoughClusterList          houghClusters;
        
        findHoughClusters(hitPairListPtr, pca, nLoops, rhoThetaAccumulatorBinMap, houghClusters);
        
        // **********************************************************************
        // Part II: Go through the clusters to find the peak bins
        // **********************************************************************
        
        // We need to use a set so we can be sure to have unique hits
        reco::HitPairListPtr                clusterHitsList;
        std::set<const reco::ClusterHit3D*> masterHitPtrList;
        std::set<const reco::ClusterHit3D*> peakBinPtrList;
        
        size_t firstPeakCount(0);
        
        // Loop through the list of all clusters found above
        for(auto& houghCluster : houghClusters)
        {
            BinIndex peakBin   = houghCluster.front();
            size_t   peakCount = 0;
            size_t   totalHits = 0;
            
            // Make a local (to this cluster) set of of hits
            std::set<const reco::ClusterHit3D*> localHitPtrList;
            
            // Now loop through the bins that were attached to this cluster
            for(auto& binIndex : houghCluster)
            {
                // An even more local list so we can keep track of peak values
                std::set<const reco::ClusterHit3D*> tempHitPtrList;
                
                // Recover the hits associated to this cluster
                for(auto& hitItr : rhoThetaAccumulatorBinMap[binIndex].getAccumulatorValues())
                {
                    reco::HitPairListPtr::const_iterator hit3DItr = hitItr.getHitIterator();
                    
                    tempHitPtrList.insert(*hit3DItr);
                }
                
                // count hits before we remove any
                totalHits += tempHitPtrList.size();
                
                // Trim out any hits already used by a bigger/better cluster
                std::set<const reco::ClusterHit3D*> tempHit3DSet;
                
                std::set_difference(tempHitPtrList.begin(),   tempHitPtrList.end(),
                                    masterHitPtrList.begin(), masterHitPtrList.end(),
                                    std::inserter(tempHit3DSet, tempHit3DSet.end()) );
                
                tempHitPtrList = tempHit3DSet;
                
                size_t binCount = tempHitPtrList.size();
                
                if (peakCount < binCount)
                {
                    peakCount      = binCount;
                    peakBin        = binIndex;
                    peakBinPtrList = tempHitPtrList;
                }
                
                // Add this to our local list
                localHitPtrList.insert(tempHitPtrList.begin(),tempHitPtrList.end());
            }
            
            if (localHitPtrList.size() < m_minimum3DHits) continue;
            
            if (!firstPeakCount) firstPeakCount = peakCount;
            
            // If the peak counts are significantly less than the first cluster's peak then skip
            if (peakCount < firstPeakCount / 10) continue;
            
            // **********************************************************************
            // Part III: Make a list of hits from the total number associated
            // **********************************************************************
            
            hitPairListPtrList.push_back(reco::HitPairListPtr());
            
            hitPairListPtrList.back().resize(localHitPtrList.size());
            std::copy(localHitPtrList.begin(), localHitPtrList.end(), hitPairListPtrList.back().begin());
            
            // We want to remove the hits which have been used from further contention
            masterHitPtrList.insert(localHitPtrList.begin(),localHitPtrList.end());
            
            if (hitPairListPtr.size() - masterHitPtrList.size() < m_minimum3DHits) break;
        } // end loop over hough clusters
    }
    
    return true;
}
    
    
//------------------------------------------------------------------------------
void HoughSeedFinderAlg::LineFit2DHits(std::set<const reco::ClusterHit2D*>& hit2DSet,
                                       double                               XOrigin,
                                       TVector3&                            Pos,
                                       TVector3&                            Dir,
                                       double&                              ChiDOF) const
{
    // The following is lifted from Bruce Baller to try to get better
    // initial parameters for a candidate Seed. It is slightly reworked
    // which is why it is included here instead of used as is.
    //
    // Linear fit using X as the independent variable. Hits to be fitted
    // are passed in the hits vector in a pair form (X, WireID). The
    // fitted track position at XOrigin is returned in the Pos vector.
    // The direction cosines are returned in the Dir vector.
    //
    // SVD fit adapted from $ROOTSYS/tutorials/matrix/solveLinear.C
    // Fit equation is w = A(X)v, where w is a vector of hit wires, A is
    // a matrix to calculate a track projected to a point at X, and v is
    // a vector (Yo, Zo, dY/dX, dZ/dX).
    //
    // Note: The covariance matrix should also be returned
    // B. Baller August 2014
    
    // assume failure
    ChiDOF = -1;
    
    if(hit2DSet.size() < 4) return;
    
    const unsigned int nvars = 4;
    unsigned int       npts  = hit2DSet.size();
    
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    unsigned short iht(0), cstat, tpc, ipl;
    double x, cw, sw, off;
    
    // Loop over unique 2D hits from the input list of 3D hits
    for (const auto& hit : hit2DSet)
    {
        geo::WireID wireID = hit->getHit().WireID();
    
        cstat = wireID.Cryostat;
        tpc   = wireID.TPC;
        ipl   = wireID.Plane;
    
        // get the wire plane offset
        off = m_geometry->WireCoordinate(0, 0, ipl, tpc, cstat);
    
        // get the "cosine-like" component
        cw  = m_geometry->WireCoordinate(1, 0, ipl, tpc, cstat) - off;
    
        // the "sine-like" component
        sw  = m_geometry->WireCoordinate(0, 1, ipl, tpc, cstat) - off;
    
        x = hit->getXPosition() - XOrigin;
        
        A[iht][0] = cw;
        A[iht][1] = sw;
        A[iht][2] = cw * x;
        A[iht][3] = sw * x;
        w[iht]    = wireID.Wire - off;
        
        ++ninpl[ipl];
        
        // need at least two points in a plane
        if(ninpl[ipl] == 2) ++nok;
    
        iht++;
    }
    
    // need at least 2 planes with at least two points
    if(nok < 2) return;
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    
    ChiDOF = 0;
    
    // not enough points to calculate Chisq
    if(npts <= 4) return;
    
    double ypr, zpr, diff;
    
    for (const auto& hit : hit2DSet)
    {
        geo::WireID wireID = hit->getHit().WireID();
        
        cstat = wireID.Cryostat;
        tpc   = wireID.TPC;
        ipl   = wireID.Plane;
        off   = m_geometry->WireCoordinate(0, 0, ipl, tpc, cstat);
        cw    = m_geometry->WireCoordinate(1, 0, ipl, tpc, cstat) - off;
        sw    = m_geometry->WireCoordinate(0, 1, ipl, tpc, cstat) - off;
        x     = hit->getXPosition() - XOrigin;
        ypr   = tVec[0] + tVec[2] * x;
        zpr   = tVec[1] + tVec[3] * x;
        diff  = ypr * cw + zpr * sw - (wireID.Wire - off);
        ChiDOF += diff * diff;
    }

    
    float werr2 = m_geometry->WirePitch() * m_geometry->WirePitch();
    ChiDOF /= werr2;
    ChiDOF /= (float)(npts - 4);
    
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    Dir[0] = 1 / norm;
    Dir[1] = tVec[2] / norm;
    Dir[2] = tVec[3] / norm;
    
    Pos[0] = XOrigin;
    Pos[1] = tVec[0];
    Pos[2] = tVec[1];
    
} // TrkLineFit()
    

} // namespace lar_cluster3d
