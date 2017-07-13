/**
 *  @file   PCASeedFinderAlg.cxx
 * 
 *  @brief  Implementation of the Seed Finder Algorithm
 *
 *          The intent of this algorithm is to take an input list of 3D space points and from those
 *          to find candidate track start points and directions
 */

// The main include
#include "dune/RecoAlgDUNE/Cluster3DAlgs/PCASeedFinderAlg.h"
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

PCASeedFinderAlg::PCASeedFinderAlg(fhicl::ParameterSet const &pset) :
     m_gapDistance(5.),
     m_numSeed2DHits(80),
     m_minAllowedCosAng(0.7),
     m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->reconfigure(pset);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    //    auto const* detectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    m_geometry = &*geometry;
    //    m_detector = detectorProperties->provider();
}

//------------------------------------------------------------------------------------------------------------------------------------------

PCASeedFinderAlg::~PCASeedFinderAlg()
{
}
    
void PCASeedFinderAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_gapDistance      = pset.get<double>("GapDistance",       5.);
    m_numSeed2DHits    = pset.get<size_t>("NumSeed2DHits",     80);
    m_minAllowedCosAng = pset.get<double>("MinAllowedCosAng", 0.7);
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
}    
    
bool PCASeedFinderAlg::findTrackSeeds(reco::HitPairListPtr&      inputHitPairListPtr,
                                      reco::PrincipalComponents& inputPCA,
                                      SeedHitPairListPairVec&    seedHitPairVec) const
{
    bool foundGoodSeed(false);
    
    // Make sure we are using the right pca
    reco::HitPairListPtr hitPairListPtr = inputHitPairListPtr;
    
    // Make a local copy of the input pca
    reco::PrincipalComponents pca = inputPCA;

    // We also require that there be some spread in the data, otherwise not worth running?
    double eigenVal0 = 3. * sqrt(pca.getEigenValues()[0]);
    double eigenVal1 = 3. * sqrt(pca.getEigenValues()[1]);
    
    if (eigenVal0 > 5. && eigenVal1 > 0.001)
    {
        // Presume CR muons will be "downward going"...
        if (pca.getEigenVectors()[0][1] > 0.) pca.flipAxis(0);
        
        // Use the following to set the 3D doca and arclength for each hit
        m_pcaAlg.PCAAnalysis_calc3DDocas(hitPairListPtr, pca);
        
        // Use this info to sort the hits along the principle axis
        hitPairListPtr.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
        //hitPairListPtr.sort(PcaSort3DHitsByAbsArcLen3D());
        
        // Make a local copy of the 3D hits
        reco::HitPairListPtr hit3DList;
        
        hit3DList.resize(hitPairListPtr.size());
        
        std::copy(hitPairListPtr.begin(), hitPairListPtr.end(), hit3DList.begin());
        
        reco::PrincipalComponents seedPca = pca;
        
        if (getHitsAtEnd(hit3DList, seedPca))
        {
            // We can use the same routine to check hits at the opposite end to make sure
            // we have consistency between both ends of the track.
            // So... follow the same general program as above
            reco::HitPairListPtr checkList;
            
            checkList.resize(hitPairListPtr.size());
            
            std::copy(hitPairListPtr.begin(), hitPairListPtr.end(), checkList.begin());
            
            std::reverse(checkList.begin(), checkList.end());
            
            reco::PrincipalComponents checkPca = pca;
            
            if (getHitsAtEnd(checkList, checkPca))
            {
                // We want our seed to be in the middle of our seed points
                // First step to getting there is to reorder the hits along the
                // new pca axis
                m_pcaAlg.PCAAnalysis_calc3DDocas(hit3DList, seedPca);
            
                // Use this info to sort the hits along the principle axis - note by absolute value of arc length
                hit3DList.sort(SeedFinderAlgBase::Sort3DHitsByAbsArcLen3D());
            
                // Now translate the seedCenter by the arc len to the first hit
                double seedDir[3]   = {seedPca.getEigenVectors()[0][0], seedPca.getEigenVectors()[0][1], seedPca.getEigenVectors()[0][2]};
                double seedStart[3] = {hit3DList.front()->getX(), hit3DList.front()->getY(), hit3DList.front()->getZ()};
                
                if (hit3DList.size() > 10)
                {
                    TVector3 newSeedPos;
                    TVector3 newSeedDir;
                    double   chiDOF;
                
                    LineFit2DHits(hit3DList, seedStart[0], newSeedPos, newSeedDir, chiDOF);
                
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
            
                for(const auto& hit3D : hit3DList) hit3D->setStatusBit(0x40000000);
            
                seedHitPairVec.emplace_back(std::pair<recob::Seed, reco::HitPairListPtr>(recob::Seed(seedStart, seedDir), hit3DList));
                
                foundGoodSeed = true;
            }
        }
    }
    
    return foundGoodSeed;
}
    
bool PCASeedFinderAlg::getHitsAtEnd(reco::HitPairListPtr& hit3DList, reco::PrincipalComponents& seedPca) const
{
    bool foundGoodSeedHits(false);
    
    // Use a set to count 2D hits by keeping only a single copy
    std::set<const reco::ClusterHit2D*> hit2DSet;
    
    // Look for numSeed2DHits which are continuous
    double lastArcLen = hit3DList.front()->getArclenToPoca();
    
    reco::HitPairListPtr::iterator startItr = hit3DList.begin();
    reco::HitPairListPtr::iterator lastItr  = hit3DList.begin();
    
    while(++lastItr != hit3DList.end())
    {
        const reco::ClusterHit3D* hit3D  = *lastItr;
        double                    arcLen = hit3D->getArclenToPoca();
        
        if (fabs(arcLen-lastArcLen) > m_gapDistance)
        {
            startItr = lastItr;
            hit2DSet.clear();
        }
        
        for(const auto& hit2D : hit3D->getHits())
        {
            hit2DSet.insert(hit2D);
        }
        
        if (hit2DSet.size() > m_numSeed2DHits) break;
        
        lastArcLen = arcLen;
    }
    
    if (hit2DSet.size() > m_numSeed2DHits)
    {
        if (startItr != hit3DList.begin()) hit3DList.erase(hit3DList.begin(), startItr);
        if (lastItr  != hit3DList.end()  ) hit3DList.erase(lastItr,           hit3DList.end());
        
        // On input, the seedPca will contain the original values so we can recover the original axis now
        TVector3 planeVec0(seedPca.getEigenVectors()[0][0],seedPca.getEigenVectors()[0][1],seedPca.getEigenVectors()[0][2]);
        
        m_pcaAlg.PCAAnalysis_3D(hit3DList, seedPca, true);
        
        if (seedPca.getSvdOK())
        {
            // Still looking to point "down"
            if (seedPca.getEigenVectors()[0][1] > 0.) seedPca.flipAxis(0);
            
            // Check that the seed PCA we have found is consistent with the input PCA
            TVector3 primarySeedAxis(seedPca.getEigenVectors()[0][0],seedPca.getEigenVectors()[0][1],seedPca.getEigenVectors()[0][2]);
            
            double cosAng = primarySeedAxis.Dot(planeVec0);
            
            // If the proposed seed axis is not relatively aligned with the input PCA then
            // we should not be using this method to return seeds. Check that here
            if (cosAng > m_minAllowedCosAng) foundGoodSeedHits = true;
        }
    }
    
    return foundGoodSeedHits;
}
    
//------------------------------------------------------------------------------
void PCASeedFinderAlg::LineFit2DHits(const reco::HitPairListPtr& hit3DList,
                                     double                      XOrigin,
                                     TVector3&                   Pos,
                                     TVector3&                   Dir,
                                     double&                     ChiDOF) const
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
    
    // Make a set of 2D hits so we consider only unique hits
    std::set<const reco::ClusterHit2D*> hit2DSet;
    
    for(const auto& hit3D : hit3DList)
    {
        for(const auto& hit2D : hit3D->getHits()) hit2DSet.insert(hit2D);
    }
    
    if(hit2DSet.size() < 4) return;
    
    const unsigned int nvars = 4;
    unsigned int       npts  = 3 * hit3DList.size();
    
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    unsigned short iht(0), cstat, tpc, ipl;
    double x, cw, sw, off;
    
    // Loop over the 2D hits in the above
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

    /*
    //------------------------------------------------------------------------------
    void PCASeedFinderAlg::LineFit2DHits(const reco::HitPairListPtr& hit3DList,
                                         double                      XOrigin,
                                         TVector3&                   Pos,
                                         TVector3&                   Dir,
                                         double&                     ChiDOF) const
    {
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
        
        if(hit3DList.size() < 4) return;
        
        const unsigned int nvars = 4;
        unsigned int       npts  = 3 * hit3DList.size();
        
        TMatrixD A(npts, nvars);
        // vector holding the Wire number
        TVectorD w(npts);
        unsigned short ninpl[3] = {0};
        unsigned short nok = 0;
        unsigned short iht(0), cstat, tpc, ipl;
        double x, cw, sw, off;
        
        for (const auto& hit3D : hit3DList)
        {
            // Inner loop over the 2D hits in the above
            for (const auto& hit : hit3D->getHits())
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
        
        for(const auto& hit3D : hit3DList)
        {
            for (const auto& hit : hit3D->getHits())
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
    */

} // namespace lar_cluster3d
