/**
 *  @file   Cluster3D_module.cc
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes

#include "dune/RecoAlgDUNE/Cluster3DAlgs/PrincipalComponentsAlg.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// ROOT includes
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

PrincipalComponentsAlg::PrincipalComponentsAlg(fhicl::ParameterSet const &pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PrincipalComponentsAlg::~PrincipalComponentsAlg()
{
}
    
void PrincipalComponentsAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    art::ServiceHandle<geo::Geometry>            geometry;
    
    m_parallel = pset.get<double>("ParallelLines", 0.00001);
    m_geometry = &*geometry;
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}
    
void PrincipalComponentsAlg::getHit2DPocaToAxis(const TVector3&            axisPos,
                                                const TVector3&            axisDir,
                                                const reco::ClusterHit2D* hit2D,
                                                TVector3&                  poca,
                                                double&                    arcLenAxis,
                                                double&                    arcLenWire,
                                                double&                    doca)
{
    // Step one is to set up to determine the point of closest approach of this 2D hit to
    // the cluster's current axis.
    // Get this wire's geometry object
    const geo::WireID&  hitID     = hit2D->getHit().WireID();
    const geo::WireGeo& wire_geom = m_geometry->WireIDToWireGeo(hitID);
        
    // From this, get the parameters of the line for the wire
    double wirePos[3] = {0.,0.,0.};
    TVector3 wireDir(wire_geom.Direction());
        
    wire_geom.GetCenter(wirePos);
        
    // Correct the wire position in x to set to correspond to the drift time
    double hitPeak(hit2D->getHit().PeakTime());
        
    wirePos[0] = m_detector->ConvertTicksToX(hitPeak, hitID.Plane, hitID.TPC, hitID.Cryostat);
        
    // Get a vector from the wire position to our cluster's current average position
    TVector3 wVec(axisPos.X()-wirePos[0], axisPos.Y()-wirePos[1], axisPos.Z()-wirePos[2]);
        
    // Get the products we need to compute the arc lengths to the distance of closest approach
    double a(axisDir.Dot(axisDir));
    double b(axisDir.Dot(wireDir));
    double c(wireDir.Dot(wireDir));
    double d(axisDir.Dot(wVec));
    double e(wireDir.Dot(wVec));
        
    double den(a*c - b*b);
    
    // Parallel lines is a special case
    if (fabs(den) > m_parallel)
    {
        arcLenAxis = (b*e - c*d) / den;
        arcLenWire = (a*e - b*d) / den;
    }
    else
    {
        mf::LogDebug("Cluster3D") << "** Parallel lines, need a solution here" << std::endl;
        arcLenAxis = 0.;
        arcLenWire = 0.;
    }
        
    // Now get the hit position we'll use for the pca analysis
    poca =      TVector3(wirePos[0] + arcLenWire * wireDir[0],
                         wirePos[1] + arcLenWire * wireDir[1],
                         wirePos[2] + arcLenWire * wireDir[2]);
    TVector3 axisPocaPos(axisPos[0] + arcLenAxis  * axisDir[0],
                         axisPos[1] + arcLenAxis  * axisDir[1],
                         axisPos[2] + arcLenAxis  * axisDir[2]);
    
    double deltaX(poca.X() - axisPocaPos.X());
    double deltaY(poca.Y() - axisPocaPos.Y());
    double deltaZ(poca.Z() - axisPocaPos.Z());
    double doca2(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    
    doca = sqrt(doca2);
    
    return;
}


struct Sort3DHitsByDocaToAxis
{
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        return left->getDocaToAxis() < right->getDocaToAxis();
    }
    
};
    
struct Sort3DHitsByArcLen3D
{
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        return left->getArclenToPoca() < right->getArclenToPoca();
    }
    
};

struct Sort3DHitsByAbsArcLen3D
{
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
    {
        return fabs(left->getArclenToPoca()) < fabs(right->getArclenToPoca());
    }
    
};
    
void PrincipalComponentsAlg::PCAAnalysis(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, double doca3DScl) const
{
    // This is the controlling outside function for running
    // a Principal Components Analysis on the hits in our
    // input cluster.
    // There is a bootstrap process to be followed
    // 1) Get the initial results from all 3D hits associated
    //    with the cluster
    // 2) Refine the axis by using only the 2D hits on wires

    // Get the first axis from all 3D hits
    PCAAnalysis_3D(hitPairVector, pca);
    
    // Make sure first pass was good
    if (!pca.getSvdOK()) return;

    // First attempt to refine it using only 2D information
    reco::PrincipalComponents pcaLoop = pca;

    PCAAnalysis_2D(hitPairVector, pcaLoop);
    
    // If valid result then go to next steps
    if (pcaLoop.getSvdOK())
    {
        // Let's check the angle between the original and the updated axis
        double cosAngle = pcaLoop.getEigenVectors()[0][0] * pca.getEigenVectors()[0][0]
                        + pcaLoop.getEigenVectors()[0][1] * pca.getEigenVectors()[0][1]
                        + pcaLoop.getEigenVectors()[0][2] * pca.getEigenVectors()[0][2];
        
        // Set the scale factor for the outlier rejection
        double sclFctr(3.);
        
        // If we had a significant change then let's do some outlier rejection, etc.
        if (cosAngle < 1.0)   // pretty much everyone takes a turn
        {
            int    maxIterations(3);
            double maxRange = 3.*sqrt(pcaLoop.getEigenValues()[1]);
            double aveDoca  = pcaLoop.getAveHitDoca();                 // was 0.2
            
            maxRange = sclFctr * 0.5*(maxRange+aveDoca); // was std::max(maxRange, aveDoca);
            
            int numRejHits = PCAAnalysis_reject2DOutliers(hitPairVector, pcaLoop, maxRange);
            int totalRejects(numRejHits);
            int maxRejects(0.4*hitPairVector.size());
            
            // Try looping to see if we improve things
            while(maxIterations-- && numRejHits > 0 && totalRejects < maxRejects)
            {
                // Run the PCA
                PCAAnalysis_2D(hitPairVector, pcaLoop, true);
                
                maxRange = sclFctr * 0.5*(3.*sqrt(pcaLoop.getEigenValues()[1])+pcaLoop.getAveHitDoca());
                
                numRejHits = PCAAnalysis_reject2DOutliers(hitPairVector, pcaLoop, maxRange);
            }

        }
        
        // Ok at this stage copy the latest results back into the cluster
        pca = pcaLoop;

        // Now we make one last pass through the 3D hits to reject outliers there
        PCAAnalysis_reject3DOutliers(hitPairVector, pca, doca3DScl * pca.getAveHitDoca());
    }
    else pca = pcaLoop;
    
    return;
}
    
void PrincipalComponentsAlg::PCAAnalysis_3D(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, bool skeletonOnly) const
{
    // We want to run a PCA on the input TkrVecPoints...
    // The steps are:
    // 1) do a mean normalization of the input vec points
    // 2) compute the covariance matrix
    // 3) run the SVD
    // 4) extract the eigen vectors and values
    // see what happens
    
    // Run through the HitPairList and get the mean position of all the hits
    double meanPos[] = {0.,0.,0.};
    int    numPairsInt(0);
    
    for (const auto& hit : hitPairVector)
    {
        if (skeletonOnly && !((hit->getStatusBits() & 0x10000000) == 0x10000000)) continue;
        
        meanPos[0] += hit->getPosition()[0];
        meanPos[1] += hit->getPosition()[1];
        meanPos[2] += hit->getPosition()[2];
        numPairsInt++;
    }
    
    double numPairs = double(numPairsInt);
    
    meanPos[0] /= numPairs;
    meanPos[1] /= numPairs;
    meanPos[2] /= numPairs;
    
    // Define elements of our covariance matrix
    double xi2  = 0.;
    double xiyi = 0.;
    double xizi = 0.;
    double yi2  = 0.;
    double yizi = 0.;
    double zi2  = 0.;
    
    // Back through the hits to build the matrix
    for (const auto& hit : hitPairVector)
    {
        if (skeletonOnly && !((hit->getStatusBits() & 0x10000000) == 0x10000000)) continue;
        
        double x = hit->getPosition()[0] - meanPos[0];
        double y = hit->getPosition()[1] - meanPos[1];
        double z = hit->getPosition()[2] - meanPos[2];
        
        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z * z;
    }
    
    // Create the actual matrix
    TMatrixD sigma(3, 3);
    
    sigma(0,0) = xi2;
    sigma(0,1) = sigma(1,0) = xiyi;
    sigma(0,2) = sigma(2,0) = xizi;
    sigma(1,1) = yi2;
    sigma(1,2) = sigma(2,1) = yizi;
    sigma(2,2) = zi2;
    
    // Scale by number of pairs
    sigma *= (1./(numPairs - 1.));
    
    // Set up the SVD
    TDecompSVD rootSVD(sigma);
    
    // run the decomposition
    bool svdOk(false);
    
    try
    {
        svdOk = rootSVD.Decompose();
    }
    catch(...)
    {
        svdOk = false;
        mf::LogDebug("Cluster3D") << "PCA decompose failure, numPairs = " << numPairs << std::endl;
    }
    
    if (svdOk)
    {
        // Extract results
        TVectorD eigenVals = rootSVD.GetSig();
        TMatrixD eigenVecs = rootSVD.GetU();
        
        // Get the eigen values
        double recobEigenVals[] = {eigenVals[0], eigenVals[1], eigenVals[2]};
        
        // Grab the principle axes
        reco::PrincipalComponents::EigenVectors recobEigenVecs;
        std::vector<double> tempVec;
        
        // Get the first column vector
        tempVec.push_back(eigenVecs(0, 0));
        tempVec.push_back(eigenVecs(1, 0));
        tempVec.push_back(eigenVecs(2, 0));
        recobEigenVecs.push_back(tempVec);
        
        // Now the second
        tempVec.clear();
        tempVec.push_back(eigenVecs(0, 1));
        tempVec.push_back(eigenVecs(1, 1));
        tempVec.push_back(eigenVecs(2, 1));
        recobEigenVecs.push_back(tempVec);
        
        // And the last
        tempVec.clear();
        tempVec.push_back(eigenVecs(0, 2));
        tempVec.push_back(eigenVecs(1, 2));
        tempVec.push_back(eigenVecs(2, 2));
        recobEigenVecs.push_back(tempVec);
        
        // Store away
        pca = reco::PrincipalComponents(svdOk, numPairsInt, recobEigenVals, recobEigenVecs, meanPos);
    }
    else
    {
        pca = reco::PrincipalComponents();
    }
    
    return;
}
    
void PrincipalComponentsAlg::PCAAnalysis_2D(const reco::HitPairListPtr& hitPairVector, reco::PrincipalComponents& pca, bool updateAvePos) const
{
    // Once an axis has been found our goal is to refine it by using only the 2D hits
    // We'll get 3D information for each of these by using the axis as a reference and use
    // the point of closest approach as the 3D position
    
    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double xizi(0.);
    double yi2(0.);
    double yizi(0.);
    double zi2(0.);
    
    double   aveHitDoca(0.);
    TVector3 avePosUpdate(0.,0.,0.);
    int      nHits(0);
    
    // Recover existing line parameters for current cluster
    const reco::PrincipalComponents& inputPca      = pca;
    TVector3                         avePosition(inputPca.getAvePosition()[0], inputPca.getAvePosition()[1], inputPca.getAvePosition()[2]);
    TVector3                         axisDirVec(inputPca.getEigenVectors()[0][0], inputPca.getEigenVectors()[0][1], inputPca.getEigenVectors()[0][2]);
    
    // We double loop here so we can use this method for both the first time through
    // and a second time through where we re-calculate the mean position
    // So, we need to keep track of the poca which we do with a double vector
    std::vector<TVector3> hitPosVec;
    
    // Outer loop over 3D hits
    for (const auto& hit3D : hitPairVector)
    {
        // Inner loop over 2D hits
        for (const auto& hit : hit3D->getHits())
        {
            // Step one is to set up to determine the point of closest approach of this 2D hit to
            // the cluster's current axis.
            // Get this wire's geometry object
            const geo::WireID&  hitID     = hit->getHit().WireID();
            const geo::WireGeo& wire_geom = m_geometry->WireIDToWireGeo(hitID);
            
            // From this, get the parameters of the line for the wire
            double wirePosArr[3] = {0.,0.,0.};
            wire_geom.GetCenter(wirePosArr);

            TVector3 wireCenter(wirePosArr[0], wirePosArr[1], wirePosArr[2]);
            TVector3 wireDirVec(wire_geom.Direction());
            
            // Correct the wire position in x to set to correspond to the drift time
            double hitPeak(hit->getHit().PeakTime());
            
            TVector3 wirePos(m_detector->ConvertTicksToX(hitPeak, hitID.Plane, hitID.TPC, hitID.Cryostat), wireCenter[1], wireCenter[2]);
            
            // Compute the wire plane normal for this view
            TVector3 xAxis(1.,0.,0.);
            TVector3 planeNormal = xAxis.Cross(wireDirVec);   // This gives a normal vector in +z for a Y wire

            double docaInPlane(wirePos[0] - avePosition[0]);
            double arcLenToPlane(0.);
            double cosAxisToPlaneNormal = axisDirVec.Dot(planeNormal);
            
            TVector3 axisPlaneIntersection = wirePos;
            TVector3 hitPosTVec            = wirePos;
            
            if (fabs(cosAxisToPlaneNormal) > 0.)
            {
                TVector3 deltaPos = wirePos - avePosition;
                
                arcLenToPlane         = deltaPos.Dot(planeNormal) / cosAxisToPlaneNormal;
                axisPlaneIntersection = avePosition + arcLenToPlane * axisDirVec;
                docaInPlane           = wirePos[0] - axisPlaneIntersection[0];
                
                TVector3 axisToInter  = axisPlaneIntersection - wirePos;
                double   arcLenToDoca = axisToInter.Dot(wireDirVec);
                
                hitPosTVec += arcLenToDoca * wireDirVec;
            }

            // Get a vector from the wire position to our cluster's current average position
            TVector3 wVec = avePosition - wirePos;
            
            // Get the products we need to compute the arc lengths to the distance of closest approach
            double a(axisDirVec.Dot(axisDirVec));
            double b(axisDirVec.Dot(wireDirVec));
            double c(wireDirVec.Dot(wireDirVec));
            double d(axisDirVec.Dot(wVec));
            double e(wireDirVec.Dot(wVec));
            
            double den(a*c - b*b);
            double arcLen1(0.);
            double arcLen2(0.);
            
            // Parallel lines is a special case
            if (fabs(den) > m_parallel)
            {
                arcLen1 = (b*e - c*d) / den;
                arcLen2 = (a*e - b*d) / den;
            }
            else
            {
                mf::LogDebug("Cluster3D") << "** Parallel lines, need a solution here" << std::endl;
                break;
            }
            
            // Now get the hit position we'll use for the pca analysis
            //double hitPos[]  = {wirePos[0]     + arcLen2 * wireDirVec[0],
            //                    wirePos[1]     + arcLen2 * wireDirVec[1],
            //                    wirePos[2]     + arcLen2 * wireDirVec[2]};
            //double axisPos[] = {avePosition[0] + arcLen1 * axisDirVec[0],
            //                    avePosition[1] + arcLen1 * axisDirVec[1],
            //                    avePosition[2] + arcLen1 * axisDirVec[2]};
            TVector3 hitPos  = wirePos + arcLen2 * wireDirVec;
            TVector3 axisPos = avePosition + arcLen1 * axisDirVec;
            double   deltaX  = hitPos[0] - axisPos[0];
            double   deltaY  = hitPos[1] - axisPos[1];
            double   deltaZ  = hitPos[2] - axisPos[2];
            double   doca2   = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            double   doca    = sqrt(doca2);
            
            docaInPlane = doca;

            aveHitDoca += fabs(docaInPlane);

            //TVector3 deltaPos  = hitPos - hitPosTVec;
            //double   deltaDoca = doca - docaInPlane;
            
            //if (fabs(deltaPos[0]) > 1. || fabs(deltaPos[1]) > 1. || fabs(deltaPos[2]) > 1. || fabs(deltaDoca) > 2.)
            //{
            //    std::cout << "**************************************************************************************" << std::endl;
            //    std::cout << "Diff in doca: " << deltaPos[0] << "," << deltaPos[1] << "," << deltaPos[2] << ", deltaDoca: " << deltaDoca << std::endl;
            //    std::cout << "-- HitPosTVec: " << hitPosTVec[0] << "," << hitPosTVec[1] << "," << hitPosTVec[2] << std::endl;
            //    std::cout << "-- hitPos:     " << hitPos[0] << "," << hitPos[1] << "," << hitPos[2] << std::endl;
            //    std::cout << "-- WirePos:    " << wirePos[0] << "," << wirePos[1] << "," << wirePos[2] << std::endl;
            //}

            // Set the hit's doca and arclen
            hit->setDocaToAxis(fabs(docaInPlane));
            hit->setArcLenToPoca(arcLenToPlane);
            
            // If this point is considered an outlier then we skip
            // the accumulation for average position and covariance
            if (hit->getStatusBits() & 0x80)
            {
                continue;
            }

            //hitPosTVec = hitPos;
            
            avePosUpdate += hitPosTVec;
            
            hitPosVec.push_back(hitPosTVec);
            
            nHits++;
        }
    }
    
    // Get updated average position
    avePosUpdate[0] /= double(nHits);
    avePosUpdate[1] /= double(nHits);
    avePosUpdate[2] /= double(nHits);
    
    // Get the average hit doca
    aveHitDoca /= double(nHits);
    
    if (updateAvePos)
    {
        avePosition = avePosUpdate;
    }
    
    // Now loop through the hits and build out the covariance matrix
    for(auto& hitPos : hitPosVec)
    {
        // And increment the values in the covariance matrix
        double x = hitPos[0] - avePosition[0];
        double y = hitPos[1] - avePosition[1];
        double z = hitPos[2] - avePosition[2];
            
        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z * z;
    }
    
    // Accumulation done, now do the actual work
               
    // Create the actual matrix
    TMatrixD sigma(3, 3);
               
    sigma(0,0) = xi2;
    sigma(0,1) = sigma(1,0) = xiyi;
    sigma(0,2) = sigma(2,0) = xizi;
    sigma(1,1) = yi2;
    sigma(1,2) = sigma(2,1) = yizi;
    sigma(2,2) = zi2;
               
    // Scale by number of pairs
    sigma *= (1./(nHits - 1.));
               
    // Set up the SVD
    TDecompSVD rootSVD(sigma);
               
    // run the decomposition
    bool svdOk(false);
    
    try
    {
        svdOk = rootSVD.Decompose();
    }
    catch(...)
    {
        svdOk = false;
        mf::LogDebug("Cluster3D") << "PCA decompose failure, nhits = " << nHits << std::endl;
    }
    
    if (svdOk)
    {
        // Extract results
        TVectorD eigenVals = rootSVD.GetSig();
        TMatrixD eigenVecs = rootSVD.GetU();
                   
        // Get the eigen values
        double recobEigenVals[] = {eigenVals[0], eigenVals[1], eigenVals[2]};
                   
        // Grab the principle axes
        reco::PrincipalComponents::EigenVectors recobEigenVecs;
        std::vector<double> tempVec;
                   
        // Get the first column vector
        tempVec.push_back(eigenVecs(0, 0));
        tempVec.push_back(eigenVecs(1, 0));
        tempVec.push_back(eigenVecs(2, 0));
        recobEigenVecs.push_back(tempVec);
                   
        // Now the second
        tempVec.clear();
        tempVec.push_back(eigenVecs(0, 1));
        tempVec.push_back(eigenVecs(1, 1));
        tempVec.push_back(eigenVecs(2, 1));
        recobEigenVecs.push_back(tempVec);
                   
        // And the last
        tempVec.clear();
        tempVec.push_back(eigenVecs(0, 2));
        tempVec.push_back(eigenVecs(1, 2));
        tempVec.push_back(eigenVecs(2, 2));
        recobEigenVecs.push_back(tempVec);
        
        // Save the average position
        double avePosToSave[] = {avePosition[0],avePosition[1],avePosition[2]};
        
        // Store away
        pca = reco::PrincipalComponents(svdOk, nHits, recobEigenVals, recobEigenVecs, avePosToSave, aveHitDoca);
    }
    else
    {
        pca = reco::PrincipalComponents();
    }
    
    return;
}
    
void PrincipalComponentsAlg::PCAAnalysis_calc3DDocas(const reco::HitPairListPtr& hitPairVector,
                                                     const reco::PrincipalComponents& pca) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.
    
    // We'll need the current PCA axis to determine doca and arclen
    TVector3 avePosition(pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);
    TVector3 axisDirVec(pca.getEigenVectors()[0][0], pca.getEigenVectors()[0][1], pca.getEigenVectors()[0][2]);
    
    // We want to keep track of the average
    double aveDoca3D(0.);
    
    // Outer loop over views
    for (const auto* clusterHit3D : hitPairVector)
    {
        // Always reset the existing status bit
        clusterHit3D->clearStatusBits(0x80);
        
        // Now we need to calculate the doca and poca...
        // Start by getting this hits position
        TVector3 clusPos(clusterHit3D->getPosition()[0],clusterHit3D->getPosition()[1],clusterHit3D->getPosition()[2]);
        
        // Form a TVector from this to the cluster average position
        TVector3 clusToHitVec = clusPos - avePosition;
        
        // With this we can get the arclength to the doca point
        double arclenToPoca = clusToHitVec.Dot(axisDirVec);
        
        // Get the coordinates along the axis for this point
        TVector3 docaPos = avePosition + arclenToPoca * axisDirVec;
        
        // Now get doca and poca
        TVector3 docaPosToClusPos = clusPos - docaPos;
        double   docaToAxis       = docaPosToClusPos.Mag();
        
        aveDoca3D += docaToAxis;
        
        // Ok, set the values in the hit
        clusterHit3D->setDocaToAxis(docaToAxis);
        clusterHit3D->setArclenToPoca(arclenToPoca);
    }
    
    // Compute the average and store
    aveDoca3D /= double(hitPairVector.size());
    
    pca.setAveHitDoca(aveDoca3D);
    
    return;
}
    
void PrincipalComponentsAlg::PCAAnalysis_calc2DDocas(const reco::Hit2DListPtr&        hit2DListPtr,
                                                     const reco::PrincipalComponents& pca         ) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.
    
    // We'll need the current PCA axis to determine doca and arclen
    TVector3 avePosition(pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);
    TVector3 axisDirVec(pca.getEigenVectors()[0][0], pca.getEigenVectors()[0][1], pca.getEigenVectors()[0][2]);
    
    // Recover the principle eigen value for range constraints
    double maxArcLen = 4.*sqrt(pca.getEigenValues()[0]);

    // We want to keep track of the average
    double aveHitDoca(0.);
    
    // Outer loop over views
    for (const auto* hit : hit2DListPtr)
    {
        // Step one is to set up to determine the point of closest approach of this 2D hit to
        // the cluster's current axis. We do that by finding the point of intersection of the
        // cluster's axis with a plane defined by the wire the hit is associated with.
        // Get this wire's geometry object
        const geo::WireID&  hitID     = hit->getHit().WireID();
        const geo::WireGeo& wire_geom = m_geometry->WireIDToWireGeo(hitID);
        
        // From this, get the parameters of the line for the wire
        double wirePosArr[3] = {0.,0.,0.};
        wire_geom.GetCenter(wirePosArr);
        
        TVector3 wireCenter(wirePosArr[0], wirePosArr[1], wirePosArr[2]);
        TVector3 wireDirVec(wire_geom.Direction());
        
        // Correct the wire position in x to set to correspond to the drift time
        TVector3 wirePos(hit->getXPosition(), wireCenter[1], wireCenter[2]);

        // Compute the wire plane normal for this view
        TVector3 xAxis(1.,0.,0.);
        TVector3 planeNormal = xAxis.Cross(wireDirVec);   // This gives a normal vector in +z for a Y wire
        
        double arcLenToPlane(0.);
        double docaInPlane(wirePos[0] - avePosition[0]);
        double cosAxisToPlaneNormal = axisDirVec.Dot(planeNormal);
        
        TVector3 axisPlaneIntersection = wirePos;

        // If current cluster axis is not parallel to wire plane then find intersection point
        if (fabs(cosAxisToPlaneNormal) > 0.)
        {
            TVector3 deltaPos = wirePos - avePosition;
            
            arcLenToPlane         = std::min(deltaPos.Dot(planeNormal) / cosAxisToPlaneNormal, maxArcLen);
            axisPlaneIntersection = avePosition + arcLenToPlane * axisDirVec;
            TVector3 axisToInter  = axisPlaneIntersection - wirePos;
            double   arcLenToDoca = axisToInter.Dot(wireDirVec);

            // If the arc length along the wire to the poca is outside the TPC then reset
            if (fabs(arcLenToDoca) > wire_geom.HalfL()) arcLenToDoca = wire_geom.HalfL();

            // If we were successful in getting to the wire plane then the doca is simply the
            // difference in x coordinates... but we hvae to worry about the special cases so
            // we calculate a 3D doca based on arclengths above...
            TVector3 docaVec = axisPlaneIntersection - (wirePos + arcLenToDoca * wireDirVec);
            docaInPlane = docaVec.Mag();
        }
        
        aveHitDoca += fabs(docaInPlane);
        
        // Set the hit's doca and arclen
        hit->setDocaToAxis(fabs(docaInPlane));
        hit->setArcLenToPoca(arcLenToPlane);
    }
    
    // Compute the average and store
    aveHitDoca /= double(hit2DListPtr.size());
    
    pca.setAveHitDoca(aveHitDoca);
    
    return;
}
    
int PrincipalComponentsAlg::PCAAnalysis_reject2DOutliers(const reco::HitPairListPtr& hitPairVector,
                                                         reco::PrincipalComponents&  pca,
                                                         double                      maxDocaAllowed) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.
    
    // First get the average doca scaled by some appropriate factor
    int    numRejHits(0);
    
    // Outer loop over views
    for (const auto& hit3D : hitPairVector)
    {
        // Inner loop over hits in this view
        for (const auto& hit : hit3D->getHits())
        {
            // Always reset the existing status bit
            hit->clearStatusBits(0x80);
            
            if (hit->getDocaToAxis() > maxDocaAllowed)
            {
                hit->setStatusBit(0x80);
                numRejHits++;
            }
        }
    }
    
    return numRejHits;
}
    
int PrincipalComponentsAlg::PCAAnalysis_reject3DOutliers(const reco::HitPairListPtr&      hitPairVector,
                                                         const reco::PrincipalComponents& pca,
                                                         double                           maxDocaAllowed) const
{
    // Our mission, should we choose to accept it, is to scan through the 2D hits and reject
    // any outliers. Basically, any hit outside a scaled range of the average doca from the
    // first pass is marked by setting the bit in the status word.
        
    // First get the average doca scaled by some appropriate factor
    int    numRejHits(0);
    
    // We'll need the current PCA axis to determine doca and arclen
    TVector3 avePosition(pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);
    TVector3 axisDirVec(pca.getEigenVectors()[0][0], pca.getEigenVectors()[0][1], pca.getEigenVectors()[0][2]);
    
    // Outer loop over views
    for (const auto* clusterHit3D : hitPairVector)
    {
        // Always reset the existing status bit
        clusterHit3D->clearStatusBits(0x80);
        
        // Now we need to calculate the doca and poca...
        // Start by getting this hits position
        TVector3 clusPos(clusterHit3D->getPosition()[0],clusterHit3D->getPosition()[1],clusterHit3D->getPosition()[2]);
        
        // Form a TVector from this to the cluster average position
        TVector3 clusToHitVec = clusPos - avePosition;
        
        // With this we can get the arclength to the doca point
        double arclenToPoca = clusToHitVec.Dot(axisDirVec);
        
        // Get the coordinates along the axis for this point
        TVector3 docaPos = avePosition + arclenToPoca * axisDirVec;
        
        // Now get doca and poca
        TVector3 docaPosToClusPos = clusPos - docaPos;
        double   docaToAxis       = docaPosToClusPos.Mag();
        
        // Ok, set the values in the hit
        clusterHit3D->setDocaToAxis(docaToAxis);
        clusterHit3D->setArclenToPoca(arclenToPoca);
        
        // Check to see if this is a keeper
        if (clusterHit3D->getDocaToAxis() > maxDocaAllowed)
        {
            clusterHit3D->setStatusBit(0x80);
            numRejHits++;
        }
    }
        
    return numRejHits;
}

    

} // namespace lar_cluster3d
