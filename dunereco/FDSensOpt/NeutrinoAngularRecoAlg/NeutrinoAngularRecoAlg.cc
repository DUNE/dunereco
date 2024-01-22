/**
*  @file   dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.cc

*  @brief  Implementation file for the neutrino angle reconstruction algorithm.
*  Written by Henrique Souza (hvsouza@apc.in2p3.fr)
*  $Log: $
*/

//STL
#include <limits>
#include <algorithm>
//ROOT
//ART
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
//DUNE
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"

#include "dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h"

namespace dune
{
NeutrinoAngularRecoAlg::NeutrinoAngularRecoAlg(fhicl::ParameterSet const& pset, const std::string &trackLabel, 
    const std::string &showerLabel, const std::string &hitLabel, const std::string wireLabel, 
    const std::string &trackToHitLabel, const std::string &showerToHitLabel, const std::string &hitToSpacePointLabel) :
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fTrackLabel(trackLabel),
    fShowerLabel(showerLabel),
    fHitLabel(hitLabel),
    fWireLabel(wireLabel),
    fTrackToHitLabel(trackToHitLabel),
    fShowerToHitLabel(showerToHitLabel),
    fHitToSpacePointLabel(hitToSpacePointLabel),
    fCalorimetryLabel(pset.get<std::string>("CalorimetryLabel")),
    fDistanceToWallThreshold(pset.get<float>("DistanceToWallThreshold"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event)
{
    if (!pMuonTrack.isAvailable() || pMuonTrack.isNull())
    {
        mf::LogWarning("NeutrinoAngularRecoAlg") << " Cannot access the muon track which is needed for this angular reconstruction method.\n"
        << "Not able to calculate neutrino direction" << std::endl;

        Point3D_t vertex(0,0,0);
        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }

    Point3D_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

    TVector3 dir_start;
    dir_start = pMuonTrack->VertexDirection<TVector3>();

    Direction_t direction(dir_start.X(), dir_start.Y(), dir_start.Z());
    AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kMuon);
    return this->ReturnNeutrinoAngle(angularRecoInputHolder);

}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const art::Ptr<recob::Shower> &pElectronShower, 
    const art::Event &event)
{
    if (!pElectronShower.isAvailable() || pElectronShower.isNull())
    {
        mf::LogWarning("NeutrinoAngularRecoAlg") 
        << " Cannot access the electron shower which is needed for this angular reconstruction method.\n"
        << "Not able to calculate neutrino direction" << std::endl;
        Point3D_t vertex(0,0,0);
        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }


    Point3D_t vertex(pElectronShower->ShowerStart().X(), pElectronShower->ShowerStart().Y(), pElectronShower->ShowerStart().Z());

    TVector3 const& dir_start = pElectronShower->Direction();
    Direction_t direction(dir_start.X(), dir_start.Y(), dir_start.Z());
    AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kElectron);
    return this->ReturnNeutrinoAngle(angularRecoInputHolder);
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const std::vector<art::Ptr<recob::Track>> &pTracks,
                                                                       const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                                                       const std::vector<art::Ptr<recob::Shower>> &pShowers,
                                                                       const art::Event &event)
{
    //Using all the reco particles to determine the angle
    Point3D_t vertex(0,0,0); //No meaning vertex value is filled with this method

    if (pTracks.size() == 0 && pShowers.size() == 0)
    {
        mf::LogWarning("NeutrinoAngularRecoAlg") << " No track or shower selected.\n"
        << "Not able to calculate neutrino direction" << std::endl;

        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }

    Momentum_t mom_showers = ComputeShowersMomentum(pShowers);
    Momentum_t mom_tracks = ComputeTracksMomentum(pTracks, tracksPID, event);
    Momentum_t mom_tot(mom_showers.X() + mom_tracks.X(), mom_showers.Y() + mom_tracks.Y(), mom_showers.Z() + mom_tracks.Z());
    Momentum_t mom_direction = mom_tot/sqrt(mom_tot.Mag2());

    Direction_t direction(mom_direction.X(), mom_direction.Y(), mom_direction.Z());
    AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoParticles);
    return this->ReturnNeutrinoAngle(angularRecoInputHolder);

}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const art::Ptr<recob::Track> &pMuonTrack,
                                                                       const std::vector<art::Ptr<recob::Track>> &pTracks,
                                                                       const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                                                       const std::vector<art::Ptr<recob::Shower>> &pShowers,
                                                                       const art::Event &event)
{
    //Using all the reco particles to determine the angle, but using the longest track as muon
    Point3D_t vertex(0,0,0); //No meaning vertex value is filled with this method

    if (!pMuonTrack.isAvailable() || pMuonTrack.isNull())
    {
        mf::LogWarning("NeutrinoAngularRecoAlg")
        << " Cannot access the muon track which is needed for this angular reconstruction method.\n"
        << "Not able to calculate neutrino direction" << std::endl;

        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }

    //This method is very similar to the one above. The only difference is that we force set the pid of the longest track to 13.
    std::map<art::Ptr<recob::Track>, int> tracksPIDUpdated(tracksPID);
    bool foundMatch = false;
    for(auto &[pTrack, pid] : tracksPIDUpdated){
        if(pTrack == pMuonTrack){
            pid = 13;
            foundMatch = true;
            break;
        }
    }

    if(!foundMatch){
        mf::LogWarning("NeutrinoAngularRecoAlg") << "Did not find a match for the track. This is unexpected!!" << std::endl;
    }

    dune::AngularRecoOutput output = CalculateNeutrinoAngle(pTracks, tracksPIDUpdated, pShowers, event);
    if(output.recoMethodUsed != kRecoMethodNotSet){
        output.recoMethodUsed = kMuonRecoParticles;
    }
    return output;
}

//------------------------------------------------------------------------------------------------------------------------------------------


dune::AngularRecoOutput NeutrinoAngularRecoAlg::ReturnNeutrinoAngle(const AngularRecoInputHolder &angularRecoInputHolder)
{
    dune::AngularRecoOutput output;
    output.recoMethodUsed = angularRecoInputHolder.fAngularRecoMethod;
    output.fRecoVertex = angularRecoInputHolder.fVertex;
    output.fRecoDirection = angularRecoInputHolder.fNuDirection;
    return output;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Momentum_t NeutrinoAngularRecoAlg::ComputeShowersMomentum(const std::vector<art::Ptr<recob::Shower>> &pShowers) const{
    Momentum_t momentum(0, 0, 0);

    for(const art::Ptr<recob::Shower> &pShower : pShowers){
        float E = GetShowerEnergy(pShower)*1e-3; //Converts from MeV to GeV
        const TVector3 &dir = pShower->Direction();
        momentum.SetX(momentum.X() + E*dir.X());
        momentum.SetY(momentum.Y() + E*dir.Y());
        momentum.SetZ(momentum.Z() + E*dir.Z());
    }

    return momentum;
};

//------------------------------------------------------------------------------------------------------------------------------------------

Momentum_t NeutrinoAngularRecoAlg::ComputeTracksMomentum(const std::vector<art::Ptr<recob::Track>> &pTracks,
                                                         const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                                         const art::Event &event) const{
    Momentum_t momentum(0, 0, 0);
    trkf::TrackMomentumCalculator TrackMomCalc;

    art::FindManyP<anab::Calorimetry> fmCal(pTracks, event, fCalorimetryLabel);

    for(uint iTrack = 0; iTrack < pTracks.size(); iTrack++){
        const art::Ptr<recob::Track> &pTrack = pTracks[iTrack];
        int pid = 0;
        if (tracksPID.count(pTrack) == 1){
            pid = tracksPID.at(pTrack);
        }

        float momentum_norm = 0;

        if(!IsTrackContained(pTrack)){
            //TODO: In the future change the momentum evaluation method. Requires GetMomentumMultiScatter to be consolidated first
            // and some checks about its validity for multiple PID particles
        }

        switch (abs(pid))
        {
        case 13:
            momentum_norm = TrackMomCalc.GetTrackMomentum(pTrack->Length(), pid);
            break;
        case 2212:
            momentum_norm = TrackMomCalc.GetTrackMomentum(pTrack->Length(), pid);
            break;
        case 211:
            momentum_norm = GetTrackKE(fmCal.at(iTrack)); //This is in fact only KE for now
            momentum_norm = sqrt(momentum_norm*(2*fPION_MASS + momentum_norm))*1e-3; //Converts MeV to GeV
            break;
        default:
            //If we have no relevant PID assigned but still a track, we assume a massless particle, this should be better than nothing
            momentum_norm = GetTrackKE(fmCal.at(iTrack))*1e-3; //Converts MeV to GeV
            break;
        }
        const recob::Track::Vector_t &dir = pTrack->StartDirection();
        momentum.SetX(momentum.X() + momentum_norm*dir.X());
        momentum.SetY(momentum.Y() + momentum_norm*dir.Y());
        momentum.SetZ(momentum.Z() + momentum_norm*dir.Z());
    }

    return momentum;
};

//------------------------------------------------------------------------------------------------------------------------------------------

float NeutrinoAngularRecoAlg::GetShowerEnergy(const art::Ptr<recob::Shower> &pShower) const{
    //Using the energy in the best plane
    return pShower->Energy().at(pShower->best_plane());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float NeutrinoAngularRecoAlg::GetTrackKE(const std::vector<art::Ptr<anab::Calorimetry>> &calos) const {
    //Estimates the KE from the plane with the highest deposited energy in MeV
    float maxKE = 0;
    for (const art::Ptr<anab::Calorimetry> &cal : calos){
        if(!cal.isAvailable() || cal.isNull())
            continue;
        if(!cal->PlaneID().isValid)
            continue;

        maxKE = std::max(maxKE, cal->KineticEnergy());
    }
    return maxKE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoAngularRecoAlg::IsTrackContained(const art::Ptr<recob::Track> &pTrack) const {
    art::ServiceHandle<geo::Geometry> fGeometry;
    double minX(std::numeric_limits<double>::max());
    double maxX(std::numeric_limits<double>::lowest());
    double minY(std::numeric_limits<double>::max());
    double maxY(std::numeric_limits<double>::lowest());
    double minZ(std::numeric_limits<double>::max());
    double maxZ(std::numeric_limits<double>::lowest());
    for (geo::TPCGeo const& tpc: fGeometry->Iterate<geo::TPCGeo>()) {
        minX = std::min(minX,tpc.MinX());
        maxX = std::max(maxX,tpc.MaxX());
        minY = std::min(minY,tpc.MinY());
        maxY = std::max(maxY,tpc.MaxY());
        minZ = std::min(minZ,tpc.MinZ());
        maxZ = std::max(maxZ,tpc.MaxZ());
    } // for all TPC
    minX += fDistanceToWallThreshold;
    maxX -= fDistanceToWallThreshold;
    minY += fDistanceToWallThreshold;
    maxY -= fDistanceToWallThreshold;
    minZ += fDistanceToWallThreshold;
    maxZ -= fDistanceToWallThreshold;

    TVector3 start = pTrack->Start<TVector3>();
    TVector3 end = pTrack->End<TVector3>();

    if (start.X() - minX < -1.*std::numeric_limits<double>::epsilon() ||
        start.X() - maxX > std::numeric_limits<double>::epsilon() ||
        end.X() - minX < -1.*std::numeric_limits<double>::epsilon() ||
        end.X() - maxX > std::numeric_limits<double>::epsilon()
        )
        return false;
    if (start.Y() - minY < -1.*std::numeric_limits<double>::epsilon() ||
        start.Y() - maxY > std::numeric_limits<double>::epsilon() ||
        end.Y() - minY < -1.*std::numeric_limits<double>::epsilon() ||
        end.Y() - maxY > std::numeric_limits<double>::epsilon())
        return false;
    if (start.Z() - minZ < -1.*std::numeric_limits<double>::epsilon() ||
        start.Z() - maxZ > std::numeric_limits<double>::epsilon() ||
        end.Z() - minZ < -1.*std::numeric_limits<double>::epsilon() ||
        end.Z() - maxZ > std::numeric_limits<double>::epsilon())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} //dune
