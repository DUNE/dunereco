/**
*  @file   dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.cc

*  @brief  Implementation file for the neutrino angle reconstruction algorithm.
*
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
    fHitToSpacePointLabel(hitToSpacePointLabel)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event)
{
    if (!pMuonTrack.isAvailable() || pMuonTrack.isNull())
    {
        mf::LogWarning("NeutrinoAngularRecoAlg") << " Cannot access the muon track which is needed for this angular reconstruction method.\n"
        << "Not able to calculate neutrino direction" << std::endl;

        Point_t vertex(0,0,0);
        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }

    Point_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

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
        Point_t vertex(0,0,0);
        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }


    Point_t vertex(pElectronShower->ShowerStart().X(), pElectronShower->ShowerStart().Y(), pElectronShower->ShowerStart().Z());

    TVector3 const& dir_start = pElectronShower->Direction();
    Direction_t direction(dir_start.X(), dir_start.Y(), dir_start.Z());
    AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kElectron);
    return this->ReturnNeutrinoAngle(angularRecoInputHolder);
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

} //dune
