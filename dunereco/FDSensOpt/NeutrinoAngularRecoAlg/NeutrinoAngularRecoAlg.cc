/**
*  @file   dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.cc

*  @brief  Implementation file for the neutrino angle reconstruction algorithm.
*  Written by Pierre Granger (pierre.granger@cern.ch) & Henrique Souza (hvsouza@apc.in2p3.fr)
*  $Log: $
*/

//STL
#include <limits>
#include <algorithm>
#include <complex>
//ROOT
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
//ART
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
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

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const std::vector<art::Ptr<recob::Track>> &pTracks,
                                                                       const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                                                       const std::vector<art::Ptr<recob::Shower>> &pShowers,
                                                                       const art::Event &event)
{
    //Using all the reco particles to determine the angle
    Point_t vertex(0,0,0); //No meaning vertex value is filled with this method

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
    Point_t vertex(0,0,0); //No meaning vertex value is filled with this method

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

dune::AngularRecoOutput NeutrinoAngularRecoAlg::CalculateNeutrinoAngle(const art::Event &event, const Point_t &vertex){
    std::vector<PolarFitOutput> output_vec;
    //Getting all the reconstructed 2D hits
    const std::vector<art::Ptr<recob::Hit>> hits(dune_ana::DUNEAnaEventUtils::GetHits(event, fHitLabel));
    //Sorting all the hits per view
    //Taking as reference the views for positive drift -> U- == V+ and V- == U+
    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit>>> hits_per_view;
    for(art::Ptr<recob::Hit> const& hit: hits){
        geo::View_t view = GetTargetView(hit);

        if(hits_per_view.count(view) == 0){
            hits_per_view[view] = {};
        }
        hits_per_view[view].push_back(hit);
    }

    for(auto const& [view, hits] : hits_per_view){
        PolarFitOutput output = FitViewHits(event, view, hits, vertex);
        output_vec.push_back(output);
    }

    const std::vector<std::pair<double, double>> start_values = { //Changng the starting point seems to help convergence in some cases
        {0.5, 0.5},
        {0.5, 2},
        {2, 0.5},
        {2, 2},
    };

    //Combining together the infos from the views
    CalorimetricDirectionFitter const caloFitter{output_vec, GetViewTheta(geo::kU), GetViewTheta(geo::kV)};
    ROOT::Math::Functor fitFunc([&caloFitter](double const* xs){return caloFitter.Chi2(xs);}, 2);

    if(caloFitter.NViews() < 2){
        Point_t vertex(0,0,0);
        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kHits);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }


    ROOT::Minuit2::Minuit2Minimizer minimizer{ROOT::Minuit2::kCombined};
    
    bool mstatus = false;

    for(const std::pair<double, double> &start_value : start_values){
        minimizer.Clear();
        minimizer.SetFunction(fitFunc);
        minimizer.SetVariable(0, "theta", start_value.first, 0.01);
        minimizer.SetVariable(1, "phi", start_value.second, 0.01);
        minimizer.SetMaxFunctionCalls(1e9);
        minimizer.SetMaxIterations(1e9);
        minimizer.SetTolerance(0.01);
        minimizer.SetStrategy(2);
        minimizer.SetErrorDef(1);
        mstatus = minimizer.Minimize();

        if(mstatus){
            break; //If fit converged, no need to test alternative stating values
        }
    }

    if(!mstatus){ //Did not manage to converge, all hope is lost
        Point_t vertex(0,0,0);
        Direction_t direction(0,0,0);
        AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kRecoMethodNotSet);
        return this->ReturnNeutrinoAngle(angularRecoInputHolder);
    }

    const double* pars = minimizer.X();
    double theta = pars[0];
    double phi = pars[1];
    Direction_t direction(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    AngularRecoInputHolder angularRecoInputHolder(vertex, direction, kHits);
    return this->ReturnNeutrinoAngle(angularRecoInputHolder);
}

//------------------------------------------------------------------------------------------------------------------------------------------

geo::View_t NeutrinoAngularRecoAlg::GetTargetView(const art::Ptr<recob::Hit> &hit) const{
    const geo::TPCID hit_tpcID(hit->WireID());
    if(fGeometry->TPC(hit_tpcID).DriftAxisWithSign().sign == geo::DriftSign::Positive){ //"normal" drift direction, we don't change anything
        return hit->View();
    }
    else{ //We invert U and V
        if(hit->View() == geo::kU){
            return geo::kV;
        }
        else if(hit->View() == geo::kV){
            return geo::kU;
        }
        else{
            return hit->View();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float NeutrinoAngularRecoAlg::GetViewTheta(geo::View_t view) const{
    for (geo::TPCGeo const& TPC: fGeometry->Iterate<geo::TPCGeo>()) {
        if(TPC.DriftAxisWithSign().sign == geo::DriftSign::Positive){ //Trying to find a TPC with the right drift direction
            for (unsigned int p = 0; p < fWireReadout->Get().Nplanes(TPC.ID()); ++p) {
                geo::PlaneGeo const& plane = fWireReadout->Get().Plane(geo::PlaneID(TPC.ID(), p));
                if (plane.View() == view){
                    return 0.5f * M_PI - plane.ThetaZ();
                }
            }
        }
    } // for TPCs

    throw cet::exception("NeutrinoAngularRecoAlg") << "Could not find a TPC with the right drift direction for view " << view << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PolarFitOutput NeutrinoAngularRecoAlg::FitViewHits(const art::Event &event, const geo::View_t &view, const std::vector<art::Ptr<recob::Hit>> &hits, const Point_t &vertex) const{
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    // const geo::TPCID tpc(0, 0);

    if(view != geo::View_t::kU && view != geo::View_t::kV && view != geo::View_t::kW){
        mf::LogWarning("NeutrinoAngularRecoAlg") << "Unknown view " << view << ". Skipping!";
        PolarFitOutput output = {};
        output.success = false;
        return output;
    }

    std::vector<float> x_vec;
    std::vector<float> y_vec;
    std::vector<float> adc_vec;
    
    //Getting the wire angle for this view\TPC wire
    const float theta = GetViewTheta(view);

    for(art::Ptr<recob::Hit> const &hit : hits){
        const geo::WireID hit_WireID(hit->WireID());

        //Getting X from the drift
        float hitX = detProp.ConvertTicksToX(hit->PeakTime(), hit_WireID);
        //Shifting X by the vertex position
        hitX -= vertex.X();

        //Getting hit"Y" from the view coordinate
        
        auto const xyz = fWireReadout->Get().Wire(hit_WireID).GetCenter();
        
        const float hitY((xyz.Z() - vertex.Z())*cos(theta) - (xyz.Y() - vertex.Y())*sin(theta));

        //Getting the hit ADC
        const float adc = hit->Integral()*dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, hit);

        x_vec.push_back(hitX);
        y_vec.push_back(hitY);
        adc_vec.push_back(adc);

    }

    PolarFitOutput output(PolarFit(x_vec, y_vec, adc_vec));
    output.view = view;
    output.nhits = hits.size();

    return output;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PolarFitOutput NeutrinoAngularRecoAlg::PolarFit(const std::vector<float> &xvec, const std::vector<float> &yvec, const std::vector<float> &weights) const{
    std::complex<float> sum = 0;
    PolarFitOutput output = {};

    const uint npoints = xvec.size();

    if(npoints < 1){
        //Impossible to make a fit
        output.success = false;
        return output;
    }

    for(uint i = 0; i < npoints; i++){
        std::complex<float> new_direction(xvec[i], yvec[i]); //Using complex numbers to represent 2D vector from vertex to point
        new_direction /= sqrtf(std::norm(new_direction)); //Normalizing the direction
        sum += new_direction*weights[i]; //Add this direction weighted by the correct weight
    }

    output.pdrift = sum.real();
    output.pview = sum.imag();
    
    output.success = true;

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

        if(dir.X() < -1){ //Direction is not defined (-999, -999, -999),
            continue; //We skip this shower
        }

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

        const recob::Track::Vector_t &dir = pTrack->StartDirection();

        if(dir.X() < -1){ //Direction is not defined (-999, -999, -999),
            continue; //We skip this shower
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


CalorimetricDirectionFitter::CalorimetricDirectionFitter(std::vector<PolarFitOutput>& observed_data, double thetaU, double thetaV) :
    _thetaU(thetaU), _thetaV(thetaV)
{

    for(const PolarFitOutput& view_fit : observed_data){
        if(!view_fit.success){
            continue;
        }

        switch (view_fit.view)
        {
            case geo::View_t::kU:
                _pU = view_fit.pview;
                _pxU = view_fit.pdrift;
                _nhitsU = view_fit.nhits;
                break;
            case geo::View_t::kV:
                _pV = view_fit.pview;
                _pxV = view_fit.pdrift;
                _nhitsV = view_fit.nhits;
                break;
            case geo::View_t::kW:
                _pW = view_fit.pview;
                _pxW = view_fit.pdrift;
                _nhitsW = view_fit.nhits;
                break;
            default:
                break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------


double CalorimetricDirectionFitter::Chi2(double const* x) const{
    double theta = x[0];
    double phi = x[1];

    double px = sin(theta)*cos(phi);
    double py = sin(theta)*sin(phi);
    double pz = cos(theta);

    double pu = cos(_thetaU)*pz - sin(_thetaU)*py;
    double pv = cos(_thetaV)*pz - sin(_thetaV)*py;

    double u_theta_meas = atan2(_pU, _pxU);
    double v_theta_meas = atan2(_pV, _pxV);
    double w_theta_meas = atan2(_pW, _pxW);

    double u_theta_fit = atan2(pu, px);
    double v_theta_fit = atan2(pv, px);
    double w_theta_fit = atan2(pz, px);

    double chi2 = 0;
    if (_nhitsU > 0){
        chi2 += (u_theta_fit - u_theta_meas)*(u_theta_fit - u_theta_meas);
    }
    if(_nhitsV > 0){
        chi2 += (v_theta_fit - v_theta_meas)*(v_theta_fit - v_theta_meas);
    }
    if(_nhitsW > 0)
    {
        chi2 += (w_theta_fit - w_theta_meas)*(w_theta_fit - w_theta_meas);
    }

    return chi2;

}

//------------------------------------------------------------------------------------------------------------------------------------------

int CalorimetricDirectionFitter::NViews() const{
    return (_nhitsU > 0) + (_nhitsV > 0) + (_nhitsW > 0);
}

} //dune
