///////////////////////////////////////////////////////////////////////////////////////
//
// EmLikeHits class
// 
// robert.sulej@cern.ch
//  
// This module produces new hit container with hadron/muon-like hits subtracted. The
// goal is to provide EM cascade-like hits as input to clustering algorithms that
// expect EM showers not tracks.
// 
// Module needs reconstructed and tagged tracks on its input. For the moment tagging
// is encoded in the track ID - TO BE CHANGED TO A SIMPLE DATA PRODUCT ASSIGNED TO
// TRACK.
// 
// Also hits unmatched to any track are checked if they are not very close to any of
// hadron/muon tracks. The distance is calculated in 2D projection, using track
// trajectory points. As soon as pma::Tracks are able to produce assiciated data
// product with node description the distance will use NODES INSTEAD OF TRAJECTORY
// POINTS since many trajectory points may belong to single linear segment.
//
///////////////////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

namespace dune {

class EmLikeHits : public art::EDProducer {

public:

  explicit EmLikeHits(fhicl::ParameterSet const & p);

  EmLikeHits(EmLikeHits const &) = delete;

  EmLikeHits(EmLikeHits &&) = delete;

  EmLikeHits & operator = (EmLikeHits const &) = delete;

  EmLikeHits & operator = (EmLikeHits &&) = delete;

  void reconfigure(fhicl::ParameterSet const& p);

  void produce(art::Event & e) override;

private:

  void removeHitsAssignedToTracks(
	std::vector< art::Ptr<recob::Hit> >& hitlist,
	const std::vector<recob::Track>& tracks,
	const art::FindManyP< recob::Hit >& fbp);

  void removeUnmatchedHitsCloseToTracks(
        detinfo::DetectorPropertiesData const& detProp,
	std::vector< art::Ptr<recob::Hit> >& hitlist,
	const std::vector<recob::Track>& tracks,
	const art::FindManyP< recob::Hit >& fbp);

  bool isCloseToTrack(
	TVector2 p, const recob::Track& trk,
	unsigned int view,
	unsigned int tpc,
	unsigned int cryo);

  static double getDist2(
	const TVector2& psrc,
	const TVector2& p0,
	const TVector2& p1);

// ******************** parameters **********************

  std::string fHitModuleLabel;
  std::string fTrk3DModuleLabel;

};
// ------------------------------------------------------

EmLikeHits::EmLikeHits(fhicl::ParameterSet const & p) : EDProducer{p}
{
        this->reconfigure(p);
        produces< std::vector<recob::Hit> >();
}
// ------------------------------------------------------

void EmLikeHits::reconfigure(fhicl::ParameterSet const& pset)
{
        fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
        fTrk3DModuleLabel = pset.get< std::string >("Trk3DModuleLabel");
}
// ------------------------------------------------------

void EmLikeHits::removeHitsAssignedToTracks(
	std::vector< art::Ptr<recob::Hit> >& hitlist,
	const std::vector<recob::Track>& tracks,
	const art::FindManyP< recob::Hit >& fbp)
{
	for (size_t t = 0; t < tracks.size(); t++)
		if (!(tracks[t].ID() & 0x10000))
	{
		std::vector< art::Ptr<recob::Hit> > v = fbp.at(t);
		mf::LogVerbatim("EmLikeHits") << "   track-like trajectory: " << v.size() << std::endl;
		size_t ih = 0;
		while (ih < hitlist.size())
		{
			bool found = false;
			for (size_t it = 0; it < v.size(); it++)
				if (v[it].key() == hitlist[ih].key())
			{
				found = true; break;
			}
			if (found) hitlist.erase(hitlist.begin() + ih);
			else ih++;
		}
	}
}
// ------------------------------------------------------

bool EmLikeHits::isCloseToTrack(TVector2 p, const recob::Track& trk,
	unsigned int view, unsigned int tpc, unsigned int cryo)
{
	art::ServiceHandle<geo::Geometry> geom;
        geo::PlaneID const planeID{cryo, tpc, view};
        double wirePitch = geom->Plane(planeID).WirePitch();

	//double driftPitch = detProp.GetXTicksCoefficient(tpc, cryo);

	double max_d2_d = 0.3 * 0.3;
	double max_d2_w = (wirePitch + 0.1) * (wirePitch + 0.1);

	bool isClose = false;
	for (size_t i = 0; i < trk.NumberTrajectoryPoints() - 1; ++i)
	{
		TVector2 p0 = pma::GetVectorProjectionToPlane(trk.LocationAtPoint<TVector3>(i), view, tpc, cryo);
		TVector2 p1 = pma::GetVectorProjectionToPlane(trk.LocationAtPoint<TVector3>(i + 1), view, tpc, cryo);
		double d2 = getDist2(p, p0, p1);

		double dpx = fabs(p0.X() - p1.X());
		double dpy = fabs(p0.Y() - p1.Y());

		if (((dpx > 0.5 * dpy) && (d2 < max_d2_d)) || (d2 < max_d2_w))
		{
			isClose = true; break;
		}
	}
	return isClose;
}
// ------------------------------------------------------

void EmLikeHits::removeUnmatchedHitsCloseToTracks(
        detinfo::DetectorPropertiesData const& detProp,
	std::vector< art::Ptr<recob::Hit> >& hitlist,
	const std::vector<recob::Track>& tracks,
	const art::FindManyP< recob::Hit >& fbp)
{
	size_t ih = 0;
	while (ih < hitlist.size())
	{
		bool unmatched = true, close = false;;
		for (size_t t = 0; t < tracks.size(); t++)
		{
			std::vector< art::Ptr<recob::Hit> > v = fbp.at(t);
			for (size_t it = 0; it < v.size(); it++)
				if (v[it].key() == hitlist[ih].key())
				{
					unmatched = false; break;
				}
			if (!unmatched) break;
		}
		if (unmatched)
		{
			unsigned int plane = hitlist[ih]->WireID().Plane;
			unsigned int tpc = hitlist[ih]->WireID().TPC;
			unsigned int cryo = hitlist[ih]->WireID().Cryostat;

                        TVector2 hcm = pma::WireDriftToCm(detProp,
				hitlist[ih]->WireID().Wire, hitlist[ih]->PeakTime(), plane, tpc, cryo);

			for (size_t t = 0; t < tracks.size(); t++)
				if (!(tracks[t].ID() & 0x10000) &&
				    isCloseToTrack(hcm, tracks[t], plane, tpc, cryo))
			{
				close = true; break;
			}
		}
		if (close) hitlist.erase(hitlist.begin() + ih);
		else ih++;
	}
}
// ------------------------------------------------------

void EmLikeHits::produce(art::Event& evt)
{
	std::unique_ptr< std::vector< recob::Hit > > not_track_hits(new std::vector< recob::Hit >);

	std::vector< art::Ptr<recob::Hit> > hitlist;
	auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	auto trkListHandle = evt.getHandle< std::vector<recob::Track> >(fTrk3DModuleLabel);

	if (hitListHandle && trkListHandle)
	{
		art::FindManyP< recob::Hit > fbp(trkListHandle, evt, fTrk3DModuleLabel);

		art::fill_ptr_vector(hitlist, hitListHandle);
		mf::LogVerbatim("EmLikeHits") << "all hits: " << hitlist.size() << std::endl;

		removeHitsAssignedToTracks(hitlist, *trkListHandle, fbp);
                auto const detProp =
                  art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
                removeUnmatchedHitsCloseToTracks(detProp, hitlist, *trkListHandle, fbp);

		for (auto const& hit : hitlist) not_track_hits->push_back(recob::Hit(*hit));
		mf::LogVerbatim("EmLikeHits") << "remaining not track-like hits: " << not_track_hits->size() << std::endl;
	}
	evt.put(std::move(not_track_hits));
}
// ------------------------------------------------------

double EmLikeHits::getDist2(const TVector2& psrc, const TVector2& p0, const TVector2& p1)
{
	TVector2 v0(psrc); v0 -= p0;
	TVector2 v1(p1);   v1 -= p0;

	TVector2 v2(psrc); v2 -= p1;
	TVector2 v3(v1);   v3 *= -1.0;

	double v0Norm2 = v0.Mod2();
	double v1Norm2 = v1.Mod2();

	double eps = 1.0E-7; // 0.001mm
	if (v1Norm2 > eps)
	{
		double mag01 = sqrt(v0Norm2 * v1Norm2);
		double cosine01 = 0.0;
		if (mag01 != 0.0) cosine01 = v0 * v1 / mag01;

		double v2Norm2 = v2.Mod2();
		double mag23 = sqrt(v2Norm2 * v3.Mod2());
		double cosine23 = 0.0;
		if (mag23 != 0.0) cosine23 = v2 * v3 / mag23;

		double result = 0.0;
		if ((cosine01 > 0.0) && (cosine23 > 0.0))
		{
			result = (1.0 - cosine01 * cosine01) * v0Norm2;
		}
		else
		{
			if (cosine01 <= 0.0) result = v0Norm2;
			else result = v2Norm2;
		}

		if (result >= 0.0) return result;
		else return 0.0;
	}
	else
	{
		v1 = p0; v1 += p1; v1 *= 0.5;
		return pma::Dist2(v1, psrc);
	}
}
// ------------------------------------------------------

DEFINE_ART_MODULE(EmLikeHits)

} // namespace dune
