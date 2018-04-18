////////////////////////////////////////////////////////////////////////
// Classes:     IniSegAlg, Hit2D, bDistCentLess2D
// File:        IniSegAlg.h
//
// dorota.stefan@cern.ch, robert.sulej@cern.ch, tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef IniSegAlg_h
#define IniSegAlg_h

#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "lardataobj/RecoBase/Track.h"

namespace fhicl {
	class ParameterSet;
}

namespace dunefd {
	class Hit2D;
	class IniSegAlg;
	class bDistCentLess2D;
}

class dunefd::Hit2D
{
	public:
	Hit2D(TVector2 point2d, size_t key);
	TVector2 const & GetPointCm(void) const { return fPoint;}
	size_t const & GetKey(void) const { return fKey; }

	private:
	TVector2 fPoint;
	size_t fKey;
};

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif
class dunefd::IniSegAlg
{
	public:
	IniSegAlg(std::map<size_t, std::vector<dunefd::Hit2D> > clusters); 
	IniSegAlg(std::vector< art::Ptr<recob::Track> > const & tracks, TVector3 const & mcvtx); 

	void FeedwithMc(TVector2 const & vtx, TVector2 const & dir, TVector3 const & dir3d);
	void FeedwithMc(TVector3 const & dir3d);
	std::map<size_t, std::vector< dunefd::Hit2D > > const & GetSelectedCl() const { return fSelCls; } 
	std::vector< dunefd::Hit2D > const & GetCl() const { return fCl; }

	const bool IsFound() { return fFound; }
	art::Ptr<recob::Track> const & GetTrk() const { return fTrk; }

	float const & GetDist() const { return fDistVtxCl; } 
	float const & GetCos() const { return fCos; }

	private:
	void FindClustersInRad(); 
	void FindCluster();
	void SortLess(); 

	void Find3dTrack();

	TVector2 ClusterDir(std::vector< Hit2D > const & hits);

	std::map<size_t, std::vector<dunefd::Hit2D> > fClusters;
	std::map<size_t, std::vector<dunefd::Hit2D> > fSelCls;
	//

	std::vector< art::Ptr<recob::Track> > fSelTrks;
	art::Ptr<recob::Track> fTrk;
	bool fFound;
	TVector3 fDir3d;

	TVector3 fFront;
	TVector3 fBack;
	
	//
	std::vector< dunefd::Hit2D > fCl;
	TVector2 fMcVtx;
	TVector3 fMcVtx3d;
	TVector2 fDir;

	float fRadius;
	float fDistVtxCl;
	float const fThrcos;
	float fCos;
};
#if defined __clang__
  #pragma clang diagnostic pop
#endif

class dunefd::bDistCentLess2D :
	public std::binary_function< Hit2D, Hit2D, bool >
	{
		public:
		bDistCentLess2D(const TVector2& c) : center(c) {}

		bool operator() (Hit2D p1, Hit2D p2)
		{
			double dist1 = pma::Dist2(p1.GetPointCm(), center);
			double dist2 = pma::Dist2(p2.GetPointCm(), center);

			return dist1 < dist2;
		}

		private:
		TVector2 center; 
	};

#endif //IniSegAlg_h
