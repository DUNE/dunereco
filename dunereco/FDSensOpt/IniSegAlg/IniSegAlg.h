////////////////////////////////////////////////////////////////////////
// Classes:     IniSegAlg, Hit2D, bDistCentLess2D
// File:        IniSegAlg.h
//
// dorota.stefan@cern.ch, robert.sulej@cern.ch, tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef IniSegAlg_h
#define IniSegAlg_h

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

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

class dunefd::IniSegAlg
{
	public:
	IniSegAlg(std::map<size_t, std::vector<dunefd::Hit2D> > clusters); 

	void FeedwithMc(TVector2 const & vtx, TVector2 const & dir);
	std::map<size_t, std::vector< dunefd::Hit2D > > const & GetSelectedCl() const { return fSelCls; } 
	std::vector< dunefd::Hit2D > const & GetCl() const { return fCl; }
	float const & GetDist() const { return fDistVtxCl; } 

	private:
	void FindClustersInRad(); 
	void FindCluster();
	void SortLess(); 
	TVector2 ClusterDir(std::vector< Hit2D > const & hits);

	std::map<size_t, std::vector<dunefd::Hit2D> > fClusters;
	std::map<size_t, std::vector<dunefd::Hit2D> > fSelCls;
	std::vector< dunefd::Hit2D > fCl;
	TVector2 fMcVtx;
	TVector2 fDir;
	float fRadius;
	float fDistVtxCl;
};

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
