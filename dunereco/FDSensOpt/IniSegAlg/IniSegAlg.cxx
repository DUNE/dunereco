////////////////////////////////////////////////////////////////////////
// Classes:     IniSegAlg, Hit2D, bDistCentLess2D
// File:        IniSegAlg.cxx
//
// dorota.stefan@cern.ch, robert.sulej@cern.ch, tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include "IniSegAlg.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT libraries
#include "TMath.h"

// C/C++ libraries
#include <memory>

dunefd::Hit2D::Hit2D(TVector2 point2d, size_t key) :
fPoint(point2d),
fKey(key)
{
}

dunefd::IniSegAlg::IniSegAlg(std::map<size_t, std::vector<dunefd::Hit2D> > clusters) :
fClusters(clusters),
fRadius(10.0),
fDistVtxCl(0.0F),
fThrcos(0.9),
fCos(0.0F)
{
}

dunefd::IniSegAlg::IniSegAlg(std::vector< art::Ptr<recob::Track> > const & tracks, TVector3 const & mcvtx) :
fRadius(10.0),
fDistVtxCl(0.0F),
fThrcos(0.9),
fCos(-9999.0F)
{
	fSelTrks = tracks;
	fMcVtx3d = mcvtx;
	fFound = false;
}

void dunefd::IniSegAlg::FeedwithMc(TVector2 const & vtx, TVector2 const & dir, TVector3 const & dir3d)
{
	fMcVtx = vtx;
	fDir = dir;
	fDir3d = dir3d;
	FindClustersInRad();
	SortLess();
	FindCluster();
}

void dunefd::IniSegAlg::FeedwithMc(TVector3 const & dir3d)
{
	fDir3d = dir3d;
	Find3dTrack();
}

void dunefd::IniSegAlg::FindClustersInRad()
{
	for (auto& cl: fClusters)
		for (size_t h = 0; h < cl.second.size(); ++h)
			if (pma::Dist2(fMcVtx, cl.second[h].GetPointCm()) < fRadius*fRadius)
			{
				fSelCls[cl.first] = cl.second;
				break;
			}
}

void dunefd::IniSegAlg::SortLess()
{
	for (auto& cl: fSelCls)
		std::sort(cl.second.begin(), cl.second.end(), bDistCentLess2D(fMcVtx));
}

TVector2 dunefd::IniSegAlg::ClusterDir(std::vector< Hit2D > const & hits)
{
	TVector2 p1 = hits[0].GetPointCm();
	size_t id = 0;
	for (size_t i = 0; i < hits.size(); ++i)
		if (pma::Dist2(hits[0].GetPointCm(), hits[i].GetPointCm()) < 1.0)
		{
			id = i; 
		}
		else break;

	TVector2 p2 = hits[id].GetPointCm(); 
	TVector2 dir = (p2 - p1) * (1 / (p2 - p1).Mod());

	return dir;	
}

void dunefd::IniSegAlg::FindCluster()
{
	double maxcos = 0;
	std::map<size_t, std::vector<dunefd::Hit2D> > chosen;
	for (auto& cl: fSelCls)
		if (cl.second.size() > 2)
		{
			TVector2 dir = ClusterDir(cl.second);
			double cos = fDir * dir;
			if (cos > maxcos) 
			{
				maxcos = cos;
				chosen[cl.first] = cl.second;
				fCl.clear();
				fCl = cl.second;
				fDistVtxCl = std::sqrt(pma::Dist2(fCl[0].GetPointCm(), fMcVtx));
			}
		}

	fSelCls.clear();
	fSelCls = chosen;
}

void dunefd::IniSegAlg::Find3dTrack()
{
	double larStart[3] = {0.0, 0.0, 0.0};
  	double larEnd[3] = {0.0, 0.0, 0.0};

	double maxcos = 0.0; 
	const double thrdist = 10; // cm

	for (auto& trk: fSelTrks)
	{
		trk->Direction(larStart,larEnd); // correct, check both directions
		TVector3 dir(larStart[0], larStart[1], larStart[2]);
	
		if (!trk->NumberTrajectoryPoints()) continue;
		TVector3 pos_p = trk->LocationAtPoint(0);
		TVector3 pos_end = trk->LocationAtPoint(trk->NumberTrajectoryPoints()-1);

		double dist = std::sqrt(pma::Dist2(pos_p, fMcVtx3d));		

		if (dist > std::sqrt(pma::Dist2(pos_end, fMcVtx3d))) continue;

		double cos = fDir3d * dir;
		
		if ((cos > maxcos) && (dist < thrdist))
		{		
			maxcos = cos;
			fTrk = trk;
			fCos = cos;
			fFound = true;
		}
	}
}


