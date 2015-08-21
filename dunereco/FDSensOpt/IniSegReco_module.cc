////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       IniSegAlg
// Module Type: producer
// File:        IniSegReco_module.cc
// Authors:      dorota.stefan@cern.ch, robert.sulej@cern.ch, tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Vertex.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/TimeService.h"
#include "MCCheater/BackTracker.h"
#include "SimulationBase/MCTruth.h"
#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"
#include "IniSegAlg/IniSegAlg.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

// C/C++ libraries
#include <memory>
#include <utility>

namespace dunefd {
	class IniSegReco;
	class Hit2D;
	class IniSegAlg;
}

class dunefd::IniSegReco : public art::EDProducer {
public:
  explicit IniSegReco(fhicl::ParameterSet const & p);

  IniSegReco(IniSegReco const &) = delete;
  IniSegReco(IniSegReco &&) = delete;
  IniSegReco & operator = (IniSegReco const &) = delete;
  IniSegReco & operator = (IniSegReco &&) = delete;

  void beginJob() override;

  void reconfigure(fhicl::ParameterSet const& p) override;

  void produce(art::Event & e) override;


private:
	void ResetVars(void);

	bool insideFidVol(TLorentzVector const & pvtx) const;
	bool insideFidVolandTPC(TLorentzVector const & pvtx, size_t tpc) const;

	TVector2 getMCvtx2d(TVector3 const & mcvtx3d, 
											const size_t cryo, const size_t tpc, const size_t plane) const;
	TVector2 getMCdir2d(TVector3 const & mcvtx3d, TVector3 const & mcdir3d, 
											const size_t cryo, const size_t tpc, const size_t plane) const;	
	TVector3 findElDir(art::Ptr<simb::MCTruth> const mctruth) const;
	std::vector< TVector3 > findDirs(art::Ptr<simb::MCTruth> const mctruth, int pdg) const;
 
	void collectCls(art::Event const & evt, art::Ptr<simb::MCTruth> const mctruth);
	std::vector< dunefd::Hit2D > reselectCls(std::map<size_t, std::vector< dunefd::Hit2D > > const & cls, 
																					art::Ptr<simb::MCTruth> const mctruth,
																					const size_t cryo, const size_t tpc, const size_t plane) const;
	void make3dseg(art::Event const & evt, std::vector< std::vector<Hit2D> > const & src, TVector3 const & primary);

	recob::Track convertFrom(pma::Track3D const & src);

	std::vector< pma::Track3D* > pmatracks;

  TTree *fTree;
  
  int run;
  int subrun;
  int event;
	short isdata;
	int trkindex;        // track index in the event

	Float_t lep_dedx;  
	Float_t lep_dist;
	int ngamma;
	Float_t convdist;

	std::string fHitsModuleLabel;
	std::string fClusterModuleLabel;
	std::string fGenieGenModuleLabel;

	double fFidVolCut;
	pma::ProjectionMatchingAlg fProjectionMatchingAlg;

};
// ------------------------------------------------------

dunefd::IniSegReco::IniSegReco(fhicl::ParameterSet const & pset)
: fProjectionMatchingAlg(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"))
{
	this->reconfigure(pset);
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< art::Assns<recob::Track, recob::Hit> >();
	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
}
// ------------------------------------------------------

void dunefd::IniSegReco::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("nueana","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
	fTree->Branch("lep_dedx",&lep_dedx,"lep_dedx/F");
	fTree->Branch("lep_dist",&lep_dist,"lep_dist/F");
	fTree->Branch("ngamma", &ngamma, "ngamma/I");
	fTree->Branch("convdist", &convdist, "convdist/F");
}

void dunefd::IniSegReco::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitsModuleLabel     =   pset.get< std::string >("HitsModuleLabel");
	fClusterModuleLabel  =   pset.get< std::string >("ClusterModuleLabel");
	fGenieGenModuleLabel =   pset.get< std::string >("GenieGenModuleLabel");

	fProjectionMatchingAlg.reconfigure(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"));
	fFidVolCut           =   pset.get< double >("FidVolCut");
	return;
}

void dunefd::IniSegReco::ResetVars()
{
	ngamma = 0;
	convdist = -10.0F;
	pmatracks.clear();
	return;
}

void dunefd::IniSegReco::produce(art::Event& evt)
{
	ResetVars();
	art::ServiceHandle<geo::Geometry> geom;
	run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
	isdata = evt.isRealData();

	std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
	std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);

	std::unique_ptr< art::Assns< recob::Track, recob::Hit > > trk2hit(new art::Assns< recob::Track, recob::Hit >);
	std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
	std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);

	TVector3 primary(0, 0, 0);
	if (!isdata)
	{
		 // * MC truth information
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    if (mclist.size()){
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
			collectCls(evt, mctruth);
			const TLorentzVector& pvtx = mctruth->GetNeutrino().Nu().Position();
			primary = TVector3(pvtx.X(), pvtx.Y(), pvtx.Z());

			art::ServiceHandle<cheat::BackTracker> bt;
			const sim::ParticleList& plist = bt->ParticleList();

			bool photon = false;
			for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
			{
				simb::MCParticle* particle = ipar->second;
				TLorentzVector mom = particle->Momentum();
				TVector3 momvec(mom.Px(), mom.Py(), mom.Pz());
			
				if ((particle->PdgCode() == 22) && (momvec.Mag() > 0.030))
				{		
					ngamma++; photon = true;
					TLorentzVector conversion = particle->EndPosition();
					TVector3 convec(conversion.X(), conversion.Y(), conversion.Z());
					convdist = std::sqrt(pma::Dist2(primary, convec));				

					for (int i = 0; i < particle->NumberDaughters(); i++)
					{		
						std::cout << " pdg = " << bt->TrackIDToParticle(particle->Daughter(i))->PdgCode() << std::endl;
						std::cout << " no of particles " << particle->NumberDaughters() << std::endl;
					}
				}
			}
			if (!photon) ngamma = -10;

			std::cout << " ****** end ****** " << std::endl;
		}
	}

	if (pmatracks.size())
	{
		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6];
		for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

		trkindex = 0;
		for (auto trk : pmatracks)
		{
			lep_dedx = fProjectionMatchingAlg.selectInitialHits(*trk, geo::kZ);
			if (lep_dedx == 0) lep_dedx = fProjectionMatchingAlg.selectInitialHits(*trk, geo::kV);
			if (lep_dedx == 0) lep_dedx = fProjectionMatchingAlg.selectInitialHits(*trk, geo::kU);

			lep_dist = std::sqrt(pma::Dist2(primary, trk->front()->Point3D()));

		
			fTree->Fill();

			tracks->push_back(convertFrom(*trk));
			trkindex++;

			std::vector< art::Ptr< recob::Hit > > hits2d;
			art::PtrVector< recob::Hit > sp_hits;
			spStart = allsp->size();
			for (int h = trk->size() - 1; h >= 0; h--)
			{
				pma::Hit3D* h3d = (*trk)[h];
				hits2d.push_back(h3d->Hit2DPtr());

				if ((h == 0) ||
					      (sp_pos[0] != h3d->Point3D().X()) ||
					      (sp_pos[1] != h3d->Point3D().Y()) ||
					      (sp_pos[2] != h3d->Point3D().Z()))
				{
					if (sp_hits.size()) // hits assigned to the previous sp
					{
						util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
						sp_hits.clear();
					}
					sp_pos[0] = h3d->Point3D().X();
					sp_pos[1] = h3d->Point3D().Y();
					sp_pos[2] = h3d->Point3D().Z();
					allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
				}
					sp_hits.push_back(h3d->Hit2DPtr());
			}
			if (sp_hits.size()) // hits assigned to the last sp
			{
				util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
			}
			spEnd = allsp->size();

			if (hits2d.size())
			{
					util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);
					util::CreateAssn(*this, evt, *tracks, hits2d, *trk2hit);
			}
		}
		// data prods done, delete all pma::Track3D's
		for (size_t t = 0; t < pmatracks.size(); ++t) delete pmatracks[t];

	}
	

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(sp2hit));
}


/***********************************************************************/

void dunefd::IniSegReco::collectCls(art::Event const & evt, art::Ptr<simb::MCTruth> const mctruth)
{
	art::ServiceHandle<geo::Geometry> geom;
	
	// * clusters
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
 	std::vector<art::Ptr<recob::Cluster> > clusterlist;

	if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
	{
   	art::fill_ptr_vector(clusterlist, clusterListHandle);
		art::FindManyP< recob::Hit > hc(clusterListHandle, evt, fClusterModuleLabel);		

		for (size_t c = 0; c < geom->Ncryostats(); ++c) // iter...
		{
			const geo::CryostatGeo& cryo = geom->Cryostat(c);
			for (size_t t = 0; t < cryo.NTPC(); ++t)
			{
				const geo::TPCGeo& tpc = cryo.TPC(t);
				std::vector< std::vector< dunefd::Hit2D > > clinput;
				for (size_t p = 0; p < tpc.Nplanes(); ++p) 
				{	
					std::map<size_t, std::vector< dunefd::Hit2D > > cls; 
					for (size_t i = 0; i < clusterlist.size(); ++i)
					{
						std::vector< art::Ptr<recob::Hit> > hitscl;
						std::vector< Hit2D > hits2dcl;
						hitscl = hc.at(i);
					
						bool found = false;
						for (size_t h = 0; h < hitscl.size(); ++h)
						{
							size_t wire = hitscl[h]->WireID().Wire;
							size_t plane = hitscl[h]->WireID().Plane;
							size_t tpc = hitscl[h]->WireID().TPC;
							size_t cryo = hitscl[h]->WireID().Cryostat;

							if ((plane == p) && (tpc == t) && (cryo == c))
							{
								found = true;
								TVector2 point = pma::WireDriftToCm(wire, hitscl[h]->PeakTime(), plane, tpc, cryo);
								Hit2D hit2d(point, hitscl[h].key());
								hits2dcl.push_back(hit2d);
							}
						}

						if (found) cls[i] = hits2dcl;
					}

					if (cls.size())
							clinput.push_back(reselectCls(cls, mctruth, c, t, p));							
				}

				if (clinput.size() > 1)
				{
					const TLorentzVector& pvtx = mctruth->GetNeutrino().Nu().Position();
					TVector3 primary(pvtx.X(), pvtx.Y(), pvtx.Z());
					make3dseg(evt, clinput, primary);
				}
			}
		}

	}		
}

/***********************************************************************/

std::vector< dunefd::Hit2D > dunefd::IniSegReco::reselectCls(std::map<size_t, std::vector< dunefd::Hit2D > > const & cls, art::Ptr<simb::MCTruth> const mctruth,
																														const size_t cryo, const size_t tpc, const size_t plane) const
{
	std::vector< dunefd::Hit2D > cluster;
	if (!cls.size()) return cluster;

	const simb::MCParticle& particle = mctruth->GetNeutrino().Nu();
	const TLorentzVector& pvtx = particle.Position();
	
	if (insideFidVolandTPC(pvtx, tpc) && (abs(mctruth->GetNeutrino().Lepton().PdgCode()) == 11))
	{
		// mc
		TVector3 mcvtx3d(pvtx.X(), pvtx.Y(), pvtx.Z());
		TVector3 mcdir3d = findElDir(mctruth);
		TVector2 mcvtx2d = getMCvtx2d(mcvtx3d, cryo, tpc, plane);
		TVector2 test = pma::CmToWireDrift(mcvtx2d.X(), mcvtx2d.Y(), plane, tpc, cryo);
		TVector2 mcdir2d = getMCdir2d(mcvtx3d, mcdir3d, cryo, tpc, plane);
		
		// reco: find the best clusters to proceed with segment reconstruction
		IniSegAlg recoini(cls);
		recoini.FeedwithMc(mcvtx2d, mcdir2d);
		cluster = recoini.GetCl();
	}

	return cluster;
}

/***********************************************************************/

void dunefd::IniSegReco::make3dseg(art::Event const & evt, 
		std::vector< std::vector< dunefd::Hit2D > > const & src, 
		TVector3 const & primary)
{	
	if (src.size() < 2) return;
	if ((src[0].size() < 2) || (src[1].size() < 2)) return;

	art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

	std::vector< art::Ptr<recob::Hit> > hitsel;

	for (size_t i = 0; i < src.size(); ++i)
		for (size_t h = 0; h < src[i].size(); ++h)
			hitsel.push_back(hitlist[src[i][h].GetKey()]);

	if (hitsel.size() > 5)
	{
		pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(hitsel, primary);
		pmatracks.push_back(trk);
	}
}

/***********************************************************************/

TVector2 dunefd::IniSegReco::getMCdir2d(TVector3 const & mcvtx3d, TVector3 const & mcdir3d, 
																				const size_t cryo, const size_t tpc, const size_t plane) const
{
	TVector3 shift3d = mcvtx3d + mcdir3d;
	TVector2 shift2d = pma::GetProjectionToPlane(shift3d, plane, tpc, cryo) - getMCvtx2d(mcvtx3d, cryo, tpc, plane);
	TVector2 mcdir2d = shift2d * (1.0 / shift2d.Mod());

	return mcdir2d;
}

/***********************************************************************/

TVector2 dunefd::IniSegReco::getMCvtx2d(TVector3 const & mcvtx3d, 
																				const size_t cryo, const size_t tpc, const size_t plane) const
{
	TVector2 mcvtx2d = pma::GetProjectionToPlane(mcvtx3d, plane, tpc, cryo);

	return mcvtx2d;
}

/***********************************************************************/

TVector3 dunefd::IniSegReco::findElDir(art::Ptr<simb::MCTruth> const mctruth) const
{
	TVector3 dir(0, 0, 0);
	const simb::MCParticle& lepton = mctruth->GetNeutrino().Lepton();

	if (abs(lepton.PdgCode()) == 11)
	{
		TLorentzVector mom = lepton.Momentum();
		TVector3 momvec3(mom.Px(), mom.Py(), mom.Pz());
		dir = momvec3 * (1 / momvec3.Mag());
	}

	return dir;
}

/***********************************************************************/

std::vector< TVector3 > dunefd::IniSegReco::findDirs(art::Ptr<simb::MCTruth> const mctruth, int pdg) const
{
	std::vector< TVector3 > dirs;

	art::ServiceHandle<cheat::BackTracker> bt;
	const sim::ParticleList& plist = bt->ParticleList();

	for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
	{
		simb::MCParticle* particle = ipar->second;
		TLorentzVector mom = particle->Momentum();
		TVector3 momvec(mom.Px(), mom.Py(), mom.Pz());
			
		if ((particle->PdgCode() == pdg) && (momvec.Mag() > 0.030))
		{		
			TLorentzVector momconv = particle->EndMomentum();
			TVector3 momconvec3(momconv.Px(), momconv.Py(), momconv.Pz());			
			TVector3 dir = momconvec3 * (1 / momconvec3.Mag());
			dirs.push_back(dir);
		}
	}
	
	return dirs;
}

/***********************************************************************/

bool dunefd::IniSegReco::insideFidVol(TLorentzVector const & pvtx) const
{
	art::ServiceHandle< geo::Geometry > geom;

	double vtx[3];
	vtx[0] = pvtx.X(); vtx[1] = pvtx.Y(); vtx[2] = pvtx.Z();

	bool inside = false;
	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);
	if (geom->HasTPC(idtpc))
	{
		const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc); 
		if (((vtx[0] - tpcgeo.MinX()) > fFidVolCut) &&
				((tpcgeo.MaxX() - vtx[0]) > fFidVolCut) &&
				((vtx[1] - tpcgeo.MinY()) > fFidVolCut) &&
				((tpcgeo.MaxY() - vtx[1]) > fFidVolCut) &&
				((vtx[2] - tpcgeo.MinZ()) > fFidVolCut) &&
				((tpcgeo.MaxZ() - vtx[2]) > fFidVolCut)) inside = true;
	}
		
	return inside;
}

/***********************************************************************/

bool dunefd::IniSegReco::insideFidVolandTPC(TLorentzVector const & pvtx, size_t tpc) const
{
	art::ServiceHandle< geo::Geometry > geom;

	double vtx[3];
	vtx[0] = pvtx.X(); vtx[1] = pvtx.Y(); vtx[2] = pvtx.Z();

	bool inside = false;
	geo::TPCID const & idtpc = geom->FindTPCAtPosition(vtx);
	if (geom->HasTPC(idtpc) && (idtpc.TPC == tpc))
	{
		const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc); 
		if (((vtx[0] - tpcgeo.MinX()) > fFidVolCut) &&
				((tpcgeo.MaxX() - vtx[0]) > fFidVolCut) &&
				((vtx[1] - tpcgeo.MinY()) > fFidVolCut) &&
				((tpcgeo.MaxY() - vtx[1]) > fFidVolCut) &&
				((vtx[2] - tpcgeo.MinZ()) > fFidVolCut) &&
				((tpcgeo.MaxZ() - vtx[2]) > fFidVolCut)) inside = true;
	}
		
	return inside;
}

/***********************************************************************/

recob::Track dunefd::IniSegReco::convertFrom(pma::Track3D const & src)
{
	std::vector< TVector3 > xyz, dircos;

	for (size_t i = 0; i < src.size(); ++i)
	{
		xyz.push_back(src[i]->Point3D());

		if (i < src.size() - 1)
		{
			TVector3 dc(src[i + 1]->Point3D());
			dc -= src[i]->Point3D();
			dc *= 1.0 / dc.Mag();
			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());
	}

	if (xyz.size() != dircos.size())
	{
		mf::LogError("IniSegReco") << "pma::Track3D to recob::Track conversion problem.";
	}
	return recob::Track(xyz, dircos, std::vector< std::vector< double > >(0), std::vector< double >(2, util::kBogusD), trkindex);
}

/***********************************************************************/

DEFINE_ART_MODULE(dunefd::IniSegReco)
