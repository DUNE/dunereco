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
#include "Utilities/AssociationUtil.h"
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
	
	float t0Corr(art::Event const & evt, TLorentzVector const & pvtx);

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

	void UseClusters(art::Event const & evt);

	void UseTracks(art::Event const & evt);

	std::vector< pma::Track3D* > pmatracks;

  	TTree *fTree;
  
  	int run;
  	int subrun;
  	int event;
	short isdata;
	int trkindex;        // track index in the event

	Float_t lep_dedx; 
	Float_t dx; 
	Float_t lep_dist;
	Float_t cos;
	Float_t t0;
	Float_t vtxrecomc;
	Float_t vtxrecomcx;
	Float_t vtxrecomcy;
	Float_t vtxrecomcz;
	int ngamma;
	Float_t convdist;

	std::string fHitsModuleLabel;
	std::string fClusterModuleLabel;
	std::string fTrackModuleLabel;
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
	fTree->Branch("dx",&dx,"dx/F");
	fTree->Branch("lep_dist",&lep_dist,"lep_dist/F");
	fTree->Branch("cos", &cos, "cos/F");
	fTree->Branch("t0", &t0, "t0/F");
	fTree->Branch("vtxrecomc", &vtxrecomc, "vtxrecomc/F");
	fTree->Branch("vtxrecomcx", &vtxrecomcx, "vtxrecomcx/F");
	fTree->Branch("vtxrecomcy", &vtxrecomcy, "vtxrecomcy/F");
	fTree->Branch("vtxrecomcz", &vtxrecomcz, "vtxrecomcz/F");
	fTree->Branch("ngamma", &ngamma, "ngamma/I");
	fTree->Branch("convdist", &convdist, "convdist/F");
}

void dunefd::IniSegReco::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitsModuleLabel     =   pset.get< std::string >("HitsModuleLabel");
	fClusterModuleLabel  =   pset.get< std::string >("ClusterModuleLabel");
	fTrackModuleLabel	 = 	pset.get< std::string >("TrackModuleLabel");
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
	lep_dist = -9999; //cm
	lep_dedx = -9999;
	t0 = -9999;
	dx = -9999;
	cos = -9999;
	vtxrecomc = -9999;
	vtxrecomcx = -9999;
	vtxrecomcy = -9999;
	vtxrecomcz = -9999;
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

		if (mclist.size())
		{
			art::Ptr<simb::MCTruth> mctruth = mclist[0];
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
				}
			}
			if (!photon) ngamma = -10;
		}
	}

	UseTracks(evt);

	fTree->Fill();

	if (!isdata && pmatracks.size())
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

void dunefd::IniSegReco::UseClusters(art::Event const & evt)
{
	if (!isdata)
	{
		 // * MC truth information
    		art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    		std::vector<art::Ptr<simb::MCTruth> > mclist;
    		if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      		art::fill_ptr_vector(mclist, mctruthListHandle);

    		if (mclist.size())
		{
      		art::Ptr<simb::MCTruth> mctruth = mclist[0];
			collectCls(evt, mctruth);
		}
	}
}

/***********************************************************************/

void dunefd::IniSegReco::UseTracks(art::Event const & evt)
{
	if (!isdata)
	{
		// * hits
		art::Handle< std::vector<recob::Hit> > hitListHandle;
  		std::vector<art::Ptr<recob::Hit> > hitlist;
 		if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    			art::fill_ptr_vector(hitlist, hitListHandle);

		// * tracks
		art::Handle< std::vector<recob::Track> > trackListHandle;
  		std::vector<art::Ptr<recob::Track> > tracklist;
  		if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    			art::fill_ptr_vector(tracklist, trackListHandle); 

		// * monte carlo
    		art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    		std::vector<art::Ptr<simb::MCTruth> > mclist;
    		if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      		art::fill_ptr_vector(mclist, mctruthListHandle);

		if (mclist.size())
		{
      		art::Ptr<simb::MCTruth> mctruth = mclist[0];

			const simb::MCParticle& particle = mctruth->GetNeutrino().Nu();
			const TLorentzVector& pvtx = particle.Position();
			TVector3 primary(pvtx.X(), pvtx.Y(), pvtx.Z());
	
			if (insideFidVol(pvtx) && (abs(mctruth->GetNeutrino().Lepton().PdgCode()) == 11) && tracklist.size()) 
			{	
				// mc
				TVector3 mcvtx3d(pvtx.X(), pvtx.Y(), pvtx.Z());
				TVector3 mcdir3d = findElDir(mctruth);

				// reco
				IniSegAlg recoini(tracklist, mcvtx3d); 
				recoini.FeedwithMc(mcdir3d); // mcvtx3d == primary	

				if (recoini.IsFound())
				{
					lep_dedx = 0.0; lep_dist = 0.0;
					art::Ptr<recob::Track> recotrack = recoini.GetTrk(); 
					if (!recotrack->NumberTrajectoryPoints()) return;

					art::FindManyP< recob::Hit > fb(trackListHandle, evt, fTrackModuleLabel);
					std::vector< art::Ptr<recob::Hit> > recoinihit = fb.at(recotrack.key());
					
					// use recob::Track functionality as much as possible
					const double setlength = 2.5; double length = 0.0; // cm
							
					TVector3 pos_p = recotrack->LocationAtPoint(0);
					
					float px = pos_p.X();
					if (px > 0) px -= t0Corr(evt, pvtx);
					else px += t0Corr(evt, pvtx);
					pos_p.SetX(px);
					
					vtxrecomc = std::sqrt(pma::Dist2(primary, pos_p));
					vtxrecomcx = primary.X() - pos_p.X();
					vtxrecomcy = primary.Y() - pos_p.Y();
					vtxrecomcz = primary.Z() - pos_p.Z();

					lep_dist = std::sqrt(pma::Dist2(primary, pos_p));	
					cos = recoini.GetCos();

					size_t fp = 0; bool hitcoll = false;
					for (size_t p = 0; p < recotrack->NumberTrajectoryPoints(); ++p)
						if (recotrack->DQdxAtPoint(p, geo::kZ) > 0) {pos_p = recotrack->LocationAtPoint(p); fp = p; hitcoll = true; break;}
					
					// loop over trajectory point to get dQdx.
					if (hitcoll)
						for (size_t p = (fp+1); p < recotrack->NumberTrajectoryPoints(); ++p)
						{
							TVector3 pos = recotrack->LocationAtPoint(p);
							length += std::sqrt(pma::Dist2(pos_p, pos));
							pos_p = recotrack->LocationAtPoint(p);

							if (length > setlength) break;
							dx = length;	
							double dqdx_p = recotrack->DQdxAtPoint(p, geo::kZ);
							if (dqdx_p > 0) lep_dedx += recoinihit[p]->SummedADC();
						}

					
					if (dx > 0.0) lep_dedx /= dx;
				}
				
			}
		}
	}
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

		for (size_t c = 0; c < geom->Ncryostats(); ++c) 
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

std::vector< dunefd::Hit2D > dunefd::IniSegReco::reselectCls(std::map<size_t, std::vector< dunefd::Hit2D > > const & cls, art::Ptr<simb::MCTruth> const mctruth, const size_t cryo, const size_t tpc, const size_t plane) const
{
	std::vector< dunefd::Hit2D > cluster;
	if (!cls.size()) return cluster;

	const simb::MCParticle& particle = mctruth->GetNeutrino().Nu();
	const TLorentzVector& pvtx = particle.Position();
	
	if (insideFidVol(pvtx) && (abs(mctruth->GetNeutrino().Lepton().PdgCode()) == 11))
	{
		// mc
		TVector3 mcvtx3d(pvtx.X(), pvtx.Y(), pvtx.Z());
		TVector3 mcdir3d = findElDir(mctruth);
		TVector2 mcvtx2d = getMCvtx2d(mcvtx3d, cryo, tpc, plane);
		TVector2 mcdir2d = getMCdir2d(mcvtx3d, mcdir3d, cryo, tpc, plane);
		
		// reco: find the best clusters to proceed with segment reconstruction
		IniSegAlg recoini(cls); 
		recoini.FeedwithMc(mcvtx2d, mcdir2d, mcdir3d);
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
		pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(hitsel);
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

TVector2 dunefd::IniSegReco::getMCvtx2d(TVector3 const & mcvtx3d, const size_t cryo, const size_t tpc, const size_t plane) const
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

float dunefd::IniSegReco::t0Corr(art::Event const & evt, TLorentzVector const & pvtx) 
{
	float corrt0x = 0.0F;

	art::ServiceHandle<geo::Geometry> geom;
	auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();

	// * MC truth information
	art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    	std::vector<art::Ptr<simb::MCTruth> > mclist;
    	if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      	art::fill_ptr_vector(mclist, mctruthListHandle);

	double vtx[3] = {pvtx.X(), pvtx.Y(), pvtx.Z()};

	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

	if (geom->HasTPC(idtpc))
		if (mclist.size())
		{
			art::Ptr<simb::MCTruth> mctruth = mclist[0];
			if (mctruth->NParticles())
			{
				simb::MCParticle particle = mctruth->GetParticle(0);
				t0 = particle.T(); // ns
				corrt0x = t0 * 1.e-3 * larprop->DriftVelocity();
     	 	}
		}
	

	return corrt0x;
}

/***********************************************************************/

bool dunefd::IniSegReco::insideFidVol(TLorentzVector const & pvtx) const
{
	art::ServiceHandle<geo::Geometry> geom;
	double vtx[3] = {pvtx.X(), pvtx.Y(), pvtx.Z()};
	bool inside = false;

	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

	if (geom->HasTPC(idtpc))
	{		
		const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
		double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
		double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
		double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

		for (size_t c = 0; c < geom->Ncryostats(); c++)
		{
			const geo::CryostatGeo& cryostat = geom->Cryostat(c);
			for (size_t t = 0; t < cryostat.NTPC(); t++)
			{	
				const geo::TPCGeo& tpcg = cryostat.TPC(t);
				if (tpcg.MinX() < minx) minx = tpcg.MinX();
				if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX(); 
				if (tpcg.MinY() < miny) miny = tpcg.MinY();
				if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
				if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
				if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
			}
		}	

		//x
		double dista = fabs(minx - pvtx.X());
		double distb = fabs(pvtx.X() - maxx); 
;
		if ((pvtx.X() > minx) && (pvtx.X() < maxx) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut))
		{ 
			inside = true;
		}
		else { inside = false; }
		//y
		dista = fabs(maxy - pvtx.Y());
		distb = fabs(pvtx.Y() - miny);
		if (inside && (pvtx.Y() > miny) && (pvtx.Y() < maxy) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
		else inside = false;

		//z
		dista = fabs(maxz - pvtx.Z());
		distb = fabs(pvtx.Z() - minz);
		if (inside && (pvtx.Z() > minz) && (pvtx.Z() < maxz) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
		else inside = false;
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
