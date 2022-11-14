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
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "IniSegAlg/IniSegAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

// C/C++ libraries
#include <memory>
#include <utility>
#include <fstream>
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

  void reconfigure(fhicl::ParameterSet const& p) ;

  void produce(art::Event & e) override;


private:
	void ResetVars(void);
	void chargeParticlesatVtx(art::Event const & evt);

	bool insideFidVol(TLorentzVector const & pvtx) const;
	
        float t0Corr(art::Event const & evt,
                     detinfo::DetectorPropertiesData const& detProp,
                     TLorentzVector const & pvtx);

	TVector2 getMCvtx2d(TVector3 const & mcvtx3d, 
											const size_t cryo, const size_t tpc, const size_t plane) const;
	TVector2 getMCdir2d(TVector3 const & mcvtx3d, TVector3 const & mcdir3d, 
											const size_t cryo, const size_t tpc, const size_t plane) const;	
	TVector3 findElDir(art::Ptr<simb::MCTruth> const mctruth) const;
	std::vector< TVector3 > findPhDir() const;
	std::vector< TVector3 > findDirs(art::Ptr<simb::MCTruth> const mctruth, int pdg) const;
 
        void collectCls(art::Event const & evt,
                        detinfo::DetectorPropertiesData const& detProp,
                        art::Ptr<simb::MCTruth> const mctruth);
	std::vector< dunefd::Hit2D > reselectCls(std::map<size_t, std::vector< dunefd::Hit2D > > const & cls, 
																					art::Ptr<simb::MCTruth> const mctruth,
																					const size_t cryo, const size_t tpc, const size_t plane) const;
        void make3dseg(art::Event const & evt,
                       detinfo::DetectorPropertiesData const& detProp,
                       std::vector< std::vector<Hit2D> > const & src, TVector3 const & primary);

	recob::Track convertFrom(pma::Track3D const & src);

        void UseClusters(art::Event const & evt,
                         detinfo::DetectorPropertiesData const& detProp);

        void UseTracks(art::Event const & evt,
                       detinfo::DetectorPropertiesData const& detProp);

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
	Float_t lep_distx;
	Float_t lep_disty;
	Float_t lep_distz;
	Float_t lep_distreco;
	Float_t lep_distrecox;
	Float_t lep_distrecoy;
	Float_t lep_distrecoz;
	Float_t cos;
	Float_t coslx; Float_t cosly; Float_t coslz;
	Float_t coslrx; Float_t coslry; Float_t coslrz;
	Float_t t0;
	Float_t vtxrecomc;
	Float_t vtxrecomcx;
	Float_t vtxrecomcy;
	Float_t vtxrecomcz;
	Float_t vtxupstream;
	Float_t vtxupstreamx;
	Float_t vtxupstreamy;
	Float_t vtxupstreamz;
	int ngamma;
	Float_t convdist;

	std::string fHitsModuleLabel;
	std::string fClusterModuleLabel;
	std::string fVertexModuleLabel;
	std::string fTrackModuleLabel;
	std::string fGenieGenModuleLabel;

	double fFidVolCut;
	pma::ProjectionMatchingAlg fProjectionMatchingAlg;

        std::ofstream file;
        std::ofstream file1;
};
// ------------------------------------------------------

dunefd::IniSegReco::IniSegReco(fhicl::ParameterSet const & pset)
: EDProducer{pset}, fProjectionMatchingAlg(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"))
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
	fTree->Branch("lep_distx",&lep_distx,"lep_distx/F");
	fTree->Branch("lep_disty",&lep_disty,"lep_disty/F");
	fTree->Branch("lep_distz",&lep_distz,"lep_distz/F");
	fTree->Branch("lep_distreco",&lep_distreco,"lep_distreco/F");
	fTree->Branch("lep_distrecox",&lep_distrecox,"lep_distrecox/F");
	fTree->Branch("lep_distrecoy",&lep_distrecoy,"lep_distrecoy/F");
	fTree->Branch("lep_distrecoz",&lep_distrecoz,"lep_distrecoz/F");
	fTree->Branch("cos", &cos, "cos/F");
	fTree->Branch("coslx", &coslx, "coslx/F");
	fTree->Branch("cosly", &cosly, "cosly/F");
	fTree->Branch("coslz", &coslz, "coslz/F");
	fTree->Branch("coslrx", &coslrx, "coslrx/F");
	fTree->Branch("coslry", &coslry, "coslry/F");
	fTree->Branch("coslrz", &coslrz, "coslrz/F");
	fTree->Branch("t0", &t0, "t0/F");
	fTree->Branch("vtxrecomc", &vtxrecomc, "vtxrecomc/F");
	fTree->Branch("vtxrecomcx", &vtxrecomcx, "vtxrecomcx/F");
	fTree->Branch("vtxrecomcy", &vtxrecomcy, "vtxrecomcy/F");
	fTree->Branch("vtxrecomcz", &vtxrecomcz, "vtxrecomcz/F");
	fTree->Branch("vtxupstream", &vtxupstream, "vtxupstream/F");
	fTree->Branch("vtxupstreamx", &vtxupstreamx, "vtxupstreamx/F");
	fTree->Branch("vtxupstreamy", &vtxupstreamy, "vtxupstreamy/F");
	fTree->Branch("vtxupstreamz", &vtxupstreamz, "vtxupstreamz/F");
	fTree->Branch("ngamma", &ngamma, "ngamma/I");
	fTree->Branch("convdist", &convdist, "convdist/F");


	file.open("data.dat");
	file1.open("data1.dat");
}

void dunefd::IniSegReco::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitsModuleLabel     =   pset.get< std::string >("HitsModuleLabel");
	fClusterModuleLabel  =   pset.get< std::string >("ClusterModuleLabel");
	fTrackModuleLabel	 = 	pset.get< std::string >("TrackModuleLabel");
	fVertexModuleLabel 	 =   pset.get< std::string >("VertexModuleLabel");
	fGenieGenModuleLabel =   pset.get< std::string >("GenieGenModuleLabel");

	fFidVolCut           =   pset.get< double >("FidVolCut");
	return;
}

void dunefd::IniSegReco::ResetVars()
{
	ngamma = 0;
	convdist = -10.0F;
	pmatracks.clear();
	lep_dist = -9999; //cm
	lep_distx = -9999;
	lep_disty = -9999;
	lep_distz = -9999;
	lep_distreco = -9999;
	lep_distrecox = -9999;
	lep_distrecoy = -9999;
	lep_distrecoz = -9999;
	lep_dedx = -9999;
	t0 = -9999;
	dx = -9999;
	cos = -9999;
	coslx = -9999; cosly = -9999; coslz = -9999;
	coslrx = -9999; coslry = -9999; coslrz = -9999;
	vtxrecomc = -9999;
	vtxrecomcx = -9999;
	vtxrecomcy = -9999;
	vtxrecomcz = -9999;
	vtxupstream = -9999;
	vtxupstreamx = -9999;
	vtxupstreamy = -9999;
	vtxupstreamz = -9999;
	return;
}

void dunefd::IniSegReco::produce(art::Event& evt)
{
	ResetVars();
	const art::ServiceHandle<geo::Geometry> geom;
	
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
	bool isinside = false;
	if (!isdata)
	{
		// * MC truth information
    		std::vector<art::Ptr<simb::MCTruth> > mclist;
    		auto mctruthListHandle = evt.getHandle< std::vector<simb::MCTruth> >(fGenieGenModuleLabel);
    		if (mctruthListHandle)
      		art::fill_ptr_vector(mclist, mctruthListHandle);

		if (mclist.size()) 
		{
			art::Ptr<simb::MCTruth> mctruth = mclist[0];
			const TLorentzVector& pvtx = mctruth->GetNeutrino().Nu().Position();
			primary = TVector3(pvtx.X(), pvtx.Y(), pvtx.Z());

			if (insideFidVol(pvtx)) 
			{
				isinside = true;
				chargeParticlesatVtx(evt);

				art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
				const sim::ParticleList& plist = pi_serv->ParticleList();

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
	}

	if (isinside)
	{
                auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
                UseTracks(evt, detProp);
		fTree->Fill(); 
	}

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(sp2hit));
}

/***********************************************************************/

void dunefd::IniSegReco::UseClusters(art::Event const & evt,
                                     detinfo::DetectorPropertiesData const& detProp)
{
	if (!isdata)
	{
		 // * MC truth information
    		std::vector<art::Ptr<simb::MCTruth> > mclist;
    		auto mctruthListHandle = evt.getHandle< std::vector<simb::MCTruth> >(fGenieGenModuleLabel);
    		if (mctruthListHandle)
      		art::fill_ptr_vector(mclist, mctruthListHandle);

    		if (mclist.size())
		{
      		art::Ptr<simb::MCTruth> mctruth = mclist[0];
                        collectCls(evt, detProp, mctruth);
		}
	}
}

/***********************************************************************/

void dunefd::IniSegReco::UseTracks(art::Event const & evt,
                                   detinfo::DetectorPropertiesData const& detProp)
{
	if (!isdata)
	{
		// * hits
  		std::vector<art::Ptr<recob::Hit> > hitlist;
		auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel);
 		if (hitListHandle)
    			art::fill_ptr_vector(hitlist, hitListHandle);

		// * tracks
  		std::vector<art::Ptr<recob::Track> > tracklist;
		auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  		if (trackListHandle)
    			art::fill_ptr_vector(tracklist, trackListHandle); 

		// * vertices
  		std::vector<art::Ptr<recob::Vertex> > vtxlist;
  		auto vtxListHandle = evt.getHandle< std::vector<recob::Vertex> >(fVertexModuleLabel);
  		if (vtxListHandle)
    			art::fill_ptr_vector(vtxlist, vtxListHandle);

		// * monte carlo
    		std::vector<art::Ptr<simb::MCTruth> > mclist;
    		auto mctruthListHandle = evt.getHandle< std::vector<simb::MCTruth> >(fGenieGenModuleLabel);
    		if (mctruthListHandle)
      		art::fill_ptr_vector(mclist, mctruthListHandle);

		if (mclist.size())
		{
      		art::Ptr<simb::MCTruth> mctruth = mclist[0];

			const simb::MCParticle& particle = mctruth->GetNeutrino().Nu();
			const TLorentzVector& pvtx = particle.Position();
			TVector3 primary(pvtx.X(), pvtx.Y(), pvtx.Z());

			// search for the closest reco vertex to mc primary
			TVector3 closestvtx; TVector3 minzvtx;
			if (vtxlist.size())
			{
				double xyz[3] = {0.0, 0.0, 0.0};
				vtxlist[0]->XYZ(xyz);
				float vxreco = xyz[0];
                                if (vxreco > 0) vxreco -= t0Corr(evt, detProp, pvtx);
                                else vxreco += t0Corr(evt, detProp, pvtx);
				xyz[0] = vxreco;

				TVector3 vtxreco(xyz);
				vtxreco.SetXYZ(xyz[0], xyz[1], xyz[2]);
				closestvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);	

				double mindist2 = pma::Dist2(primary, vtxreco);								
				// loop over vertices to look for the closest to primary
				for (size_t v = 1; v < vtxlist.size(); ++v)
				{
	      			vtxlist[v]->XYZ(xyz);
					float temp = xyz[0];
                                        if (temp > 0) temp -= t0Corr(evt, detProp, pvtx);
                                        else temp += t0Corr(evt, detProp, pvtx);
					xyz[0] = temp;

	      			vtxreco.SetXYZ(xyz[0], xyz[1], xyz[2]);
	      			float dist2 = pma::Dist2(primary, vtxreco);
	      			if (dist2 < mindist2)
					{
						closestvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
						mindist2 = dist2;
					}
				}

				// loop over vertices to look the most upstream 
				double minz_xyz[3] = {0.0, 0.0, 0.0};
				double minz = 9999;
				for (size_t v = 0; v < vtxlist.size(); ++v)
				{
					vtxlist[v]->XYZ(minz_xyz);
					float temp = minz_xyz[0];
                                        if (temp > 0) temp -= t0Corr(evt, detProp, pvtx);
                                        else temp += t0Corr(evt, detProp, pvtx);
					minz_xyz[0] = temp;

					if (minz_xyz[2] < minz)
					{
						minz = minz_xyz[2];
						minzvtx.SetXYZ(minz_xyz[0], minz_xyz[1], minz_xyz[2]);						
					}
				}

				// loop over tracks
				for  (size_t t = 0; t < tracklist.size(); ++t)
				{
					if (!tracklist[t]->NumberTrajectoryPoints()) continue;
					TVector3 pos_p = tracklist[t]->LocationAtPoint<TVector3>(0);

					float temp = pos_p.X();
                                        if (temp > 0) temp -= t0Corr(evt, detProp, pvtx);
                                        else temp += t0Corr(evt, detProp, pvtx);
					pos_p.SetX(temp);

					float dist2 = pma::Dist2(primary, pos_p);
					if (dist2 < mindist2)
					{
						closestvtx.SetXYZ(pos_p.X(), pos_p.Y(), pos_p.Z());
						mindist2 = dist2;
					}	

					if (pos_p.Z() < minzvtx.Z())
					{
						minzvtx.SetXYZ(pos_p.X(), pos_p.Y(), pos_p.Z());
					}			
		
					TVector3 pos_end = tracklist[t]->LocationAtPoint<TVector3>(tracklist[t]->NumberTrajectoryPoints()-1);
					temp = pos_end.X();
                                        if (temp > 0) temp -= t0Corr(evt, detProp, pvtx);
                                        else temp += t0Corr(evt, detProp, pvtx);
					pos_end.SetX(temp);

					dist2 = pma::Dist2(primary, pos_end);
					if (dist2 < mindist2)
					{
						closestvtx.SetXYZ(pos_end.X(), pos_end.Y(), pos_end.Z());
						mindist2 = dist2;
					}

					if (pos_end.Z() < minzvtx.Z())
					{
						minzvtx.SetXYZ(pos_end.X(), pos_end.Y(), pos_end.Z());
					}
				}
							
				vtxrecomc = std::sqrt(pma::Dist2(primary, closestvtx));
				vtxrecomcx = primary.X() - closestvtx.X();
				vtxrecomcy = primary.Y() - closestvtx.Y();
				vtxrecomcz = primary.Z() - closestvtx.Z();
				vtxupstream = std::sqrt(pma::Dist2(primary, minzvtx));
				vtxupstreamx = primary.X() - minzvtx.X();
				vtxupstreamy = primary.Y() - minzvtx.Y();
				vtxupstreamz = primary.Z() - minzvtx.Z();
			}
	
			// vertex region of event
			// background events
			if (insideFidVol(pvtx) && (abs(mctruth->GetNeutrino().Lepton().PdgCode()) != 11) && tracklist.size())
			{
				// mc
				TVector3 mcvtx3d(pvtx.X(), pvtx.Y(), pvtx.Z());
				std::vector< TVector3 > mcdir3d = findPhDir();

				float minlep_dist = 9999;
				for (size_t v = 0; v < mcdir3d.size(); ++v)
				{
					// reco
					IniSegAlg recoini(tracklist, mcvtx3d); 
					recoini.FeedwithMc(mcdir3d[v]);

					if (recoini.IsFound())
					{
						art::Ptr<recob::Track> recotrack = recoini.GetTrk(); 
						if (!recotrack->NumberTrajectoryPoints()) continue;

						art::FindManyP< recob::Hit > fb(trackListHandle, evt, fTrackModuleLabel);
						std::vector< art::Ptr<recob::Hit> > recoinihit = fb.at(recotrack.key());
					
						// use recob::Track functionality as much as possible
						// const double setlength = 2.5; double length = 0.0; // cm

						TVector3 pos_p = recotrack->LocationAtPoint<TVector3>(0); 
						float px = pos_p.X();
                                                if (px > 0) px -= t0Corr(evt, detProp, pvtx);
                                                else px += t0Corr(evt, detProp, pvtx);
						pos_p.SetX(px);

						float ldist = std::sqrt(pma::Dist2(primary, pos_p));
						if (ldist < minlep_dist)
						{
							minlep_dist = ldist;
							lep_dist = ldist;
							lep_distx = primary.X() - pos_p.X();
							lep_disty = primary.Y() - pos_p.Y();
							lep_distz = primary.Z() - pos_p.Z();
							if (mclist.size())
							{
								lep_distreco = std::sqrt(pma::Dist2(closestvtx, pos_p));
								lep_distrecox = closestvtx.X() - pos_p.X();
								lep_distrecoy = closestvtx.Y() - pos_p.Y();
								lep_distrecoz = closestvtx.Z() - pos_p.Z();
							}

							cos = recoini.GetCos();
							coslrx = mcdir3d[v].X(); coslry = mcdir3d[v].Y(); coslrz = mcdir3d[v].Z();

							/*************************************************************/
							/*                          WARNING                          */
							/*************************************************************/
							/* The dQdx information in recob::Track has been deprecated  */
							/* since 2016 and in 11/2018 the recob::Track interface was  */
							/* changed and DQdxAtPoint and NumberdQdx were removed.      */
							/* Therefore the code below is now commented out             */
							/* (note that it was most likely not functional anyways).    */
							/* For any issue please contact: larsoft-team@fnal.gov       */
							/*************************************************************/
							/*
							size_t fp = 0; bool hitcoll = false;
							for (size_t p = 0; p < recotrack->NumberTrajectoryPoints(); ++p)
								if (recotrack->DQdxAtPoint(p, geo::kZ) > 0) 
								{pos_p = recotrack->LocationAtPoint<TVector3>(p); fp = p; hitcoll = true; break;}
					
							// loop over trajectory point to get dQdx.
							if (hitcoll)
								for (size_t p = (fp+1); p < recotrack->NumberTrajectoryPoints(); ++p)
								{
									TVector3 pos = recotrack->LocationAtPoint<TVector3>(p);
									length += std::sqrt(pma::Dist2(pos_p, pos));
									pos_p = recotrack->LocationAtPoint<TVector3>(p);

									if (length > setlength) break;
									dx = length;	
									double dqdx_p = recotrack->DQdxAtPoint(p, geo::kZ);
									if (dqdx_p > 0) lep_dedx += recoinihit[p]->SummedADC();
								}
					
							if (dx > 0.0) lep_dedx /= dx;
							*/
							/*************************************************************/
						}
					}
				}	
			}
			else if (insideFidVol(pvtx) && (abs(mctruth->GetNeutrino().Lepton().PdgCode()) == 11) && tracklist.size()) // from nue vertex
			{	
				// mc
				TVector3 mcvtx3d(pvtx.X(), pvtx.Y(), pvtx.Z());
				TVector3 mcdir3d = findElDir(mctruth);
				coslx = mcdir3d.X(); cosly = mcdir3d.Y(); coslz = mcdir3d.Z();
				
				// reco
				IniSegAlg recoini(tracklist, mcvtx3d); 
				recoini.FeedwithMc(mcdir3d); // mcvtx3d == primary	

				if (recoini.IsFound())
				{
					art::Ptr<recob::Track> recotrack = recoini.GetTrk(); 
					if (!recotrack->NumberTrajectoryPoints()) return;

					art::FindManyP< recob::Hit > fb(trackListHandle, evt, fTrackModuleLabel);
					std::vector< art::Ptr<recob::Hit> > recoinihit = fb.at(recotrack.key());
					
					// use recob::Track functionality as much as possible
					// const double setlength = 2.5; double length = 0.0; // cm

					TVector3 pos_p = recotrack->LocationAtPoint<TVector3>(0);
					float px = pos_p.X();
                                        if (px > 0) px -= t0Corr(evt,detProp,  pvtx);
                                        else px += t0Corr(evt, detProp, pvtx);
					pos_p.SetX(px);

					lep_dist = std::sqrt(pma::Dist2(primary, pos_p));
					lep_distx = primary.X() - pos_p.X();
					lep_disty = primary.Y() - pos_p.Y();
					lep_distz = primary.Z() - pos_p.Z();
					if (mclist.size())
					{
						lep_distreco = std::sqrt(pma::Dist2(closestvtx, pos_p));
						lep_distrecox = closestvtx.X() - pos_p.X();
						lep_distrecoy = closestvtx.Y() - pos_p.Y();
						lep_distrecoz = closestvtx.Z() - pos_p.Z();
					}				
	
					cos = recoini.GetCos();
					coslrx = mcdir3d.X(); coslry = mcdir3d.Y(); coslrz = mcdir3d.Z();

					/*************************************************************/
					/*                          WARNING                          */
					/*************************************************************/
					/* The dQdx information in recob::Track has been deprecated  */
					/* since 2016 and in 11/2018 the recob::Track interface was  */
					/* changed and DQdxAtPoint and NumberdQdx were removed.      */
					/* Therefore the code below is now commented out             */
					/* (note that it was most likely not functional anyways).    */
					/* For any issue please contact: larsoft-team@fnal.gov       */
					/*************************************************************/
					/*
					size_t fp = 0; bool hitcoll = false;
					for (size_t p = 0; p < recotrack->NumberTrajectoryPoints(); ++p)
						if (recotrack->DQdxAtPoint(p, geo::kZ) > 0) {pos_p = recotrack->LocationAtPoint<TVector3>(p); fp = p; hitcoll = true; break;}
					
					// loop over trajectory point to get dQdx.
					if (hitcoll)
						for (size_t p = (fp+1); p < recotrack->NumberTrajectoryPoints(); ++p)
						{
							TVector3 pos = recotrack->LocationAtPoint<TVector3>(p);
							length += std::sqrt(pma::Dist2(pos_p, pos));
							pos_p = recotrack->LocationAtPoint<TVector3>(p);

							if (length > setlength) break;
							dx = length;	
							double dqdx_p = recotrack->DQdxAtPoint(p, geo::kZ);
							if (dqdx_p > 0) lep_dedx += recoinihit[p]->SummedADC();
						}
					if (dx > 0.0) lep_dedx /= dx;
					*/
					/*************************************************************/
				}
				
			}
		}
	}
}

/***********************************************************************/

void dunefd::IniSegReco::collectCls(art::Event const & evt,
                                    detinfo::DetectorPropertiesData const& detProp,
                                    art::Ptr<simb::MCTruth> const mctruth)
{
	art::ServiceHandle<geo::Geometry> geom;
	
	// * clusters
 	std::vector<art::Ptr<recob::Cluster> > clusterlist;
  	auto clusterListHandle = evt.getHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);
	if (clusterListHandle)
	{
   		art::fill_ptr_vector(clusterlist, clusterListHandle);
		art::FindManyP< recob::Hit > hc(clusterListHandle, evt, fClusterModuleLabel);		

                for (auto const& tpcg : geom->Iterate<geo::TPCGeo>())
		{
                        auto const [c, t] = std::make_pair(tpcg.ID().Cryostat, tpcg.ID().TPC);
				std::vector< std::vector< dunefd::Hit2D > > clinput;
                        for (size_t p = 0; p < tpcg.Nplanes(); ++p)
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


                                                if (cryo == c && tpc == t && plane == p)
							{
								found = true;
                                                                TVector2 point = pma::WireDriftToCm(detProp, wire, hitscl[h]->PeakTime(), plane, tpc, cryo);
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
                                        make3dseg(evt, detProp, clinput, primary);
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

void dunefd::IniSegReco::make3dseg(
                art::Event const & evt,
                detinfo::DetectorPropertiesData const& detProp,
		std::vector< std::vector< dunefd::Hit2D > > const & src, 
		TVector3 const & primary)
{	
	if (src.size() < 2) return;
	if ((src[0].size() < 2) || (src[1].size() < 2)) return;

  	std::vector<art::Ptr<recob::Hit> > hitlist;
	auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel);
 	if (hitListHandle)
    		art::fill_ptr_vector(hitlist, hitListHandle);

	std::vector< art::Ptr<recob::Hit> > hitsel;

	for (size_t i = 0; i < src.size(); ++i)
		for (size_t h = 0; h < src[i].size(); ++h)
			hitsel.push_back(hitlist[src[i][h].GetKey()]);

	if (hitsel.size() > 5)
	{
                pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(detProp, hitsel);
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

void dunefd::IniSegReco::chargeParticlesatVtx(art::Event const & evt)
{
	std::cout << " chargeParticleatVtx " << std::endl;
	art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	const sim::ParticleList& plist = pi_serv->ParticleList();

	TVector3 pospri; bool primary = false;
	for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
	{
		const simb::MCParticle* particle = ipar->second;

		if (particle->Process() != "primary") continue;
		
		pospri.SetXYZ(particle->Position(0).X(), particle->Position(0).Y(), particle->Position(0).Z());
		primary = true; break;
	}

	if (!primary) return;

	size_t ninteractions = 0;
	for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
	{
		std::vector< double > temp;
		const simb::MCParticle* particle = ipar->second;

		std::cout << " pdg " << particle->PdgCode() << " " << particle->P(0) << std::endl;
		std::cout << " position " << particle->Position(0).X() << ", " << particle->Position(0).Y() << ", " << particle->Position(0).Z() << std::endl;
	
		const TLorentzVector& startpos = particle->Position(0);
		TVector3 start(startpos.X(), startpos.Y(), startpos.Z());
		const TLorentzVector& endpos = particle->EndPosition();
		TVector3 stop(endpos.X(), endpos.Y(), endpos.Z());

		double diststartpri = std::sqrt(pma::Dist2(start, pospri));
		double diststoppri = std::sqrt(pma::Dist2(stop, pospri));
		if ((diststartpri > 0.1) && (diststartpri < 1.0))
		{
			//if (diststoppri < 3.0)
			//{
				double ek = particle->E(0) - particle->Mass();
				file << evt.run() << " " << evt.id().event() << " " << particle->PdgCode() << " " << diststartpri << " " << diststoppri << " " << particle->P(0) << " " << particle->Mass() << " " << ek << std::endl;
			
				if (ek > 0.05) ninteractions++;
			//}
		}
	}

	file1 << evt.run() << " " << evt.id().event() << " " << ninteractions << std::endl;
}

/***********************************************************************/

std::vector< TVector3 > dunefd::IniSegReco::findPhDir() const
{
	std::vector< TVector3 > phdirs;	
	art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	const sim::ParticleList& plist = pi_serv->ParticleList();

	for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar)
	{
		const simb::MCParticle* particle = ipar->second;
		if (particle->Process() != "primary") continue;

		TVector3 dir(0, 0, 0);
		if (particle->PdgCode() == 111)
		{
			if (particle->NumberDaughters() != 2) continue;

			const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P(particle->Daughter(0));
			if (daughter1->PdgCode() != 22) continue;

			const simb::MCParticle* daughter2 = pi_serv->TrackIdToParticle_P(particle->Daughter(1));
			if (daughter2->PdgCode() != 22) continue; 
			
			if (daughter1->EndProcess() == "conv")
			{
				TLorentzVector mom = pi_serv->TrackIdToParticle_P(particle->Daughter(0))->Momentum();
				TVector3 momvec3(mom.Px(), mom.Py(), mom.Pz());
				dir = momvec3 * (1 / momvec3.Mag());
				phdirs.push_back(dir);
			}
 
			if (daughter2->EndProcess() == "conv")
			{
				TLorentzVector mom = pi_serv->TrackIdToParticle_P(particle->Daughter(1))->Momentum();
				TVector3 momvec3(mom.Px(), mom.Py(), mom.Pz());
				dir = momvec3 * (1 / momvec3.Mag());
				phdirs.push_back(dir);
			} 
		}
		else if (particle->PdgCode() == 22)
		{
			TLorentzVector mom = particle->Momentum();
			TVector3 momvec3(mom.Px(), mom.Py(), mom.Pz());
			dir = momvec3 * (1 / momvec3.Mag());
			phdirs.push_back(dir);
		}
	}	

	return phdirs;
}

/***********************************************************************/

std::vector< TVector3 > dunefd::IniSegReco::findDirs(art::Ptr<simb::MCTruth> const mctruth, int pdg) const
{
	std::vector< TVector3 > dirs;

	art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	const sim::ParticleList& plist = pi_serv->ParticleList();

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

float dunefd::IniSegReco::t0Corr(art::Event const & evt,
                                 detinfo::DetectorPropertiesData const& detProp,
                                 TLorentzVector const & pvtx)
{
	float corrt0x = 0.0F;

	art::ServiceHandle<geo::Geometry> geom;

	// * MC truth information
    	std::vector<art::Ptr<simb::MCTruth> > mclist;
	auto mctruthListHandle = evt.getHandle< std::vector<simb::MCTruth> >(fGenieGenModuleLabel);
    	if (mctruthListHandle)
      	  art::fill_ptr_vector(mclist, mctruthListHandle);

        geo::Point_t const vtx{pvtx.X(), pvtx.Y(), pvtx.Z()};

	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

	if (geom->HasTPC(idtpc))
		if (mclist.size())
		{
			art::Ptr<simb::MCTruth> mctruth = mclist[0];
			if (mctruth->NParticles())
			{
				simb::MCParticle particle = mctruth->GetParticle(0);
				t0 = particle.T(); // ns
				corrt0x = t0 * 1.e-3 * detProp.DriftVelocity();
     	 	}
		}
	

	return corrt0x;
}

/***********************************************************************/

bool dunefd::IniSegReco::insideFidVol(TLorentzVector const & pvtx) const
{
	art::ServiceHandle<geo::Geometry> geom;
        geo::Point_t const vtx{pvtx.X(), pvtx.Y(), pvtx.Z()};
	bool inside = false;

	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

	if (geom->HasTPC(idtpc))
	{		
		const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
		double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
		double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
		double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

		/*for (size_t c = 0; c < geom->Ncryostats(); c++)
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
		}*/	

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
	return recob::Track(recob::TrackTrajectory(recob::tracking::convertCollToPoint(xyz),
						   recob::tracking::convertCollToVector(dircos),
						   recob::Track::Flags_t(xyz.size()), false),
			    0, -1., 0, recob::tracking::SMatrixSym55(), recob::tracking::SMatrixSym55(), trkindex);
}

/***********************************************************************/

DEFINE_ART_MODULE(dunefd::IniSegReco)
