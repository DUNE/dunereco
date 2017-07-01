#ifndef HITREPEATER__H
#define HITREPEATER__H

////////////////////////////////////////////////////////////////////////
//
// Shower Hit Separator
// 
// leigh.howard.whitehead@cern.ch June 2017
//  
// This is a very simple module that takes the output from the EM
// separating CNN and applies a cut on the track-like nature
// of each hit. Two hit selections are produced, one containing
// track-like hits and one containing shower-like hits. This should
// work for any MVA method that associates some sort of weight to
// hit objects.
//
// Can optionally produce an output tree containing the information
// to visualise the separation of hits. 
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <utility> 
#include <memory>  
#include <iostream>
#include <map>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft Includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "TTree.h"

namespace shs {
  class ShowerHitSeparator : public art::EDProducer {
  public:
    explicit ShowerHitSeparator(fhicl::ParameterSet const& pset);
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt);
    void beginJob(){};
    void endJob(){};

  private:
		template<size_t N> void ReadMVA(std::vector<float> &weights, art::Event& evt);
		void SaveTree(std::vector<recob::Hit> const& shwHits, std::vector<recob::Hit> const& trkHits);

    std::string fMVALabel;
    std::string fHitLabel;
		double fMVAOutputCut;
		bool fSaveTree;	

		std::string fMVAClusterLabel; // Get this from the MVA itself.

		// Output tree and its variables
		TTree *fOutTree;
		Int_t fOutEvent;
		Int_t fOutRun;
		Int_t fOutSubrun;
		Int_t fOutTPC;
		Int_t fOutCryo;
		Int_t fOutPlane;
		Int_t fOutTime;
		Int_t fOutWire;
		bool fOutIsShw;	
		bool fOutIsTrk;	

		art::Handle<std::vector<recob::Cluster> > fMVAClusters;
  };

  // implementation

  ShowerHitSeparator::ShowerHitSeparator(fhicl::ParameterSet const& pset) {

    this->reconfigure(pset);

    // Let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
		recob::HitCollectionCreator::declare_products(*this,"showerhits",true,true);
    recob::HitCollectionCreator::declare_products(*this,"trackhits",true,true);

		// Define some histograms if we want the plots
		if(fSaveTree){

			art::ServiceHandle<art::TFileService> tfs;

			fOutTree = tfs->make<TTree>("cnnHitTree","cnnHitTree");

			fOutTree->Branch("Event",&fOutEvent,"Event/I");			
			fOutTree->Branch("Run",&fOutRun,"Run/I");			
			fOutTree->Branch("Subrun",&fOutSubrun,"Subrun/I");			
			fOutTree->Branch("TPC",&fOutTPC,"TPC/I");			
			fOutTree->Branch("Cryostat",&fOutCryo,"Cryostat/I");			
			fOutTree->Branch("Plane",&fOutPlane,"Plane/I");			
			fOutTree->Branch("Time",&fOutTime,"Time/I");			
			fOutTree->Branch("Wire",&fOutWire,"Wire/I");			
			fOutTree->Branch("IsShower",&fOutIsShw,"IsShower/O");			
			fOutTree->Branch("IsTrack",&fOutIsTrk,"IsTrack/O");			
		}
  }

  void ShowerHitSeparator::reconfigure(fhicl::ParameterSet const& p) {

		// Get the parameters from the fcl file
    fMVALabel = p.get<std::string>("MVALabel");
		fHitLabel = p.get<std::string>("HitLabel");
		fMVAOutputCut = p.get<double>("MVAOutputCut");
		fSaveTree = p.get<bool>("SaveTree");

  }

	// Access the information from the MVA method. Very similar to the code in PMAlgTrackMaker_module
	template<size_t N> void ShowerHitSeparator::ReadMVA(std::vector<float> &weights, art::Event& evt){
		// Access the MVA. It could have 4 or 3 outputs so try to see which.
    auto cluResults = anab::MVAReader<recob::Cluster,N>::create(evt,fMVALabel);
    if (!cluResults){
			// We don't seem to have an MVA.
			return;
		}

		// This is the best way to get hold of the MVA clusters and tag used to associate the clusters
		// back to the underlying hit objects.
		fMVAClusters = cluResults->dataHandle();
		fMVAClusterLabel = cluResults->dataTag();

		// Now we have the MVA, we need to access the outputs
    int trkLikeIdx = cluResults->getIndex("track");
    int emLikeIdx = cluResults->getIndex("em");
    if ((trkLikeIdx < 0) || (emLikeIdx < 0)) {
			// If the variables we want don't exist then we have a problem.
			return; 
		}

		// Actually extract the weights
    const auto & cnnOuts = cluResults->outputs();
    for (size_t i = 0; i < cnnOuts.size(); ++i){

        double trkOrEm = cnnOuts[i][trkLikeIdx] + cnnOuts[i][emLikeIdx];
				double val = 0;
        if (trkOrEm > 0){
					// Make sure output is normalised to fall between 0 and 1. 
					val = cnnOuts[i][trkLikeIdx] / trkOrEm; 
				}
				weights.push_back(val);
    }

	}

  void ShowerHitSeparator::produce(art::Event& evt) {
		
		// We want to produce two hit collections
		recob::HitCollectionCreator showerHits(*this,evt,"showerhits",true,true);
		recob::HitCollectionCreator trackHits(*this,evt,"trackhits",true,true);

		// These are the MVA weights
    std::vector<float> mvaWeights;

		// Access the MVA. It could have 4 or 3 outputs so try to see which.
		ReadMVA<4>(mvaWeights,evt);
		if(mvaWeights.size() == 0){	
			// If we didn't get any weights then try the version with three outputs instead
			ReadMVA<3>(mvaWeights,evt);
		}

		// If we get some weights then separate the hits into two collections
		if(mvaWeights.size() != 0){

    	const art::FindManyP<recob::Hit> hitsFromClusters(fMVAClusters, evt,fMVAClusterLabel);

			// We also need the hit selection that was used by the CNN
			art::Handle<std::vector<recob::Hit> > inputHits;
			evt.getByLabel(fHitLabel,inputHits);
			art::FindOneP<raw::RawDigit> rawDigits(inputHits,evt,fHitLabel);
			art::FindOneP<recob::Wire> recoWires(inputHits,evt,fHitLabel);

			// At this stage we have the hits from the clusters and the tag for each cluster.
			// Loop over the clusters and add the hits to the shower-like hit collection.
    	for (size_t c = 0; c != fMVAClusters->size(); ++c){

				// Get the hits from this cluster
    	  auto const& hits = hitsFromClusters.at(c);

				// Now we must copy each hit and reassociate it to the raw digits and wires
				for (auto const & h : hits){
					// The art pointer key gives us the position of "h" in the hit collection
					recob::Hit thisHit = (*h);
					auto digit = rawDigits.at(h.key());
					auto wire = recoWires.at(h.key());
					// Add the hit and its associated raw digit and wire to the collection.
					if(mvaWeights[c] < fMVAOutputCut){
						showerHits.emplace_back(thisHit,wire,digit);
					}
					else{
						trackHits.emplace_back(thisHit,wire,digit);
					}
				} // End loop over hits
			} // End for loop over clusters
		}

		std::cout << "CNN output splitter: " << std::endl;
		std::cout << " - Found " << showerHits.size() << " shower-like hits" << std::endl;
		std::cout << " - Found " << trackHits.size() << " track-like hits" << std::endl;

		// Fill the output tree if requested to do so.	
		if(fSaveTree){
			fOutEvent = evt.event();
			fOutRun = evt.run();
			fOutSubrun = evt.subRun();
			SaveTree(showerHits.peek(),trackHits.peek());
		}

    // Put the hit collections and associations into the event
		showerHits.put_into(evt);
		trackHits.put_into(evt);

  }

	void ShowerHitSeparator::SaveTree(std::vector<recob::Hit> const& shwHits, std::vector<recob::Hit> const& trkHits){

		// Event, Run and Subrun set outside of this function. Take care of everything else here.

		for(auto const shwHit : shwHits){
    	fOutTPC = shwHit.WireID().TPC;
    	fOutCryo = shwHit.WireID().Cryostat;
    	fOutPlane = shwHit.WireID().planeID().Plane;
    	fOutTime = shwHit.StartTick();
    	fOutWire = shwHit.WireID().Wire;
    	fOutIsShw = true;
    	fOutIsTrk = false;

			fOutTree->Fill();
		}

		for(auto const trkHit : trkHits){
    	fOutTPC = trkHit.WireID().TPC;
    	fOutCryo = trkHit.WireID().Cryostat;
    	fOutPlane = trkHit.WireID().planeID().Plane;
    	fOutTime = trkHit.StartTick();
    	fOutWire = trkHit.WireID().Wire;
    	fOutIsShw = false;
    	fOutIsTrk = true;
			
			fOutTree->Fill();
		}

	}

  DEFINE_ART_MODULE(ShowerHitSeparator)

} // end of shs namespace
#endif 

