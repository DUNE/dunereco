#ifndef disambigcheck_h
#define disambigcheck_h
////////////////////////////////////////////////////////////////////////
//
// Disambigcheck class
// 
// tjyang@fnal.gov
//  
// Based on Tom Junk's idea of 3-view matching
// HitFinder for 35t, input is undisambiguated hits, out is disambiguated hits
// Start from code written by talion@gmail.com
//
//
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <string>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"   
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "DisambigAlg35t.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// ROOT Includes 
#include "TH1D.h"
#include "TH2D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TCanvas.h"

namespace disambigcheck{
  class disambigcheck : public art::EDAnalyzer
  {
    
  public:
    
    explicit disambigcheck(fhicl::ParameterSet const& pset); 
    
    void analyze(const art::Event& evt); 
    void beginJob(); 
    void beginRun(const art::Run& run);
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& pset);                
    
     virtual ~disambigcheck();
    
  private:
    
    art::ServiceHandle<geo::Geometry> fGeom;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    
    std::string fChanHitDisambig;
    std::string fChanHitCheater;
    
    TH1D* fCorrect;  TH1D* fMissed;  TH1D* fIncorrect;
    TH2D* fXCorrect; TH2D* fXMissed; TH2D* fXIncorrect;
    TH2D* fYCorrect; TH2D* fYMissed; TH2D* fZIncorrect;
    TH2D* fZCorrect; TH2D* fZMissed; TH2D* fYIncorrect;
    TH2D* fThetaCorrect; TH2D* fThetaMissed; TH2D* fThetaIncorrect;
    TH2D* fPhiCorrect; TH2D* fPhiMissed; TH2D* fPhiIncorrect;

    //TCanvas* cXCorrect;
    
  protected: 
    
  }; 
  
  //-------------------------------------------------
  //-------------------------------------------------
  disambigcheck::disambigcheck(fhicl::ParameterSet const& pset) 
    : EDAnalyzer (pset)  
  {
    this->reconfigure(pset);
    //  produces< std::vector<recob::Hit> >();
  }
  
  disambigcheck::~disambigcheck() {}
  //-------------------------------------------------
  //-------------------------------------------------
  void disambigcheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    
    fChanHitDisambig =  pset.get< std::string >("ChanHitDisambig");
    fChanHitCheater  =  pset.get< std::string >("ChanHitCheater");
    
  }  
  
  //-------------------------------------------------
  //-------------------------------------------------
  void disambigcheck::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    
    fCorrect   = tfs->make<TH1D>("correct"  ,"Correct disambiguated"  , 100, 0, 1.01);
    fMissed    = tfs->make<TH1D>("missed"   ,"Missed hits"            , 100, 0, 1.01);
    fIncorrect = tfs->make<TH1D>("incorrect","Incorrect disambiguated", 100, 0, 1.01);

    fXCorrect   = tfs->make<TH2D>("Xcorrect"  ,"Correct disambiguated"  , 40, -1.01, 1.01, 40, 0, 1.01);
    fXMissed    = tfs->make<TH2D>("Xmissed"   ,"Missed hits"            , 40, -1.01, 1.01, 40, 0, 1.01);
    fXIncorrect = tfs->make<TH2D>("Xincorrect","Incorrect disambiguated", 40, -1.01, 1.01, 40, 0, 1.01);
    fYCorrect   = tfs->make<TH2D>("Ycorrect"  ,"Correct disambiguated"  , 40, -1.01, 1.01, 40, 0, 1.01);
    fYMissed    = tfs->make<TH2D>("Ymissed"   ,"Missed hits"            , 40, -1.01, 1.01, 40, 0, 1.01);
    fYIncorrect = tfs->make<TH2D>("Yincorrect","Incorrect disambiguated", 40, -1.01, 1.01, 40, 0, 1.01);
    fZCorrect   = tfs->make<TH2D>("Zcorrect"  ,"Correct disambiguated"  , 40, -1.01, 1.01, 40, 0, 1.01);
    fZMissed    = tfs->make<TH2D>("Zmissed"   ,"Missed hits"            , 40, -1.01, 1.01, 40, 0, 1.01);
    fZIncorrect = tfs->make<TH2D>("Zincorrect","Incorrect disambiguated", 40, -1.01, 1.01, 40, 0, 1.01); 

    fThetaCorrect   = tfs->make<TH2D>("Thetacorrect"  ,"Correct disambiguated"  , 20, 0    , 3.15, 40, 0, 1.01);
    fThetaMissed    = tfs->make<TH2D>("Thetamissed"   ,"Missed hits"            , 20, 0    , 3.15, 40, 0, 1.01);
    fThetaIncorrect = tfs->make<TH2D>("Thetaincorrect","Incorrect disambiguated", 20, 0    , 3.15, 40, 0, 1.01);
    fPhiCorrect     = tfs->make<TH2D>("Phicorrect"    ,"Correct disambiguated"  , 40, -3.15, 3.15, 40, 0, 1.01);
    fPhiMissed      = tfs->make<TH2D>("Phimissed"     ,"Missed hits"            , 40, -3.15, 3.15, 40, 0, 1.01);
    fPhiIncorrect   = tfs->make<TH2D>("Phiincorrect"  ,"Incorrect disambiguated", 40, -3.15, 3.15, 40, 0, 1.01);
  
    //cXCorrect = tfs->make<TCanvas>("XCorrect", "", 600, 500);
}
  
  //-------------------------------------------------
  //-------------------------------------------------
  void disambigcheck::beginRun(const art::Run& run){}
  void disambigcheck::endJob()
  {
    //cXCorrect -> cd(0);
    //fXCorrect -> Draw("colz");
    //cXCorrect -> Write();
  }
  

  //-------------------------------------------------
  void disambigcheck::analyze(const art::Event& evt)
  {
    
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
    art::Handle< std::vector<recob::Hit> > ChannelHitsDisambig;
    std::vector< art::Ptr<recob::Hit> >  ChHitsDisambig;
    if (evt.getByLabel(fChanHitDisambig, ChannelHitsDisambig))
      art::fill_ptr_vector(ChHitsDisambig, ChannelHitsDisambig);

    // Make unambiguous collection hits
    art::Handle< std::vector<recob::Hit> > ChannelHitsCheater;
    std::vector< art::Ptr<recob::Hit> >  ChHitsCheater;
    if(evt.getByLabel(fChanHitCheater, ChannelHitsCheater))
      art::fill_ptr_vector(ChHitsCheater, ChannelHitsCheater);
    
    int correcthits = 0;
    int incorrecthits = 0; 
    int missinghits = 0;
    // Map of MCParticle Track IDs, Vector 1 ( matched, incorrect, missed ),  Pair 1 ( Pair 2 (Theta, Phi) , Vector 2 direction cosines ( x, y, z ) ) 
    std::map< int, std::pair< TVector3, std::pair < std::pair< double,double >, TVector3 > > > ParticleMap;     
    if (ChannelHitsDisambig.isValid()) {
      for( size_t h = 0; h < ChHitsCheater.size(); h++ ){
	if(ChHitsCheater[h]->View() == geo::kZ ) continue;
	uint32_t cheatchannel = ChHitsCheater[h]->Channel();
	double cheatpeaktime = ChHitsCheater[h]->PeakTime();
	uint32_t cheatwire = ChHitsCheater[h]->WireID().Wire;

	std::vector<sim::IDE> ides;
	std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(ChHitsCheater[h]);
	double MaxE = 0;
	int TrackID;
	for(size_t e = 0; e < TrackIDs.size(); ++e){
	  if ( TrackIDs[e].energy > MaxE ) {
	    TrackID = TrackIDs[e].trackID;
	    MaxE    = TrackIDs[e].energy;
	  }
	}
	// Now have trackID, so get PdG code and T0 etc.
	const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
	if ( ParticleMap.count(TrackID) == 0) {
	  //std::cout << "Looking at key " << TrackID << " for the first time so setting the angles." << std::endl;
	  TVector3 MC_NormMom;
	  MC_NormMom[0] = ( particle->Px() / particle->P() );
	  MC_NormMom[1] = ( particle->Py() / particle->P() );
	  MC_NormMom[2] = ( particle->Pz() / particle->P() );
	  double MC_Theta = MC_NormMom.Theta();
	  double MC_Phi   = MC_NormMom.Phi();
	  ParticleMap[TrackID].second.second[0] = MC_NormMom[0];
	  ParticleMap[TrackID].second.second[1] = MC_NormMom[1];
	  ParticleMap[TrackID].second.second[2] = MC_NormMom[2];
	  ParticleMap[TrackID].second.first.first  = MC_Theta;
	  ParticleMap[TrackID].second.first.second = MC_Phi;
	} // If first instance of this trackID fill quantities
	bool found = false;
	for( size_t w = 0; w < ChHitsDisambig.size(); w++ ) {
	  uint32_t disambigchannel = ChHitsDisambig[w]->Channel();
	  double disambigpeaktime = ChHitsDisambig[w]->PeakTime();
	  uint32_t disambigwire = ChHitsDisambig[w]->WireID().Wire;
	
	  if((cheatchannel == disambigchannel)&&(TMath::Abs(cheatpeaktime - disambigpeaktime)< 0.1)) {
	  
	    found=true;
	    if(cheatwire==disambigwire){
	      correcthits++;   
	      ++ParticleMap[TrackID].first[0];
	    }
	    else{
	      incorrecthits++;
	      ++ParticleMap[TrackID].first[1];
	    }   
	  }	
	} // Loop through disambiguated hits
	if(!found){
	  missinghits++;
	  ++ParticleMap[TrackID].first[2];
	} 
	//std::cout << "Particle Map Angles...." << ParticleMap[TrackID].first[0] << " " << ParticleMap[TrackID].first[1] << " " << ParticleMap[TrackID].first[2] <<"\n"
	//	<< "And it now has " << ParticleMap[TrackID].second[0] << " correct, " << ParticleMap[TrackID].second[1] << " incorrect and " << ParticleMap[TrackID].second[2] << " missed hits." << std::endl;
      } // size of cheated hits

      int totalhits = (correcthits + incorrecthits + missinghits); 
           
      //double disambighitsFraction = ddisambighits / dtotalhits;
      if ( totalhits != 0 ) {
	fCorrect  ->Fill( (double)correcthits   / (double)totalhits );
	fMissed   ->Fill( (double)missinghits / (double)totalhits ); 
	fIncorrect->Fill( (double)incorrecthits   / (double)totalhits ); 
	if ( (double)correcthits / (double)totalhits == 0 ) std::cout << "WHY IS THIS EVENT NOT MATCHING ANY HITS!!!????" << std::endl;
      } // totalhits != 0
      std::cout << "\nTotal hits " << totalhits << ", correct hits " << correcthits << ", incorrect hits " << incorrecthits << ", missing hits " << missinghits << std::endl;

      for (std::map< int, std::pair< TVector3, std::pair < std::pair< double,double >, TVector3 > > >::iterator ii = ParticleMap.begin(); ii!=ParticleMap.end(); ++ii){
	int MapTotalHits = ii->second.first[0] + ii->second.first[1] + ii->second.first[2];
	std::cout << "Total hits for this TrackID " << ii->first << " is; " << MapTotalHits << ", Correct " << ii->second.first[0] << ", Incorrect " << ii->second.first[1] << ", Missed " << ii->second.first[2] << std::endl;
	if ( MapTotalHits != 0 ) {
	  fXCorrect   -> Fill ( ii->second.second.second[0], (double)ii->second.first[0] / (double)MapTotalHits );
	  fXMissed    -> Fill ( ii->second.second.second[0], (double)ii->second.first[1] / (double)MapTotalHits );
	  fXIncorrect -> Fill ( ii->second.second.second[0], (double)ii->second.first[2] / (double)MapTotalHits );
	  fYCorrect   -> Fill ( ii->second.second.second[1], (double)ii->second.first[0] / (double)MapTotalHits );
	  fYMissed    -> Fill ( ii->second.second.second[1], (double)ii->second.first[1] / (double)MapTotalHits );
	  fYIncorrect -> Fill ( ii->second.second.second[1], (double)ii->second.first[2] / (double)MapTotalHits );
	  fZCorrect   -> Fill ( ii->second.second.second[2], (double)ii->second.first[0] / (double)MapTotalHits );
	  fZMissed    -> Fill ( ii->second.second.second[2], (double)ii->second.first[1] / (double)MapTotalHits );
	  fZIncorrect -> Fill ( ii->second.second.second[2], (double)ii->second.first[2] / (double)MapTotalHits );

	  fThetaCorrect   -> Fill ( ii->second.second.first.first , (double)ii->second.first[0] / (double)MapTotalHits );
	  fThetaMissed    -> Fill ( ii->second.second.first.first , (double)ii->second.first[1] / (double)MapTotalHits );
	  fThetaIncorrect -> Fill ( ii->second.second.first.first , (double)ii->second.first[2] / (double)MapTotalHits );
	  fPhiCorrect     -> Fill ( ii->second.second.first.second, (double)ii->second.first[0] / (double)MapTotalHits );
	  fPhiMissed      -> Fill ( ii->second.second.first.second, (double)ii->second.first[1] / (double)MapTotalHits );
	  fPhiIncorrect   -> Fill ( ii->second.second.first.second, (double)ii->second.first[2] / (double)MapTotalHits );	

	  //Some cout statements to have a look at......
	  /*
	    if ( (double)ii->second.first[0] / (double)MapTotalHits < 0.8 ) {
	    std::cout << "\n\nThis track has efficiency " << (double)ii->second.first[0] / (double)MapTotalHits*100 
	    << "\nTheta is " << ii->second.second.first.first << ", and phi is " << ii->second.second.first.second 
	    << std::endl;
	    if ( ii->second.second.first.first < 1.7 && ii->second.second.first.first > 1.45 ) {	
	    std::cout << "THIS TRACK HAS LOWER EFFICIENCY AT THETA IS ROUGHLY PI/2." << std::endl;
	    }
	    if ( ii->second.second.first.second < -1.45 && ii->second.second.first.second > -1.7 ) {
	    std::cout << "THIS TRACK HAS LOWER EFFICIENCY AT PHI IS ROUGHLY -PI/2." << std::endl;
	    }
	    std::cout <<"\n\n" << std::endl;
	    }//*/
	} // MapTotalHits != 0
      } // Loop through Map
    } // If hits valid   
    ParticleMap.clear();
  }
  
  DEFINE_ART_MODULE(disambigcheck)
  
}

  
  
#endif
