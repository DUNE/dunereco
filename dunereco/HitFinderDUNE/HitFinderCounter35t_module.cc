#ifndef HITFINDERCOUNTER35T_H
#define HITFINDERCOUNTER35T_H

////////////////////////////////////////////////////////////////////////
//
// HitFinderCounter35t class
// 
// k.warburton@sheffield.ac.uk
//  
// Michelle Stancari had the idea to use the 35 ton counters to aid in
//  disambiguation. This module attempts to do that.
// First only hits within a window around the counter coincidence are
//  considered.
// Then an unambigious wire segment is searched for for the hit.
//
// Uses the normal 35 ton disambiguation module as a base.
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <utility> // std::move()
#include <memory> // std::unique_ptr<>


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Seed Service
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandomEngine.h"

// Want to include the CounterPositionMapFunction
#include "dune/daqinput35t/PennToOffline.h"

#include "HitLineFitAlg.h"

// ROOT Includes 
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVector2.h"
#include "TVector3.h"

namespace dune{
  class HitFinderCounter35t : public art::EDProducer {
    
  public:
  
    explicit HitFinderCounter35t(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt);
    void beginJob();
    void endJob();
  
  private:

    //bool DoesIntersect( double A0, double B0, double A1, double B1, double A2, double B2, double A3, double B3 );
    int  pnpoly( int nvert, float *vertx, float *verty, float testx, float testy );
    void FindXZGradient   ( std::vector < recob::Hit > HitVector, float &Gradient, float& Intercept );
    void FindPlaneGradient( std::vector < recob::Hit > HitVector, float &Gradient, float& Intercept );
    void MatrixGradient   ( TMatrixD X, TMatrixD Y, size_t MatSize, float &Gradient, float &Intercept );
    double DistToLine     ( float Gradient, float Intercept, double X, double Y);
    std::vector<recob::Hit> FindUniqueHits( std::vector < recob::Hit > const GoodHits, unsigned int Plane, unsigned int TPC );
    int  CheckWhichIndex( std::vector<int> CloseHits );
    std::vector<int>  ClosestDistance( std::vector<double> CloseHits );
    void OutAndClearVector ( std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InitialBad,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &NewBad,
			     std::vector < std::pair < recob::Hit, size_t > > & NewGood, size_t HitsSize );
    void CrossCollection   ( std::vector < recob::Hit > const GoodHits,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
			     std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector );
    void AdjacentWireWidth ( std::vector < recob::Hit > const GoodHits,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
			     std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector );
    void TwoDimXZLineFit   ( std::vector < recob::Hit > const GoodHits,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
			     std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector );
    void TwoDimPlaneLineFit( std::vector < recob::Hit > const GoodHits,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
			     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
			     std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector );
    
    art::ServiceHandle<geo::Geometry> fGeom;
  
    bool         fDebug;
    std::string  fHitsModuleLabel;
    std::string  fCounterModuleLabel;
    std::string  fCounterDir;
    std::string  fCounterFile;
    int          fCoincidenceTolerance;
    double       fConvCountTimeToTicks;
    double       fExtendCounters;
    unsigned int fAdjacentWireWidth;
    unsigned int fAdjacentTimeWidth;
    unsigned int fCollectionTimeWidth;

    dune::HitLineFitAlg fFitAlg;
    art::ServiceHandle<art::RandomNumberGenerator> fRng;
    art::ServiceHandle<rndm::NuRandomService> fSeed;
    bool fDoHitLineFitAlg;

    std::map< unsigned int, std::pair < TVector3, std::vector< TVector3 > > > CounterPositionMap; // The map of counter positions....
    std::vector< std::pair< unsigned int, unsigned int > > ExternalTrigIndexVec; // My vector of counter coincidence indexes...
    double UVStartEndPointsX[4];
    
    double fDriftVelocity;
    int    TrigTime=0;

    TH2F* TwoDLineHist;
  protected:
  }; // class HitFinderCounter35t
  
  
  //-------------------------------------------------
  //-------------------------------------------------
  HitFinderCounter35t::HitFinderCounter35t(fhicl::ParameterSet const& pset) 
    :
    fDebug                 (pset.get<bool>("Debug"))
    , fHitsModuleLabel     (pset.get< std::string >("HitsModuleLabel"   ))
    , fCounterModuleLabel  (pset.get< std::string >("CounterModuleLabel"))
    , fCounterDir          (pset.get< std::string >("CounterDir" ))
    , fCounterFile         (pset.get< std::string >("CounterFile"))
    , fCoincidenceTolerance(pset.get<int>("CoincidenceTolerance")) // num PENN board ticks
    , fConvCountTimeToTicks(pset.get<double>("fConvCountTimeToTicks"))
    , fExtendCounters      (pset.get<double>("ExtendCounters"))
    , fAdjacentWireWidth   (pset.get<unsigned int>("AdjacentWireWidth"))
    , fAdjacentTimeWidth   (pset.get<unsigned int>("AdjacentTimeWidth"))
    , fCollectionTimeWidth (pset.get<unsigned int>("CollectionTimeWidth"))
    , fFitAlg              (pset.get<fhicl::ParameterSet>("HitLineFitAlg"))
    , fDoHitLineFitAlg     (pset.get<bool>("DoHitLineFitAlg",false))
  {
    fSeed->createEngine(*this,"HepJamesRandom","Seed");
    recob::HitCollectionCreator::declare_products(*this);
  }
  //-------------------------------------------------
  void HitFinderCounter35t::beginJob() {
    
    // Want to know the X co-ords for the short and long induction wires, only want to calculate this once as all wires have the same X position.
    double UWireLongDriftStart[3], UWireLongDriftEnd[3], UWireShortDriftStart[3], UWireShortDriftEnd[3];
    double VWireLongDriftStart[3], VWireLongDriftEnd[3], VWireShortDriftStart[3], VWireShortDriftEnd[3];
    fGeom->WireEndPoints( 0, 0, 0, 0, UWireShortDriftStart, UWireShortDriftEnd );
    fGeom->WireEndPoints( 0, 1, 0, 0, UWireLongDriftStart , UWireLongDriftEnd  );
    fGeom->WireEndPoints( 0, 0, 1, 0, VWireShortDriftStart, VWireShortDriftEnd );
    fGeom->WireEndPoints( 0, 1, 1, 0, VWireLongDriftStart , VWireLongDriftEnd  );
    UVStartEndPointsX[0] = UWireShortDriftStart[0];
    UVStartEndPointsX[1] = UWireLongDriftStart[0];
    UVStartEndPointsX[2] = VWireShortDriftStart[0];
    UVStartEndPointsX[3] = VWireLongDriftStart[0];
    
    DAQToOffline::MakeCounterPositionMap(fCounterDir, fCounterFile, CounterPositionMap, fExtendCounters);
    
    art::ServiceHandle<art::TFileService> tfs;
    TwoDLineHist = tfs->make<TH2F>("TwoDLineHist","Plot of collection plane hit positions in the XZ plane; Z position (cm); X position (cm)", 320, -5, 155, 800, -50, 350 ); 
    TwoDLineHist->SetMarkerStyle(6);
    TwoDLineHist->SetMarkerSize(1);

  }
  //-------------------------------------------------
  void HitFinderCounter35t::endJob() {
  
  }
  //-------------------------------------------------
  void HitFinderCounter35t::produce(art::Event& evt) {
    
    if (fDebug) std::cout << "\n\nLooking at event " << evt.event() << std::endl;

    CLHEP::HepRandomEngine const & engine = fRng->getEngine("Seed");
    fFitAlg.SetSeed(engine.getSeed());

    // get raw::ExternalTriggers
    art::Handle< std::vector< raw::ExternalTrigger> > externalTriggerListHandle;
    std::vector< art::Ptr< raw::ExternalTrigger> > trigs;
    if (evt.getByLabel(fCounterModuleLabel, externalTriggerListHandle) )
      art::fill_ptr_vector(trigs,externalTriggerListHandle);

    // get raw::ExternalTriggers
    art::Handle< std::vector< recob::Hit> > HitListHandle;
    std::vector< art::Ptr< recob::Hit> > hits;
    if (evt.getByLabel(fHitsModuleLabel, HitListHandle) )
      art::fill_ptr_vector(hits,HitListHandle);

    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fDriftVelocity = detprop->DriftVelocity();  // cm / us
        
    // Loop through counter hits, to look for external trigger pairs.
    ExternalTrigIndexVec.clear();
    for (size_t CLoop=0; CLoop < trigs.size(); ++CLoop) {
      for (size_t CLoop2=0; CLoop2 < trigs.size(); ++CLoop2) {
	if ( fabs( trigs[CLoop]->GetTrigTime() - trigs[CLoop2]->GetTrigTime()) > fCoincidenceTolerance ) continue;
	if ( CLoop == CLoop2 ) continue;
	// for c2: GetTrigID() is an unsigned int and always >= 0
	//if ( (    trigs[CLoop]->GetTrigID() >= 6  && trigs[CLoop]->GetTrigID() <= 15 && trigs[CLoop2]->GetTrigID() >= 28 && trigs[CLoop2]->GetTrigID() <= 37 ) // East Lower, West Upper 
	 //    || ( trigs[CLoop]->GetTrigID() >= 0  && trigs[CLoop]->GetTrigID() <= 5  && trigs[CLoop2]->GetTrigID() >= 22 && trigs[CLoop2]->GetTrigID() <= 27 ) // South Lower, North Upper
	  //   || ( trigs[CLoop]->GetTrigID() >= 16 && trigs[CLoop]->GetTrigID() <= 21 && trigs[CLoop2]->GetTrigID() >= 38 && trigs[CLoop2]->GetTrigID() <= 43 ) // North Lower, South Upper
	if ( (    trigs[CLoop]->GetTrigID() >= 6  && trigs[CLoop]->GetTrigID() <= 15 && trigs[CLoop2]->GetTrigID() >= 28 && trigs[CLoop2]->GetTrigID() <= 37 ) // East Lower, West Upper 
	     || ( trigs[CLoop]->GetTrigID() <= 5  && trigs[CLoop2]->GetTrigID() >= 22 && trigs[CLoop2]->GetTrigID() <= 27 ) // South Lower, North Upper
	     || ( trigs[CLoop]->GetTrigID() >= 16 && trigs[CLoop]->GetTrigID() <= 21 && trigs[CLoop2]->GetTrigID() >= 38 && trigs[CLoop2]->GetTrigID() <= 43 ) // North Lower, South Upper
	     ) {
	  if (fDebug) std::cout << "I have a match..."
				<< "CLoop " << CLoop << ", ID " << trigs[CLoop]->GetTrigID() << ", time " << trigs[CLoop]->GetTrigTime() << "..."
				<< "CLoop2 " << CLoop2 << ", ID " << trigs[CLoop2]->GetTrigID() << ", Time " <<  trigs[CLoop2]->GetTrigTime()
				<< std::endl;
	  ExternalTrigIndexVec.push_back( std::make_pair( CLoop, CLoop2 ) );
	} // If a good coincidence
      } // CLoop2
    } // CLoop
    
    std::cout << "ExternalTrigIndexVec has size " << ExternalTrigIndexVec.size() << std::endl;
    for (size_t aa=0; aa< ExternalTrigIndexVec.size(); ++aa) {
      int IDA      = trigs.at(ExternalTrigIndexVec[aa].first)->GetTrigID();
      double TimeA = trigs.at(ExternalTrigIndexVec[aa].first)->GetTrigTime();
      double TickA = trigs.at(ExternalTrigIndexVec[aa].first)->GetTrigTime() / fConvCountTimeToTicks;
      int IDB      = trigs.at(ExternalTrigIndexVec[aa].second)->GetTrigID();
      double TimeB = trigs.at(ExternalTrigIndexVec[aa].second)->GetTrigTime();
      double TickB = trigs.at(ExternalTrigIndexVec[aa].second)->GetTrigTime() / fConvCountTimeToTicks;
      if (fDebug) std::cout << "ExternalTrigIndexVec["<<aa<<"] has indexes " << ExternalTrigIndexVec[aa].first << " and " << ExternalTrigIndexVec[aa].second
			    << ". They correspond to TrigIDs " << IDA << " and " << IDB << ", and times "  << TimeA << " (" << TickA << ") and " << TimeB << " (" << TickB << ")"
			    << std::endl;
    }

    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> ChannelHitRawDigits (HitListHandle, evt, fHitsModuleLabel);
    art::FindOneP<recob::Wire>   ChannelHitWires     (HitListHandle, evt, fHitsModuleLabel);
    
    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this, evt,
    				     /* doWireAssns */ ChannelHitWires.isValid(),
    				     /* doRawDigitAssns */ ChannelHitRawDigits.isValid()
    				     );
    
    // I want a vector of a vector of hits, to store my undisambiguated hits.
    std::vector < std::pair < std::vector < recob::Hit >, size_t > > UnDisambigHits;

    // FIXME:::::I am unsure of the best way to handle showers - lots of trigger indexes at the same time.
    //           The easiest ( and laziest ) thing to do, is to just discount events which have many triggers.
    if ( ExternalTrigIndexVec.size() > 1 ) {
      std::cout << "Voiding this event as I have a vector of size " << ExternalTrigIndexVec.size() << std::endl;
      hcol.put_into(evt);
      return;
    }

    

    // Loop through all the indexes to see whether to keep the hit.
    for (size_t TrigInd=0; TrigInd< ExternalTrigIndexVec.size(); ++TrigInd) {
      // Make my Corner arrays.
      std::map< unsigned int, std::pair < TVector3, std::vector< TVector3 > > >::iterator it1 = CounterPositionMap.find( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() );
      std::map< unsigned int, std::pair < TVector3, std::vector< TVector3 > > >::iterator it2 = CounterPositionMap.find( trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() );     
      float VertexX[4] = { (float)it1->second.second[0][0], (float)it1->second.second[1][0], (float)it2->second.second[1][0], (float)it2->second.second[0][0]};
      float VertexY[4] = { (float)it1->second.second[0][1], (float)it1->second.second[2][1], (float)it2->second.second[2][1], (float)it2->second.second[0][1]};
      float VertexZ[4] = { (float)it1->second.second[0][2], (float)it1->second.second[1][2], (float)it2->second.second[1][2], (float)it2->second.second[0][2]};
      
      float SmallestX = std::min( VertexX[0], std::min( VertexX[1], std::min( VertexX[2], VertexX[3] ) ) );
      float BiggestX  = std::max( VertexX[0], std::max( VertexX[1], std::max( VertexX[2], VertexX[3] ) ) );
      
      TrigTime = ( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigTime() + trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigTime() ) / (2*fConvCountTimeToTicks);
      ///*
      if (fDebug) std::cout << "\nCounter1, with TrigID " << trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() << " (" << it1->first << ") has corners......"
			    << "( " << it1->second.second[0][0] << ", " << it1->second.second[0][1] << ", " << it1->second.second[0][2] << ") and "
			    << "( " << it1->second.second[1][0] << ", " << it1->second.second[1][1] << ", " << it1->second.second[1][2] << ") and "
			    << "( " << it1->second.second[2][0] << ", " << it1->second.second[2][1] << ", " << it1->second.second[2][2] << ") and "
			    << "( " << it1->second.second[3][0] << ", " << it1->second.second[3][1] << ", " << it1->second.second[3][2] << ") "
			    << "\nCounter2, with TrigID " << trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() << " (" << it2->first << ") has corners......"
			    << "( " << it2->second.second[0][0] << ", " << it2->second.second[0][1] << ", " << it2->second.second[0][2] << ") and "
			    << "( " << it2->second.second[1][0] << ", " << it2->second.second[1][1] << ", " << it2->second.second[1][2] << ") and "
			    << "( " << it2->second.second[2][0] << ", " << it2->second.second[2][1] << ", " << it2->second.second[2][2] << ") and "
			    << "( " << it2->second.second[3][0] << ", " << it2->second.second[3][1] << ", " << it2->second.second[3][2] << ") "
			    <<"\nVertexX has co-ords " << VertexX[0] << " " << VertexX[1] << " " << VertexX[2] << " " << VertexX[3]
			    <<"\nVertexY has co-ords " << VertexY[0] << " " << VertexY[1] << " " << VertexY[2] << " " << VertexY[3]
			    <<"\nVertexZ has co-ords " << VertexZ[0] << " " << VertexZ[1] << " " << VertexZ[2] << " " << VertexZ[3]
			    << std::endl;
      
      bool doHitLineFitAlg = fDoHitLineFitAlg;
      int trignum = -1;
      if ( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() >= 6  && trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() <= 15 && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() >= 28 && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() <= 37 ) trignum=111; // East Lower, West Upper
      // for c2: GetTrigID() is an unsigned int and always >= 0
      //else if ( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() >= 0  && trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() <= 5  && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() >= 22 && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() <= 27 ) trignum=112; // South Lower, North Upper
      else if ( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() <= 5  && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() >= 22 && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() <= 27 ) trignum=112; // South Lower, North Upper
      else if ( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() >= 16 && trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() <= 21 && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() >= 38 && trigs.at(ExternalTrigIndexVec[TrigInd].second)->GetTrigID() <= 43 ) trignum=113; // North Lower, South Upper
      
      float c1x = it1->second.first.X();
      float c1z = it1->second.first.Z();
      float c2x = it2->second.first.X();
      float c2z = it2->second.first.Z();
      if (trignum == 111)
	{
	  float slope = ((c1x-c2x)/(c1z-c2z));
	  fFitAlg.SetParameter(0,10,-50,250);
	  fFitAlg.SetParameter(1,1,slope-0.2,slope+0.2);
	  fFitAlg.SetParameter(2,0,-0.0002,0.0002);
	  fFitAlg.SetHorizVertRanges(-10,170,-50,250);
	}
      else if (trignum == 112 || trignum == 113)
	{
	  float slope = ((c1z-c2z)/(c1x-c2x));
	  fFitAlg.SetParameter(0,10,-10,170);
	  fFitAlg.SetParameter(1,1,slope-0.2,slope+0.2);
	  fFitAlg.SetParameter(2,0,-0.0002,0.0002);
	  fFitAlg.SetHorizVertRanges(-50,250,-10,170);
	}
      else
	{
	  doHitLineFitAlg = false;
	}
      if (ExternalTrigIndexVec.size() != 1) doHitLineFitAlg = false;

      std::vector<dune::HitLineFitAlg::HitLineFitData> fitdata;
      dune::HitLineFitAlg::HitLineFitResults bestfit;
      std::map<unsigned int, size_t> index_convert;
      unsigned int fitdataindex = 0;

      //*/
      // ------------------------- Get Collection wires -------------------------
      // Loop through collection plane hits first, as can use that information to help decide which TPC an induction plane hit should be on.
      for( size_t HitInd = 0; HitInd < hits.size(); HitInd++ ) {
	if ( hits[HitInd]->PeakTime()-TrigTime < 0 ) continue; // If in PreTriggerTick window then can't deduce the X position.
	if( hits[HitInd]->View() != geo::kZ ) continue;

	// Collection plane wires mainly discriminated against using XZ plane
	bool KeepHit = false;
	//art::Ptr<recob::Wire> wire = ChannelHitWires.at(HitInd);
	//art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(HitInd);
	
	double WireStart[3], WireEnd[3];
	fGeom->WireEndPoints( hits[HitInd]->WireID(), WireStart, WireEnd );
	float DriftDist = (0.5 *(hits[HitInd]->PeakTime()-TrigTime) * fDriftVelocity );
	if ( hits[HitInd]->WireID().TPC % 2 == 0) { // If short TPC
	  DriftDist = -DriftDist;
	}
	float HitXPos = DriftDist + WireStart[0];
	
	KeepHit = pnpoly( 4, VertexX, VertexZ, HitXPos, (float)WireEnd[2] );
	if ( KeepHit ) {
	  //std::cout << "\nGood collection plane hit on Channel " << hits[HitInd]->Channel() << ", Wire " << hits[HitInd]->WireID().Wire << ", Plane " << hits[HitInd]->WireID().Plane 
	  //	      << ", TPC" << hits[HitInd]->WireID().TPC << " at time " <<  hits[HitInd]->PeakTime()-TrigTime << " has X pos " << HitXPos << " and Z pos " << (float)WireEnd[2]
	  //	      << "\nI used the following cuts... X " << VertexX[0] << " " << VertexX[1] << " " << VertexX[2] << " " << VertexX[3] 
	  //	      << ", and Z " << VertexZ[0] << " " << VertexZ[1] << " " << VertexZ[2] << " " << VertexZ[3]
	  //	      << std::endl;
	  
	  dune::HitLineFitAlg::HitLineFitData hlfd;
	  if (trignum == 111)
            {
              hlfd.hitHoriz = (float)WireEnd[2];
              hlfd.hitVert = HitXPos;
              hlfd.hitHorizErrLo = 0.25;
              hlfd.hitHorizErrHi = 0.25;
              hlfd.hitVertErrLo = 0.25;
              hlfd.hitVertErrHi = 0.25;
            }
          else if (trignum == 112 || trignum == 113)
            {
              hlfd.hitHoriz = HitXPos;
              hlfd.hitVert = (float)WireEnd[2];
              hlfd.hitHorizErrLo = 0.25;
              hlfd.hitHorizErrHi = 0.25;
              hlfd.hitVertErrLo = 0.25;
              hlfd.hitVertErrHi = 0.25;
            }
          hlfd.hitREAL = false;
          fitdata.push_back(hlfd);
	  index_convert[fitdataindex] = HitInd;
	  ++fitdataindex;

	  //hcol.emplace_back(*hits[HitInd], wire, rawdigits);
	} 
      }

      int retval = -1;
      if (doHitLineFitAlg) retval = fFitAlg.FitLine(fitdata,bestfit);

      if (!doHitLineFitAlg)
	{
	  for (unsigned int i_fd = 0; i_fd < fitdata.size(); ++i_fd)
            {
	      art::Ptr<recob::Wire> wire = ChannelHitWires.at(index_convert[i_fd]);
	      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(index_convert[i_fd]);
	      hcol.emplace_back(*hits[index_convert[i_fd]],wire,rawdigits);
	    }
	}
      else if (retval != 1 && doHitLineFitAlg) 
	{
	  std::cout << "No line fit could be found, or not enough hits in counter shadow to make a line" << std::endl;
	  hcol.put_into(evt);
	  return;
	}
      else
	{
	  for (unsigned int i_fd = 0; i_fd < fitdata.size(); ++i_fd)
	    {
	      if (fitdata[i_fd].hitREAL)
		{
		  art::Ptr<recob::Wire> wire = ChannelHitWires.at(index_convert[i_fd]);
		  art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(index_convert[i_fd]);
		  hcol.emplace_back(*hits[index_convert[i_fd]],wire,rawdigits);
		}
	    }
	}

      std::cout << "After collection wires for TrigInd " << TrigInd << ", hcol has size " << hcol.size() << std::endl;
      // Got all the collection hits, so now loop through the induction plane hits.
      // ------------------------- Get Collection wires -------------------------

      // ------------------------- Get Induction wires -------------------------
      for( size_t HitInd = 0; HitInd < hits.size(); HitInd++ ) {
	if ( hits[HitInd]->PeakTime()-TrigTime < 0 ) continue; // If in PreTriggerTick window then can't deduce the X position.
	if ( hits[HitInd]->View() == geo::kZ ) continue; // If a collection wire continue...
	//if ( hits[HitInd]->PeakAmplitude() < 8 ) continue;


	// An induction signal can be in either long or short drift volume...Work out X position for each of them...
	float DriftDist = (0.5 *(hits[HitInd]->PeakTime()-TrigTime) * fDriftVelocity );
	float HitPosShortDrift = UVStartEndPointsX[ hits[HitInd]->View()*2 ]  - DriftDist;
	float HitPosLongDrift  = UVStartEndPointsX[ 1 + hits[HitInd]->View()*2 ]  + DriftDist;
	// If outside the hit X window for both long and short, then no point considering this hit...
	bool CouldBeShortDrift = true;
	bool CouldBeLongDrift = true;	  
	if ( HitPosShortDrift < SmallestX ) CouldBeShortDrift = false;
	if ( HitPosLongDrift  > BiggestX  ) CouldBeLongDrift  = false;
	if ( !CouldBeShortDrift && !CouldBeLongDrift ) continue; 
	/*
	std::cout << "\nLooking at an induction wire hit at time " << hits[HitInd]->PeakTime()-TrigTime << ", DriftDist is " << DriftDist 
		  << " HitPosShort is " << HitPosShortDrift << " ("<<UVStartEndPointsX[ hits[HitInd]->View()*2 ]<<"), and HitPosLong is " << HitPosLongDrift << " ("<<UVStartEndPointsX[ 1 + hits[HitInd]->View()*2 ]<<")"
		  << "\nI used the following X cut " << VertexX[0] << " " << VertexX[1] << " " << VertexX[2] << " " << VertexX[3] << " min " << SmallestX << ", max " << BiggestX
		  << "\nCan hit be in the short drift? " << CouldBeShortDrift << ", and the long drift? " << CouldBeLongDrift
		  << std::endl;
	*/
	// Have deduced that my hit can be in either the short drift, the long drift or both...
	// Now want to evaluate which wire segments for this channel are in my 3D window.
	// I want to consider all wire segments because the original hit may be assigned to the wrong wire segment.
	std::vector< recob::Hit > TempHitsVec;
	int Channel = hits[HitInd]->Channel();
	std::vector< geo::WireID > WireList = fGeom->ChannelToWire( Channel );
	for (unsigned int ch=0; ch<WireList.size(); ++ch) {
	  if ( !CouldBeShortDrift && (WireList[ch].TPC%2 == 0) ) continue; // If have ruled out short drift, then no point considering any short drift TPCs
	  if ( !CouldBeLongDrift  && (WireList[ch].TPC%2 != 0) ) continue; // If have ruled out long drift, then no point considering any long drift TPCs
	  // Work out start / end positions for this wire segment.
	  double WStart[3], WEnd[3];
	  fGeom->WireEndPoints( WireList[ch], WStart, WEnd );
	  //double WireLength = pow ( (pow(WStart[1]-WEnd[1],2)) + (pow(WStart[2]-WEnd[2],2)) , 0.5 );
	  //std::cout << "Wire " << ch << " for channel " << Channel << " is in TPC " << WireList[ch].TPC << " " << " Plane " << WireList[ch].Plane << " Wire " << WireList[ch].Wire << ", length " << WireLength << std::endl;
	    
	  // I want to iterate through the wire looking if points along the wire are within the YZ projection.
	  int NPoints = 5;
	  for (int WPass=0; WPass < NPoints; ++WPass) {
	    double ThisXPos = HitPosLongDrift;
	    if (WireList[ch].TPC%2 == 0) ThisXPos = HitPosShortDrift;
	    double ThisYPos = WStart[1] + ( WPass * (WEnd[1]-WStart[1])/(NPoints-1) );
	    double ThisZPos = WStart[2] + ( WPass * (WEnd[2]-WStart[2])/(NPoints-1) );

	    bool GoodWire = false;
	    // First test if it is in the XZ window
	    GoodWire = pnpoly ( 4, VertexX, VertexZ, ThisXPos, ThisZPos );
	    if (!GoodWire) continue; // If outside XZ window
	    
	    // If have EW coincidence, want to compare in YZ.
	    // If have NS coincidence, want to compare in XY
	    if ( trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() >= 6  && trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() <= 15 ) {
	      GoodWire = pnpoly ( 4, VertexY, VertexZ, ThisYPos, ThisZPos ); // If EW do this
	    } else { // If NS do this
	      GoodWire = pnpoly ( 4, VertexX, VertexY, ThisXPos, ThisYPos );
	    }
	    if (!GoodWire) continue; // If outside window
	    	    
	    // If this hit is the original one, then add it to TempHitsVec
	    if ( WireList[ch].Wire == hits[HitInd]->WireID().Wire && WireList[ch].TPC == hits[HitInd]->WireID().TPC ) {
	      TempHitsVec.push_back( *hits[HitInd] );
	    } else { // If it isn't the origianl hit, then I need to make a new one.
	      recob::Hit hit( hits[HitInd]->Channel(), hits[HitInd]->StartTick(), hits[HitInd]->EndTick(), hits[HitInd]->PeakTime(),
			      hits[HitInd]->SigmaPeakTime(), hits[HitInd]->RMS(), hits[HitInd]->PeakAmplitude(),
			      hits[HitInd]->SigmaPeakAmplitude(), hits[HitInd]->SummedADC(), hits[HitInd]->Integral(),
			      hits[HitInd]->SigmaIntegral(), hits[HitInd]->Multiplicity(), hits[HitInd]->LocalIndex(),
			      hits[HitInd]->GoodnessOfFit(), hits[HitInd]->DegreesOfFreedom(), hits[HitInd]->View(),
			      hits[HitInd]->SignalType(), WireList[ch]
			      );
	      TempHitsVec.push_back( hit );
	    }
	    break; // Only want one hit per wire segment.
	  } // Test X positions along the wire to see if within 3D volume
	} // Loop through all possible wires for this channel
	if (TempHitsVec.size() == 0) continue; // No Wire crossings, so throw this hit out.
	if (TempHitsVec.size() == 1) { // Only one possible hit!!
	  art::Ptr<recob::Wire> wire = ChannelHitWires.at(HitInd);
	  art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(HitInd);
	  hcol.emplace_back(TempHitsVec[0], wire, rawdigits);
	} else { // If have more than 1 possible wire hit then I need to work out which one to keep...
	  UnDisambigHits.push_back( std::make_pair(TempHitsVec, HitInd) );
	} // If TempHitsVec is large
      } // Loop through hits
      std::cout << "After Induction planes for TrigInd " << TrigInd << " hcol now has size " << hcol.size() << std::endl;
    } // Loop through trigger indexes
    std::cout << "After all trigger indexes the vectors have sizes: hits -> " << hits.size() << ", hcol -> " << hcol.size() << ", UnDisambigHits -> " << UnDisambigHits.size() << std::endl;
    // ------------------------- Get Induction wires -------------------------
    
    std::vector<recob::Hit> const &Peak = hcol.peek();
    int TPC1, TPC2, TPC3, TPC4, TPC5, TPC6, TPC7, TPC0;
    TPC1 = TPC2 = TPC3 = TPC4 = TPC5 = TPC6 = TPC7 = TPC0 = 0;
    for (size_t qq=0; qq<Peak.size(); ++qq) {
      if (Peak[qq].WireID().TPC == 0) ++TPC0;
      else if (Peak[qq].WireID().TPC == 1) ++TPC1;
      else if (Peak[qq].WireID().TPC == 2) ++TPC2;
      else if (Peak[qq].WireID().TPC == 3) ++TPC3;
      else if (Peak[qq].WireID().TPC == 4) ++TPC4;
      else if (Peak[qq].WireID().TPC == 5) ++TPC5;
      else if (Peak[qq].WireID().TPC == 6) ++TPC6;
      else if (Peak[qq].WireID().TPC == 7) ++TPC7;
    }
    if (fDebug) std::cout << "\n\nI have these hits in each TPC " << TPC0 << " " << TPC1 << " " << TPC2 << " " << TPC3 << " " << TPC4 << " " << TPC5 << " " << TPC6 << " " << TPC7 << "\n\n" << std::endl;
    // Now that I have gone through the trigger indexes I want to see if I can fix some of the undisambiguated hits.
    // 1) Look for any hits where there were collection plane hits at the same time.
    // 2) Look for any hits on adjacent induction plane wires at the same time.
    
    std::vector < std::pair < std::vector <recob::Hit>, size_t > > StillUnDisambigHits; // Want to clear this each time
    std::vector < std::pair < recob::Hit, size_t > > NowGoodHits;                       // Want to clear this each time

    // ------------------------- Make a 2D line in XZ and fit --------------------------
    std::cout << "\nNow to do my next step....fit hits to a 2D line in XZ and find which ambiguous hit is closer..." << std::endl;
    std::vector < recob::Hit > const &NextGoodHits = hcol.peek();

    TwoDimXZLineFit( NextGoodHits, UnDisambigHits, StillUnDisambigHits, NowGoodHits );
    for ( size_t NowDisambig=0; NowDisambig < NowGoodHits.size(); ++NowDisambig ) {
      size_t WhichRawHit = NowGoodHits[NowDisambig].second;
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(WhichRawHit);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(WhichRawHit);
      hcol.emplace_back(NowGoodHits[NowDisambig].first, wire, rawdigits);
    }
    OutAndClearVector( UnDisambigHits, StillUnDisambigHits, NowGoodHits, hcol.size() );
    // ------------------------- Make a 2D line in XZ and fit --------------------------

    // ------------------------- Make a 2D line in XZ and fit --------------------------
    std::cout << "\nNow to do my next step....fit hits to a 2D line in each TPC / Plane combination and find which ambiguous hit is closer..." << std::endl;
    std::vector < recob::Hit > const &NewerGoodHits = hcol.peek();

    TwoDimPlaneLineFit( NewerGoodHits, UnDisambigHits, StillUnDisambigHits, NowGoodHits );
    for ( size_t NowDisambig=0; NowDisambig < NowGoodHits.size(); ++NowDisambig ) {
      size_t WhichRawHit = NowGoodHits[NowDisambig].second;
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(WhichRawHit);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(WhichRawHit);
      hcol.emplace_back(NowGoodHits[NowDisambig].first, wire, rawdigits);
    }
    OutAndClearVector( UnDisambigHits, StillUnDisambigHits, NowGoodHits, hcol.size() );
    // ------------------------- Make a 2D line in XZ and fit --------------------------

    // ------------------------- Collection Wire Crosses -------------------------
    std::cout << "\nNow to do my next step....if induction wire crosses a collection wire with a good hit..." << std::endl;
    std::vector < recob::Hit > const &NewHits = hcol.peek();
    
    CrossCollection( NewHits, UnDisambigHits, StillUnDisambigHits, NowGoodHits );
    for ( size_t NowDisambig=0; NowDisambig < NowGoodHits.size(); ++NowDisambig ) {
      size_t WhichRawHit = NowGoodHits[NowDisambig].second;
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(WhichRawHit);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(WhichRawHit);
      hcol.emplace_back(NowGoodHits[NowDisambig].first, wire, rawdigits);
    }
    OutAndClearVector( UnDisambigHits, StillUnDisambigHits, NowGoodHits, hcol.size() );
    // ------------------------- Collection Wire Crosses -------------------------

    // --------------------------- Any adjacent wires? --------------------------- Do any of the questionable hits have real hits next to them at roughly the same time?
    std::cout << "\nNow to do my next step....if any adjacent wires have hits on them..." << std::endl;
    std::vector < recob::Hit > const &NewGoodHits = hcol.peek();

    AdjacentWireWidth( NewGoodHits, UnDisambigHits, StillUnDisambigHits, NowGoodHits );
    for ( size_t NowDisambig=0; NowDisambig < NowGoodHits.size(); ++NowDisambig ) {
      size_t WhichRawHit = NowGoodHits[NowDisambig].second;
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(WhichRawHit);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(WhichRawHit);
      hcol.emplace_back(NowGoodHits[NowDisambig].first, wire, rawdigits);
    }
    OutAndClearVector( UnDisambigHits, StillUnDisambigHits, NowGoodHits, hcol.size() );
    // --------------------------- Any adjacent wires? ---------------------------

    std::cout << "\nAfter all that the vectors have sizes: hits -> " << hits.size() << ", hcol -> " << hcol.size() << std::endl;
    hcol.put_into(evt);    
  } // The produce function
  //-------------------------------------------------
  /*
  bool HitFinderCounter35t::DoesIntersect( double A0, double B0, double A1, double B1, double A2, double B2, double A3, double B3 ) {
    // Given are two lines l1=((A0, B0), (A1, B1)) and l2=((A2, B2), (A3, B3)); Ax, Bx are integers and (Ax, Bx) specify the starts and ends of the lines.
    // We want to check whether both endpoins of l1 are on different sides of l2, and both endpoints of l2 are on opposite sides of l1.
    // To check on which side of l1=((A0, B0), (A1, B1)) a point (A, B) lies, we take:
    //    an arbitrary normal vector N perpendicular to the line; one such vector is (B1-B0, A1-A0)
    //    the vector P from the start of the line to the point (A, B), which is (A-A0, B-B0)
    // We then compute the dot product:
    //     N · P = (A-A0, B-B0) · (B1-B0, A1-A0) = (A-A0) * (B1-B0) + (B-B0) * (A1-A0)
    // We're only interested in the sign: if it's positive, the point is on one side of the line; if it's negative, it's on the other.
    
    std::cout << "My counter line goes from (" << A0 << ", " << B0 << ") to (" << A1 << ", " << B1 <<"). "
	      << "My wire line goes from ("    << A2 << ", " << B2 << ") to (" << A3 << ", " << B3 <<")."
	      << std::endl;
      
    if ( ((A2-A0)*(B1-B0) - (B2-B0)*(A1-A0)) * ((A3-A0)*(B1-B0) - (B3-B0)*(A1-A0)) < 0 ) {
      if ( ((A0-A2)*(B3-B2) - (B0-B2)*(A3-A2)) * ((A1-A2)*(B3-B2) - (B1-B2)*(A3-A2)) < 0 ) {
	return true;
      } 
    }
    return false;
  } // Work out if two lines intersect each other
  */
  //-------------------------------------------------
  int HitFinderCounter35t::pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
      if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	   (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
	c = !c;
    }
    return c;
  } // Work out whether a point is within a polygon
  //-------------------------------------------------
  void HitFinderCounter35t::FindXZGradient( std::vector < recob::Hit > HitVector, float &Gradient, float& Intercept ) {
    // I want to find the line of best fit using Matrices.
    TMatrixD X(HitVector.size(), 2); // Want to reserve a size of HitVector.size()
    TMatrixD Y(HitVector.size(), 1); // Want to reserve a size of HitVector.size()
    for ( size_t HitLoop=0; HitLoop<HitVector.size(); ++HitLoop) {
      double WireEnd[3], WireStart[3];
      fGeom->WireEndPoints( HitVector[HitLoop].WireID(), WireStart, WireEnd );
      float DriftDist = (0.5 *(HitVector[HitLoop].PeakTime()-TrigTime) * fDriftVelocity );
      if ( HitVector[HitLoop].WireID().TPC % 2 == 0) { // If short TPC
	DriftDist = -DriftDist;
      }
      // Want to fit to Z on X, and X on Y
      X[HitLoop][0] = WireStart[2];
      X[HitLoop][1] = 1;
      Y[HitLoop].Assign(DriftDist + WireStart[0]);
      TwoDLineHist->Fill(WireStart[2], DriftDist + WireStart[0]);
    } // HitLoop
    MatrixGradient( X, Y, HitVector.size(), Gradient, Intercept );
  } // Find the gradient of the hits for the XZ configuration
  //-------------------------------------------------
  void HitFinderCounter35t::FindPlaneGradient( std::vector < recob::Hit > HitVector, float &Gradient, float& Intercept ) {
    TMatrixD X(HitVector.size(), 2); // Want to reserve a size of HitVector.size()
    TMatrixD Y(HitVector.size(), 1); // Want to reserve a size of HitVector.size()
    for ( size_t HitLoop=0; HitLoop<HitVector.size(); ++HitLoop) {
      X[HitLoop][0] = HitVector[HitLoop].WireID().Wire;
      X[HitLoop][1] = 1;
      Y[HitLoop][0] = HitVector[HitLoop].PeakTime();
    }
    MatrixGradient( X, Y, HitVector.size(), Gradient, Intercept );
  } // Find the gradient of hits for the TPC/Plane configuration
  //-------------------------------------------------
  void HitFinderCounter35t::MatrixGradient( TMatrixD X, TMatrixD Y, size_t MatSize, float &Gradient, float &Intercept ) {
    if (MatSize < 2) {
      std::cout << "Can't construct a line, so giving null values." << std::endl;
      Gradient = Intercept = 0;
    } else {
      //We need the transpose of X
      TMatrixD X_T(2,MatSize);
      X_T.Transpose(X);
      
      //We also need the inverse of X_T*X
      TMatrixD X_T_x_X_inv = X_T*X;
      X_T_x_X_inv.Invert();
    
      //The final equation for the params is (X^T*X)^-1 * X^T*Y
      TMatrixD params = X_T_x_X_inv * X_T * Y;
      Gradient  = params[0][0];
      Intercept = params[1][0];
      if (fDebug) std::cout << "My line is of form: y = " << Gradient << "x + " << Intercept << std::endl;
    }
  } // Work out the gradient of a line using matricies
  //-------------------------------------------------
  double HitFinderCounter35t::DistToLine ( float Gradient, float Intercept, double X, double Y) {
    // Shortest distance to a line is d = | Am + Bn + C| / sqrt( A^2 + B^2 ) where the line is Ax + By + C = 0, and the point is ( m, n )
    double Dist = fabs( (Gradient*X) - Y + Intercept ) / pow( 1 + (Gradient*Gradient), 0.5);
    return Dist;
  } // Caluclate the distance of a point to a line
  //-------------------------------------------------
  std::vector<recob::Hit> HitFinderCounter35t::FindUniqueHits( std::vector < recob::Hit > const GoodHits, unsigned int Plane, unsigned int TPC ) {
    // I want to get a vector of only Collection wire hits.
    std::vector<recob::Hit> PlaneHits;
    for ( size_t HitLoop=0; HitLoop<GoodHits.size(); ++HitLoop) {
      if ( GoodHits[HitLoop].WireID().Plane != Plane ) continue;
      if ( TPC != 10 && GoodHits[HitLoop].WireID().TPC != TPC ) continue;
      PlaneHits.emplace_back(GoodHits[HitLoop]);
    }
    // I want to make a vector of only single occupancy collection plane wires.
    std::vector<recob::Hit> UniqPlaneHits;
    for ( size_t HitLoop=0; HitLoop<PlaneHits.size(); ++HitLoop) {
      unsigned int ThisChannel = PlaneHits[HitLoop].Channel();
      int NHits = 0;
      for ( size_t HitLoop2=0; HitLoop2<PlaneHits.size(); ++HitLoop2) {
	if (PlaneHits[HitLoop2].Channel() == ThisChannel) ++NHits;
      }
      if (NHits==1) UniqPlaneHits.emplace_back(PlaneHits[HitLoop]);
    }
    
    int NLong=0;
    for (size_t Q=0; Q<UniqPlaneHits.size(); ++Q)
      if ( UniqPlaneHits[Q].WireID().TPC%2 != 0 ) ++NLong;
    double rat = (double)UniqPlaneHits.size() / (double)NLong;
    if ( rat > 0.8 ) {
      std::vector<recob::Hit> LongHit;
      for (size_t W=0; W<UniqPlaneHits.size(); ++W)
	if ( UniqPlaneHits[W].WireID().TPC%2 != 0 )
	  LongHit.emplace_back(UniqPlaneHits[W]);
      return LongHit;
    } else
      return UniqPlaneHits;
  } // Return only the hits which are unique on a given channel
  //-------------------------------------------------
  int HitFinderCounter35t::CheckWhichIndex( std::vector<int> CloseHits ) {
    int GoodHit = -1;
    // Work out if one and only one of the wire possibilities crosses collection plane hits.
    for (size_t CloseLoop=0; CloseLoop < CloseHits.size(); ++CloseLoop) {
      if (CloseHits[CloseLoop] == 0) continue;
      GoodHit = CloseLoop;
      for (size_t CloseLoop2=0; CloseLoop2 < CloseHits.size(); ++CloseLoop2) {
	if (CloseLoop == CloseLoop2 ) continue;
	if (CloseHits[CloseLoop2]!=0) GoodHit = -1;
      } //CloseLoop2
    } // CloseLoop
    return GoodHit;
  } // Work out whether one and only one index has good hits
  //-------------------------------------------------
  std::vector<int> HitFinderCounter35t::ClosestDistance( std::vector<double> CloseHits ) {
    int GoodHit = -1;
    double MinDist = DBL_MAX;
    // Work out if there is a closest wire
    for (size_t CloseLoop=0; CloseLoop < CloseHits.size(); ++CloseLoop) {
      if ( CloseHits[CloseLoop] == MinDist ) {
	GoodHit = -1;
      }
      if ( CloseHits[CloseLoop] < MinDist ) {
	MinDist = CloseHits[CloseLoop];
	GoodHit = CloseLoop;
      }
    }
    // If there is a closest then return a vector of size 1
    std::vector<int> ReturnVec;
    if (GoodHit != -1)
      ReturnVec.push_back( GoodHit );
    else { // If there are more, make a vector of the closest hits.
      for (size_t CloseLoop=0; CloseLoop < CloseHits.size(); ++CloseLoop)
	if ( CloseHits[CloseLoop] == MinDist )
	  ReturnVec.push_back( CloseLoop );
    }
    return ReturnVec;
  } // Work out which index is the closest to a line
  //-------------------------------------------------
  void HitFinderCounter35t::OutAndClearVector ( std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InitialBad,
						    std::vector < std::pair < std::vector < recob::Hit >, size_t > > &NewBad,
						    std::vector < std::pair < recob::Hit, size_t > > &NewGood, size_t HitsSize ) {
    std::cout << "The vectors have sizes: UnDisambigHits -> " << InitialBad.size() 
	      << ", StillUndisambigHits -> " << NewBad.size()
	      << ", NowGoodHits -> " << NewGood.size() << ", hcol -> " << HitsSize << std::endl;
    InitialBad = NewBad;
    NewBad.clear();
    NewGood.clear();
  } // Produce some output and clear vectors
  //-------------------------------------------------
  void HitFinderCounter35t::CrossCollection( std::vector < recob::Hit > const GoodHits,
					    std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
					    std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
					    std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector ) {
    if (!InputBadVector.size()) return; // If no bad hits are given, then I'm done already!
    for (size_t QuestHit=0; QuestHit<InputBadVector.size(); ++QuestHit) {
      std::vector<int> CloseHits;
      std::vector< recob::Hit > Hitting = InputBadVector[QuestHit].first;
      for (size_t ThisQuest=0; ThisQuest<Hitting.size(); ++ThisQuest) {
	int ThisNearHit = 0;
	for (size_t HitLoop=0; HitLoop<GoodHits.size(); ++HitLoop) {
	  if ( GoodHits[HitLoop].View() != geo::kZ ) continue;
	  if ( Hitting[ThisQuest].WireID().TPC != GoodHits[HitLoop].WireID().TPC ) continue;
	  if ( fabs( Hitting[ThisQuest].PeakTime() - GoodHits[HitLoop].PeakTime() ) >= fCollectionTimeWidth ) continue;
	  // If the collection hit is within the time window on the correct TPC etc, want to work out whether the wires in question cross...
	  geo::WireIDIntersection Intersect;
	  if (!fGeom->WireIDsIntersect( Hitting[ThisQuest].WireID(), GoodHits[HitLoop].WireID(), Intersect ) ) continue;
	  ++ThisNearHit;
	} // Loop through certain hits....HitLoop
	CloseHits.push_back(ThisNearHit);
      } // Loop through UnDisambigHits[QuestHit].....ThisQuest
      int GoodHit = CheckWhichIndex( CloseHits );
      // If only one has wire has adjacent hits then it is a good hit.
      if (GoodHit != -1) {
	OutputGoodVector.emplace_back( std::make_pair(Hitting[GoodHit], InputBadVector[QuestHit].second ) );
      } else {
	OutputBadVector.push_back(InputBadVector[QuestHit]);
      }	
    } // Loop through StillUnDisambigHits......QuestHit
  } // Checking if the wires cross collection planes with hits
  //-------------------------------------------------
  void HitFinderCounter35t::AdjacentWireWidth( std::vector < recob::Hit > const GoodHits,
					       std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
					       std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
					       std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector ) {
    if (!InputBadVector.size()) return; // If no bad hits are given, then I'm done already!
    for (size_t QuestHit=0; QuestHit<InputBadVector.size(); ++QuestHit) {
      std::vector<int> CloseHits;
      std::vector< recob::Hit > Hitting = InputBadVector[QuestHit].first;
      for (size_t ThisQuest=0; ThisQuest<Hitting.size(); ++ThisQuest) {
	int ThisNearHit = 0;
	for (size_t HitLoop=0; HitLoop<GoodHits.size(); ++HitLoop) {
	  if ( Hitting[ThisQuest].WireID().TPC != GoodHits[HitLoop].WireID().TPC ) continue;
	  if ( Hitting[ThisQuest].WireID().Plane != GoodHits[HitLoop].WireID().Plane ) continue;
	  if ( fabs( Hitting[ThisQuest].WireID().Wire - GoodHits[HitLoop].WireID().Wire ) >= fAdjacentWireWidth ) continue;
	  if ( fabs( Hitting[ThisQuest].PeakTime() - GoodHits[HitLoop].PeakTime() ) >= fAdjacentTimeWidth ) continue;
	  ++ThisNearHit;
	} // Loop through certain hits....HitLoop
	CloseHits.push_back(ThisNearHit);
      } // Loop through UnDisambigHits[QuestHit].....ThisQuest
      int GoodHit = CheckWhichIndex( CloseHits );
      // If only one has wire has adjacent hits then it is a good hit.
      if (GoodHit != -1) {
	OutputGoodVector.emplace_back( std::make_pair(Hitting[GoodHit], InputBadVector[QuestHit].second ) );
      } else {
	OutputBadVector.push_back(InputBadVector[QuestHit]);
      }
    } // Loop through UnDisambigHits......QuestHit
  } // Checking if adjacent wires have hits
  //-------------------------------------------------
  void HitFinderCounter35t::TwoDimXZLineFit( std::vector < recob::Hit > const GoodHits,
					     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
					     std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
					     std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector ) {
    if (!InputBadVector.size()) return; // If no bad hits are given, then I'm done already!

    std::vector<recob::Hit> UniqColHits = FindUniqueHits( GoodHits, 2, 10 );
    if (!UniqColHits.size()) {
      std::cout << "I have no unique collection plane wires, so can't do anything..." << std::endl;
      for (size_t q=0; q<InputBadVector.size(); ++q)
	OutputBadVector.push_back(InputBadVector[q]);
      return;
    } else {
      if (fDebug) std::cout << "I have " << UniqColHits.size() << " unique collection plane hits " << std::endl;
    }
    float Gradient, Intercept;
    FindXZGradient( UniqColHits, Gradient, Intercept );
    
    // What if I can't construct a line? Then I want to put these into my bad hits.
    if ( Gradient == 0 && Intercept == 0 ) {
      for (size_t BadLoop=0; BadLoop < InputBadVector.size(); ++BadLoop) {
	OutputBadVector.push_back( InputBadVector[BadLoop]);
      }
    } else {
      // Now use the line of best fit to determine which of the hits is the best one....
      for (size_t BadLoop=0; BadLoop < InputBadVector.size(); ++BadLoop) {
	//std::cout << "\nNow looking at BadLoop index " << BadLoop << std::endl;
	std::vector<double> CloseHits;
	std::vector< recob::Hit > Hitting = InputBadVector[BadLoop].first;
	for (size_t ThisQuest=0; ThisQuest<Hitting.size(); ++ThisQuest) {
	  //.....Determine X and Z positions.....
	  double WireEnd[3], WireStart[3];
	  fGeom->WireEndPoints( Hitting[ThisQuest].WireID(), WireStart, WireEnd );
	  float DriftDist = (0.5 *(Hitting[ThisQuest].PeakTime()-TrigTime) * fDriftVelocity );
	  if ( Hitting[ThisQuest].WireID().TPC % 2 == 0) { // If short TPC
	    DriftDist = -DriftDist;
	  }
	  double HitPosX  = DriftDist + WireStart[0];
	  double HitPosZ  = (WireStart[2] + WireEnd[2])/2; // Just take the Z of the middle of the wire...
	  double distance = DistToLine ( Gradient, Intercept, HitPosZ, HitPosX ); // I am confusingly using Z as X and X as Y....

	  CloseHits.push_back(distance);
	  //std::cout << "Looking at index " << ThisQuest << ". The hit was on channel, TPC, Plane, Wire "
	  //	  << Hitting[ThisQuest].Channel() << " " << Hitting[ThisQuest].WireID().TPC << ", "  << Hitting[ThisQuest].WireID().Plane << ", " << Hitting[ThisQuest].WireID().Wire
	  //	  << ".\nStartZ " << WireStart[2] << ", EndZ " << WireEnd[2] << ", ZPos " << HitPosZ << ", XPos = " << HitPosX << " = " << DriftDist << " + " << WireStart[0]
	  //	  << ".\nIt was " << distance << " cm from my Coll line."
	  //	  << std::endl;
	} // Each hit in BadLoop...ThisQuest
	std::vector<int> WhichIndexes = ClosestDistance(CloseHits);
	if (WhichIndexes.size() == 1 ) {
	  //std::cout << "I returned index " << WhichIndexes[0] << std::endl;
	  OutputGoodVector.emplace_back( std::make_pair(Hitting[WhichIndexes[0]], InputBadVector[BadLoop].second ) );
	} else {
	  //std::cout << "I returned indexes with size " << WhichIndexes.size() << std::endl;
	  // If haven't removed any ambiguties.
	  if ( WhichIndexes.size() == Hitting.size() ) {
	    OutputBadVector.push_back(InputBadVector[BadLoop]);
	  } else { // If I have removed at least one ambiguity.
	    std::vector<recob::Hit> BadHits;
	    for (size_t aa=0; aa<WhichIndexes.size(); ++aa) {
	      BadHits.push_back( Hitting[ WhichIndexes[aa] ] );
	    }
	    OutputBadVector.push_back( std::make_pair( BadHits, InputBadVector[BadLoop].second ) );
	  } // If removed one ambiguity
	} // If need to put in bad vector
      } // BadLoop
    } // If can make line.
  } // Two dimensional line in XZ
  //-------------------------------------------------
  void HitFinderCounter35t::TwoDimPlaneLineFit( std::vector < recob::Hit > const GoodHits,
						std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
						std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
						std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector ) {
    if (!InputBadVector.size()) return; // If no bad hits are given, then I'm done already!

    // If the ambiguties are not in the same TPC then I can't reliably tell the difference.
    // Also want to figure out which TPC, Plane combinations to do.
    std::vector< std::pair< unsigned int, unsigned int > > WhichTPCPl;
    std::vector< std::pair < std::vector < recob::Hit >, size_t > > TPCPlaneCombs;
    for (size_t HitLoop=0; HitLoop < InputBadVector.size(); ++HitLoop) {
      std::vector< recob::Hit > Hits = InputBadVector[HitLoop].first;
      unsigned int FirstTPC = Hits[0].WireID().TPC;
      for ( size_t ThisHit=0; ThisHit<Hits.size(); ++ThisHit ) {
	if ( Hits[ThisHit].WireID().TPC != FirstTPC ) {
	  if (fDebug) std::cout << "First hit is in TPC " << FirstTPC << ", whilst hit " << ThisHit << " is in TPC " << Hits[ThisHit].WireID().TPC << " giving up on this hit" << std::endl;
	  OutputBadVector.push_back( InputBadVector[HitLoop]);
	  break;
	}
      }
      // This index has all hits in same TPC, so add it to my new combination vector
      std::vector<recob::Hit> CombHits;
      for (size_t aa=0; aa<Hits.size(); ++aa) {
	CombHits.push_back( Hits[aa] );
      }
      TPCPlaneCombs.push_back( std::make_pair( CombHits, InputBadVector[HitLoop].second ) );
      // Does this represent a new TPC/Plane combination though?
      bool TPCPl = true;
      for (size_t aa=0; aa<WhichTPCPl.size(); ++aa)
	if (WhichTPCPl[aa].first == FirstTPC && WhichTPCPl[aa].second == Hits[0].WireID().Plane )
	  TPCPl = false;
      if (TPCPl)
	WhichTPCPl.push_back( std::make_pair( FirstTPC, Hits[0].WireID().Plane ) );
    } // HitLoop
    //std::cout << "After that TPCPlaneCombs has size " << TPCPlaneCombs.size() << ", OutputBadVector has size " << OutputBadVector.size() << ", and I want to look at these TPC Plane combinations." << std::endl;
    //for (size_t aa=0; aa<WhichTPCPl.size(); ++aa)
    //std::cout << "TPC " << WhichTPCPl[aa].first << ", Plane " << WhichTPCPl[aa].second << std::endl;

    for (size_t LineIndex=0; LineIndex<WhichTPCPl.size(); ++LineIndex) {
      unsigned int ThisTPC = WhichTPCPl[LineIndex].first;
      unsigned int ThisPl  = WhichTPCPl[LineIndex].second;
      std::vector<recob::Hit> UniqHits = FindUniqueHits( GoodHits, ThisPl, ThisTPC );
      if ( !UniqHits.size() ) {
	std::cout << "I have no unique hits on TPC " << ThisTPC << ", Plane " << ThisPl << std::endl;
	// I need to put these hits into my badhits vector....
	break;
      } else {
	if (fDebug) std::cout << "I have " << UniqHits.size() << " unique hits on TPC " << ThisTPC << ", Plane " << ThisPl << std::endl;
      }
      // Work out the gradient in this wire tick space.
      float Gradient, Intercept;
      FindPlaneGradient( UniqHits, Gradient, Intercept );
      
      // Got my line fit, now to work out which hit is the best.
      for (size_t HitLoop=0; HitLoop < TPCPlaneCombs.size(); ++HitLoop) {
	std::vector< recob::Hit > Hits = TPCPlaneCombs[HitLoop].first;
	// Check that I'm looking at the right TPC/Plane configuration.
	if ( Hits[0].WireID().TPC != ThisTPC || Hits[0].WireID().Plane != ThisPl ) continue;
	// What if I can't construct a line? Then I want to put these into my bad hits.
	if ( Gradient == 0 && Intercept == 0 ) {
	  OutputBadVector.push_back( TPCPlaneCombs[HitLoop]);
	  break;
	}
	int WhichIndex = -1;
	double MinDist = DBL_MAX;
	for ( size_t ThisHit=0; ThisHit<Hits.size(); ++ThisHit ) {
	  double distance = DistToLine( Gradient, Intercept, (double)Hits[ThisHit].WireID().Wire, (double)Hits[ThisHit].PeakTime() );
	  //std::cout << "Distance for " << HitLoop << " " << ThisHit << " is " << distance << std::endl;
	  if ( distance < MinDist ) { // Get the hit with the lowest distance
	    MinDist    = distance;
	    WhichIndex = ThisHit;
	  }
	} // ThisHit
	// I know which hit I want to use now.
	//std::cout << "Pushing back index " << WhichIndex << " as it was only " << MinDist << " away " << std::endl;
	OutputGoodVector.emplace_back( std::make_pair(Hits[WhichIndex], TPCPlaneCombs[HitLoop].second ) );
      } // HitLoop
    } // Loop through WhichTPCPl
  } // Two dimensional fit in wire tick space
  //-------------------------------------------------
  DEFINE_ART_MODULE(HitFinderCounter35t)
}  // end of dune namespace
#endif // COUNTERHITFINDER35T_H
