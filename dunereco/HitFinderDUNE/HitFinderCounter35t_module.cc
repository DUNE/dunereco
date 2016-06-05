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


// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RawData/raw.h"
#include "lardata/RawData/ExternalTrigger.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Wire.h"
#include "lardata/RecoBaseArt/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Want to include the CounterPositionMapFunction
#include "dune/daqinput35t/PennToOffline.h"

// ROOT Includes 
#include "TH1D.h"
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

    void MakeCounterCorners( unsigned int AuxDetindex, float Corners[4][3], double Pos[3]);
    //bool DoesIntersect( double A0, double B0, double A1, double B1, double A2, double B2, double A3, double B3 );
    int  pnpoly(int nvert, float *vertx, float *verty, float testx, float testy);
    int  CheckWhichIndex( std::vector<int> CloseHits );
    void OutputAndClearVectors ( std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InitialBad, 
				 std::vector < std::pair < std::vector < recob::Hit >, size_t > > &NewBad, 
				 std::vector < std::pair < recob::Hit, size_t > > & NewGood, size_t HitsSize );
    void CollectionView( std::vector < recob::Hit > const GoodHits,
			 std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
			 std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
			 std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector );
    void InductionView(  std::vector < recob::Hit > const GoodHits,
			 std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
			 std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
			 std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector );
  
    art::ServiceHandle<geo::Geometry> fGeom;
  
    std::string  fHitsModuleLabel;
    std::string  fCounterModuleLabel;
    std::string  fCounterDir;
    std::string  fCounterFile;
    int          fCoincidenceTolerance;
    double       fConvCountTimeToTicks;
    double       fExtendCounters;
    unsigned int fInductionWireWidth;
    unsigned int fInductionTimeWidth;
    unsigned int fCollectionTimeWidth;

    std::map< unsigned int, std::pair < TVector3, std::vector< TVector3 > > > CounterPositionMap; // The map of counter positions....
    std::vector< std::pair< unsigned int, unsigned int > > ExternalTrigIndexVec; // My vector of counter coincidence indexes...
    double UVStartEndPointsX[4];
    
    double fDriftVelocity;

  protected:
  }; // class HitFinderCounter35t
  
  
  //-------------------------------------------------
  //-------------------------------------------------
  HitFinderCounter35t::HitFinderCounter35t(fhicl::ParameterSet const& pset) 
    :
    fHitsModuleLabel       (pset.get< std::string >("HitsModuleLabel"   ))
    , fCounterModuleLabel  (pset.get< std::string >("CounterModuleLabel"))
    , fCounterDir          (pset.get< std::string >("CounterDir" ))
    , fCounterFile         (pset.get< std::string >("CounterFile"))
    , fCoincidenceTolerance(pset.get<int>("CoincidenceTolerance")) // num PENN board ticks
    , fConvCountTimeToTicks(pset.get<double>("fConvCountTimeToTicks"))
    , fExtendCounters      (pset.get<double>("ExtendCounters"))
    , fInductionWireWidth  (pset.get<unsigned int>("InductionWireWidth"))
    , fInductionTimeWidth  (pset.get<unsigned int>("InductionTimeWidth"))
    , fCollectionTimeWidth (pset.get<unsigned int>("CollectionTimeWidth"))
  {
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
  }
  //-------------------------------------------------
  void HitFinderCounter35t::endJob() {
  
  }
  //-------------------------------------------------
  void HitFinderCounter35t::produce(art::Event& evt) {
    
    std::cout << "\n\nLooking at event " << evt.event() << std::endl;

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
    fDriftVelocity = detprop->DriftVelocity(0.25, 87);  // cm / us
        
    // Loop through counter hits, to look for external trigger pairs.
    ExternalTrigIndexVec.clear();
    for (size_t CLoop=0; CLoop < trigs.size(); ++CLoop) {
      for (size_t CLoop2=0; CLoop2 < trigs.size(); ++CLoop2) {
	if ( fabs( trigs[CLoop]->GetTrigTime() - trigs[CLoop2]->GetTrigTime()) > fCoincidenceTolerance ) continue;
	if ( CLoop == CLoop2 ) continue;
	if ( (    trigs[CLoop]->GetTrigID() >= 6  && trigs[CLoop]->GetTrigID() <= 15 && trigs[CLoop2]->GetTrigID() >= 28 && trigs[CLoop2]->GetTrigID() <= 37 ) // East Lower, West Upper 
	     || ( trigs[CLoop]->GetTrigID() >= 0  && trigs[CLoop]->GetTrigID() <= 5  && trigs[CLoop2]->GetTrigID() >= 22 && trigs[CLoop2]->GetTrigID() <= 27 ) // South Lower, North Upper
	     || ( trigs[CLoop]->GetTrigID() >= 16 && trigs[CLoop]->GetTrigID() <= 21 && trigs[CLoop2]->GetTrigID() >= 38 && trigs[CLoop2]->GetTrigID() <= 43 ) // North Lower, South Upper
	     ) {
	  std::cout << "I have a match..."
		    << "CLoop " << CLoop << ", ID " << trigs[CLoop]->GetTrigID() << ", time " << trigs[CLoop]->GetTrigTime() << "..."
		    << "CLoop2 " << CLoop2 << ", ID " << trigs[CLoop2]->GetTrigID() << ", Time " <<  trigs[CLoop2]->GetTrigTime() 
		    << std::endl;
	  ExternalTrigIndexVec.push_back( std::make_pair( CLoop, CLoop2 ) );
	} // If a good coincidence
      } // CLoop2
    } // CLoop
    
    std::cout << "ExternalTrigIndexVec has size " << ExternalTrigIndexVec.size() << std::endl;
    for (size_t aa=0; aa< ExternalTrigIndexVec.size(); ++aa) {
      std::cout << "ExternalTrigIndexVec["<<aa<<"] has indexes " << ExternalTrigIndexVec[aa].first << " and " << ExternalTrigIndexVec[aa].second 
		<< ". They correspond to TrigIDs " << trigs.at(ExternalTrigIndexVec[aa].first)->GetTrigID() << " and " << trigs.at(ExternalTrigIndexVec[aa].second)->GetTrigID()
		<< ". They correspond to TrigTimes " << trigs.at(ExternalTrigIndexVec[aa].first)->GetTrigTime() << " and " << trigs.at(ExternalTrigIndexVec[aa].second)->GetTrigTime()
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
      
      int TrigTime = trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigTime() / fConvCountTimeToTicks;
      /*
      std::cout << "\nCounter1, with TrigID " << trigs.at(ExternalTrigIndexVec[TrigInd].first)->GetTrigID() << " (" << it1->first << ") has corners......"
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
      */
      // ------------------------- Get Collection wires -------------------------
      // Loop through collection plane hits first, as can use that information to help decide which TPC an induction plane hit should be on.
      for( size_t HitInd = 0; HitInd < hits.size(); HitInd++ ) {
	if ( hits[HitInd]->PeakTime()-TrigTime < 0 ) continue; // If in PreTriggerTick window then can't deduce the X position.

	// Collection plane wires mainly discriminated against using XZ plane
	if( hits[HitInd]->View() == geo::kZ ) {
	  bool KeepHit = false;
	  art::Ptr<recob::Wire> wire = ChannelHitWires.at(HitInd);
	  art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(HitInd);
	  
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
	    hcol.emplace_back(*hits[HitInd], wire, rawdigits);
	  }
	} 
      }
      std::cout << "\nAfter collection wires for TrigInd " << TrigInd << ", hcol has size " << hcol.size() << std::endl;
      // Got all the collection hits, so now loop through the induction plane hits.
      // ------------------------- Get Collection wires -------------------------

      // ------------------------- Get Induction wires -------------------------
      for( size_t HitInd = 0; HitInd < hits.size(); HitInd++ ) {
	if ( hits[HitInd]->PeakTime()-TrigTime < 0 ) continue; // If in PreTriggerTick window then can't deduce the X position.
	if ( hits[HitInd]->View() == geo::kZ ) continue; // If a collection wire continue...
	
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
	//std::cout << "The original hit is on Channel " << hits[HitInd]->Channel() << " TPC " << hits[HitInd]->WireID().TPC << " Plane " << hits[HitInd]->WireID().Plane << " Wire " << hits[HitInd]->WireID().Wire << std::endl;
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
	    double ThisYPos = WStart[1] + ( WPass * (WEnd[1]-WStart[1])/(NPoints-1) );
	    double ThisZPos = WStart[2] + ( WPass * (WEnd[2]-WStart[2])/(NPoints-1) );
	    bool GoodWire = pnpoly ( 4, VertexY, VertexZ, ThisYPos, ThisZPos );
	    if (!GoodWire) continue; // If outside YZ window
	    
	    // Now test if it is in the XZ window
	    if (WireList[ch].TPC%2 == 0) GoodWire = pnpoly ( 4, VertexX, VertexZ, HitPosShortDrift, ThisZPos );
	    else GoodWire = pnpoly ( 4, VertexX, VertexZ, HitPosLongDrift, ThisZPos );
	    if (!GoodWire) continue; // If outside XZ window
	    //std::cout << "This hit is within the 3D window! " << WireList[ch].Wire << " " << hits[HitInd]->WireID().Wire << " " << WireList[ch].TPC << " " << hits[HitInd]->WireID().TPC << std::endl;;
	    
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
      std::cout << "\nAfter Induction planes for TrigInd " << TrigInd << " hcol now has size " << hcol.size() << std::endl;
    } // Loop through trigger indexes
    std::cout << "\nAfter all trigger indexes the vectors have sizes: hits -> " << hits.size() << ", hcol -> " << hcol.size() << ", UnDisambigHits -> " << UnDisambigHits.size() << std::endl;
    // ------------------------- Get Induction wires -------------------------
    
    // Now that I have gone through the trigger indexes I want to see if I can fix some of the undisambiguated hits.
    // 1) Look for any hits where there were collection plane hits at the same time.
    // 2) Look for any hits on adjacent induction plane wires at the same time.

    std::vector < std::pair < std::vector <recob::Hit>, size_t > > StillUnDisambigHits; // Want to clear this each time
    std::vector < std::pair < recob::Hit, size_t > > NowGoodHits;                       // Want to clear this each time

    // ------------------------- Collection Wire Crosses -------------------------
    std::cout << "\nNow to do my next step....if induction wire crosses a collection wire with a good hit..." << std::endl;
    std::vector < recob::Hit > const &NewHits = hcol.peek();
    
    CollectionView( NewHits, UnDisambigHits, StillUnDisambigHits, NowGoodHits );
    for ( size_t NowDisambig=0; NowDisambig < NowGoodHits.size(); ++NowDisambig ) {
      size_t WhichRawHit = NowGoodHits[NowDisambig].second;
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(WhichRawHit);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(WhichRawHit);
      hcol.emplace_back(NowGoodHits[NowDisambig].first, wire, rawdigits);
    }
    OutputAndClearVectors( UnDisambigHits, StillUnDisambigHits, NowGoodHits, hcol.size() );
    // ------------------------- Collection Wire Crosses -------------------------

    // --------------------------- Any adjacent wires? --------------------------- Do any of the questionable hits have real hits next to them at roughly the same time?
    std::cout << "\nNow to do my next step....if an induction wire has hits on any adjacent induction wires..." << std::endl;
    std::vector < recob::Hit > const &NewGoodHits = hcol.peek();

    InductionView( NewGoodHits, UnDisambigHits, StillUnDisambigHits, NowGoodHits );
    for ( size_t NowDisambig=0; NowDisambig < NowGoodHits.size(); ++NowDisambig ) {
      size_t WhichRawHit = NowGoodHits[NowDisambig].second;
      art::Ptr<recob::Wire> wire = ChannelHitWires.at(WhichRawHit);
      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(WhichRawHit);
      hcol.emplace_back(NowGoodHits[NowDisambig].first, wire, rawdigits);
    }
    OutputAndClearVectors( UnDisambigHits, StillUnDisambigHits, NowGoodHits, hcol.size() );
    // --------------------------- Any adjacent wires? ---------------------------
    
    hcol.put_into(evt);    
  }
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
  }
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
  }
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
  }
  //-------------------------------------------------
  void HitFinderCounter35t::OutputAndClearVectors ( std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InitialBad, 
						    std::vector < std::pair < std::vector < recob::Hit >, size_t > > &NewBad, 
						    std::vector < std::pair < recob::Hit, size_t > > &NewGood, size_t HitsSize ) {
    std::cout << "The vectors have sizes: UnDisambigHits -> " << InitialBad.size() 
	      << ", StillUndisambigHits -> " << NewBad.size()
	      << ", NowGoodHits -> " << NewGood.size() << ", hcol -> " << HitsSize << std::endl;
    InitialBad = NewBad;
    NewBad.clear();
    NewGood.clear();
  }
  //-------------------------------------------------
  void HitFinderCounter35t::CollectionView( std::vector < recob::Hit > const GoodHits,
					    std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
					    std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
					    std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector ) {
    for (size_t QuestHit=0; QuestHit<InputBadVector.size(); ++QuestHit) {
      std::vector<int> CloseHits;
      std::vector< recob::Hit > Hitting = InputBadVector[QuestHit].first;
      for (size_t ThisQuest=0; ThisQuest<Hitting.size(); ++ThisQuest) {
	int ThisNearHit = 0;
	for (size_t HitLoop=0; HitLoop<GoodHits.size(); ++HitLoop) {
	  if ( GoodHits[HitLoop].View() != geo::kZ ) continue;
	  if ( Hitting[ThisQuest].WireID().TPC != GoodHits[HitLoop].WireID().TPC ) continue;
	  if ( fabs( Hitting[ThisQuest].PeakTime() - GoodHits[HitLoop].PeakTime() ) > fCollectionTimeWidth ) continue;
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
  }
  //-------------------------------------------------
  void HitFinderCounter35t::InductionView( std::vector < recob::Hit > const GoodHits,
					   std::vector < std::pair < std::vector < recob::Hit >, size_t > > &InputBadVector,
					   std::vector < std::pair < std::vector < recob::Hit >, size_t > > &OutputBadVector,
					   std::vector < std::pair < recob::Hit, size_t > > &OutputGoodVector ) {
    
    std::cout << "At the start of the induction section " << GoodHits.size() <<  " " << InputBadVector.size() <<  " " << OutputBadVector.size() << " " << OutputGoodVector.size() << std::endl;

    for (size_t QuestHit=0; QuestHit<InputBadVector.size(); ++QuestHit) {
      std::vector<int> CloseHits;
      std::vector< recob::Hit > Hitting = InputBadVector[QuestHit].first;
      for (size_t ThisQuest=0; ThisQuest<Hitting.size(); ++ThisQuest) {
	int ThisNearHit = 0;
	for (size_t HitLoop=0; HitLoop<GoodHits.size(); ++HitLoop) {
	  if ( Hitting[ThisQuest].WireID().TPC != GoodHits[HitLoop].WireID().TPC ) continue;
	  if ( Hitting[ThisQuest].WireID().Plane != GoodHits[HitLoop].WireID().Plane ) continue;
	  if ( fabs( Hitting[ThisQuest].WireID().Wire - GoodHits[HitLoop].WireID().Wire ) > fInductionWireWidth ) continue;
	  if ( fabs( Hitting[ThisQuest].PeakTime() - GoodHits[HitLoop].PeakTime() ) > fInductionTimeWidth ) continue;
	  ++ThisNearHit;
	} // Loop through certain hits....HitLoop
	CloseHits.push_back(ThisNearHit);
	std::cout << "QuestHit " << QuestHit << " was on channel " << Hitting[ThisQuest].Channel() << ", Wire segment " << ThisQuest << " has " << ThisNearHit << " near hits." << std::endl;
      } // Loop through UnDisambigHits[QuestHit].....ThisQuest
      int GoodHit = CheckWhichIndex( CloseHits );
      // If only one has wire has adjacent hits then it is a good hit.
      if (GoodHit != -1) {
	OutputGoodVector.emplace_back( std::make_pair(Hitting[GoodHit], InputBadVector[QuestHit].second ) );
      } else {
	OutputBadVector.push_back(InputBadVector[QuestHit]);
      }
    } // Loop through UnDisambigHits......QuestHit
  }
  //-------------------------------------------------
  DEFINE_ART_MODULE(HitFinderCounter35t)
}  // end of dune namespace
#endif // COUNTERHITFINDER35T_H
