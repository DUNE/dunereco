////////////////////////////////////////////////////////////////////////
//
// \fileAPAGeometryAlg.cxx
//
// tylerdalion@gmail.com
//
//  Geometry interface to smooth over the awkwardness of fitting an
//  APA into LArSoft. It is geared towards making reconstruction simpler
//
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "APAGeometryAlg.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>


namespace apa{

  APAGeometryAlg::APAGeometryAlg(fhicl::ParameterSet const& pset) 
  {
    this->reconfigure(pset); 
    this->Init();
  }


  //----------------------------------------------------------
  APAGeometryAlg::APAGeometryAlg()
  {
    this->Init();
  }

  //----------------------------------------------------------
  APAGeometryAlg::~APAGeometryAlg()
  {
  }

  //----------------------------------------------------------
  void APAGeometryAlg::reconfigure(fhicl::ParameterSet const& /*p*/)
  {
  }

  //----------------------------------------------------------
  void APAGeometryAlg::Init()
  {

    // find the number of channels per APA
    uint32_t channel = 0;
    while( fGeom->ChannelToWire(channel+1)[0].TPC < 2 ) channel++;
    fChannelsPerAPA = channel + 1;

    // Step through channel c and find the view boundaries, until 
    // outside of first APA - these help optimize ChannelToAPAView
    // (very dependent on the conventions implimented in the channel map)
    fFirstU = 0;
    uint32_t c = 1;
    geo::WireID wid = (fGeom->ChannelToWire(c))[0];
    geo::WireID lastwid;
    while( wid.TPC < 2 ){

      if( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
	fLastU = c-1;
	fFirstV = c;   }
      
      if( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
	fLastV = c-1;
	fFirstZ0 = c;  }
      
      if( wid.TPC == lastwid.TPC + 1  ){
	fLastZ0 = c-1;
	fFirstZ1 = c;  }

      lastwid = wid;
      c++;
      wid = (fGeom->ChannelToWire(c))[0]; // for the while condition

    }

    fLastZ1 = c - 1;


    if( fLastZ1 + 1 != fChannelsPerAPA ) throw cet::exception("APAGeometryAlg")
					   << "Channel boundaries are inconsistent.\n";

    // some other things that will be needed
    fAPAsPerCryo = fGeom->NTPC(0)/2;
    fChannelRange[0] = (fLastU-fFirstU + 1)*fGeom->WirePitch(geo::kU);
    fChannelRange[1] = (fLastV-fFirstV + 1)*fGeom->WirePitch(geo::kV);

  }


  //----------------------------------------------------------
  void APAGeometryAlg::ChannelToAPA( uint32_t chan, 
				     unsigned int & apa, 
				     unsigned int & cryo){

    cryo  =  chan / (fAPAsPerCryo*fChannelsPerAPA);

    // Number apa uniquely across cryostats so that
    // apa to recob::Object maps are easy to work with.
    // If we decide to reset APA number per cryo, uncomment:
    //chan -= cryo*fAPAsPerCryo*fChannelsPerAPA; 
    apa   =  chan / fChannelsPerAPA;

    return;
  }

  //----------------------------------------------------------
  unsigned int APAGeometryAlg::ChannelToAPA( uint32_t chan ){

    return chan / fChannelsPerAPA;
  }


  //----------------------------------------------------------
  unsigned int APAGeometryAlg::ChannelsInView( geo::View_t geoview ){

    switch(geoview){
    default       :  
      return 0;
    case geo::kU  :
      return ChannelsInAPAView( kU );
    case geo::kV  :
      return ChannelsInAPAView( kV );
    case geo::kZ  :
      if( ChannelsInAPAView( kZ0 ) == ChannelsInAPAView( kZ1 ) ) 
	return ChannelsInAPAView( kZ0 );
      else throw cet::exception("ChannelsInView")
	     << "Both Z sides should have the same amount of channels\n";
    }

  }


  //----------------------------------------------------------
  unsigned int APAGeometryAlg::ChannelsInAPAView( APAView_t apaview ){

    switch(apaview){
    default  :  
      return 0;
    case kU  :
      return fFirstV-fFirstU;
    case kV  :
      return fFirstZ0-fFirstV;
    case kZ0 :
      return fFirstZ1-fFirstZ0;
    case kZ1 :
      return fLastZ1-fFirstZ1 + 1;
    }

  }


  //----------------------------------------------------------
  uint32_t APAGeometryAlg::FirstChannelInView( geo::View_t geoview, 
					       unsigned int apa, 
					       unsigned int cryo ){


    switch(geoview){
    default       :  
      return 0        + (uint32_t)(apa + cryo*fAPAsPerCryo)*fChannelsPerAPA;
    case geo::kU  :
      return fFirstU  + (uint32_t)(apa + cryo*fAPAsPerCryo)*fChannelsPerAPA;
    case geo::kV  :
      return fFirstV  + (uint32_t)(apa + cryo*fAPAsPerCryo)*fChannelsPerAPA;
    case geo::kZ  :
      //TODO: would need tpc number for the rest of this
      return fFirstZ0 + (uint32_t)(apa + cryo*fAPAsPerCryo)*fChannelsPerAPA;
    }

  }

  //----------------------------------------------------------
  uint32_t APAGeometryAlg::FirstChannelInView( uint32_t chan ){

    geo::View_t geoview = fGeom->View(chan);
    unsigned int apa, cryo;
    this->ChannelToAPA( chan, apa, cryo );
    return this->FirstChannelInView( geoview, chan );

  }

  //----------------------------------------------------------
  uint32_t APAGeometryAlg::FirstChannelInView( geo::View_t geoview, 
					       uint32_t chan ){

    unsigned int apa, cryo;
    this->ChannelToAPA( chan, apa, cryo );
    return this->FirstChannelInView( geoview, apa, cryo );

  }



  //----------------------------------------------------------
  APAView_t APAGeometryAlg::APAView( uint32_t chan ){
    
    // it seems trivial to do this for U and V, but this gives a side to
    // geo::kZ, unlike Geometry::View(c), as is often needed in disambiguation

    geo::View_t view = fGeom->View( chan );
    switch(view){
    default       :  
      break;
    case geo::kU  :
      return kU;
    case geo::kV  :
      return kV;
    case geo::kZ  :
      unsigned int modchan = chan % fChannelsPerAPA;
      // Channel mapping number in the order of U, V, Z0, then Z1
      if( modchan > fLastV && modchan < fFirstZ1 ) return kZ0;
      if( modchan > fLastZ0  ) return kZ1;
    }

    return kUnknown;

  }


  //----------------------------------------------------------
  std::vector<geo::WireID> APAGeometryAlg::ChanSegsPerSide(uint32_t chan, unsigned int side){

    std::vector<geo::WireID> wids = fGeom->ChannelToWire(chan);
    return this->ChanSegsPerSide(wids, side);

  }



  //----------------------------------------------------------
  std::vector<geo::WireID> APAGeometryAlg::ChanSegsPerSide(std::vector<geo::WireID> wids, unsigned int side)
  {
    // Given a vector of wireIDs and an APA side, return
    // the wireIDs the the tpc side where tpc%2 = side
    
    std::vector<geo::WireID> thisSide;
    
    for(size_t i = 0; i < wids.size(); i++)
      if( wids[i].TPC % 2 == side ) thisSide.push_back(wids[i]);
    
    return thisSide;
  }




  //----------------------------------------------------------
  geo::WireID APAGeometryAlg::NearestWireIDOnChan(  const double WorldLoc[3],  
						    uint32_t chan,
						    unsigned int const plane,  
						    unsigned int const tpc,    
						    unsigned int const cstat   )
  {

    std::vector<geo::WireID> cWids = fGeom->ChannelToWire( chan );

    if( cWids[0].Cryostat != cstat ) 
      throw cet::exception("APAGeometryAlg") << "Channel " << chan 
					     << "not in cryostat " << cstat << "\n";
    if( std::floor( cWids[0].TPC / 2 ) != std::floor( tpc / 2 ) ) 
      throw cet::exception("APAGeometryAlg") << "Channel " << chan 
					     << "not in APA " << std::floor(tpc/2) << "\n";

    // special case for vertical wires
    if(fGeom->View(chan)==geo::kZ) return fGeom->ChannelToWire(chan)[0];

    unsigned int xyzWire = fGeom->NearestWireID( WorldLoc, plane, tpc, cstat ).Wire;

    // The desired wire ID will be the only channel 
    // segment within half the channel range.
    geo::WireID wid;
    for(size_t i=0; i<cWids.size(); i++){
      if( cWids[i].TPC != tpc ) continue;
      if( std::abs((int)cWids[i].Wire - (int)xyzWire) < fChannelRange[plane]/2 ) wid=cWids[i];
    }
    
    return wid;

  }


  //----------------------------------------------------------
  bool APAGeometryAlg::LineSegChanIntersect( TVector3 xyzStart, TVector3 xyzEnd, uint32_t chan, 
					     std::vector<geo::WireID> & widsCrossed,
					     bool ExtendLine = true )
  {

    // This assumes a smooth wire numbering, and that the line seg is contained in a tpc.
    // Meant for use with the approximate line calculated
    // by matching cluster endpoints in disambiguation.

    // Find tpc, use midpoint in case start/end is on a boundary
    unsigned int tpc, cryo;
    double xyzMid[3];
    xyzMid[0] = (xyzStart[0]+xyzEnd[0])/2;
    xyzMid[1] = (xyzStart[1]+xyzEnd[1])/2;
    xyzMid[2] = (xyzStart[2]+xyzEnd[2])/2;
    fGeom->PositionToTPC(xyzMid, tpc, cryo);

    // Find the nearest wire number to the line segment endpoints
    std::vector<geo::WireID> wids = fGeom->ChannelToWire(chan);
    unsigned int startW = fGeom->NearestWire( xyzStart, wids[0].Plane, tpc, cryo );
    unsigned int endW   = fGeom->NearestWire( xyzEnd,   wids[0].Plane, tpc, cryo );

    if( startW > endW ) std::swap(startW, endW);


    // Loop through wireIDs and check for intersection, if in the right TPC
    for( size_t w = 0; w < wids.size(); w++ ){
      if( wids[w].TPC != tpc ) continue;
      if( wids[w].Cryostat != cryo ) throw cet::exception("LineSegChanIntersect")
				       << "Channel and line not in the same crostat.\n";

      // If the current wire id wire number is inbetween the start/end 
      // point wires, the line segment intersects the wireID at some point.

      // TODO: for now, extend range, but that is application specific. fix asap
      // The longer we make the range, the more conservative it is, so it is safe
      // to extend the range a bit to get hits at the ends of the line
      unsigned int ext = 0;
      if ( ExtendLine) ext = 10;

      if( fGeom->ValueInRange( wids[w].Wire*1., (startW-ext)*1., (endW+ext)*1. ) ) widsCrossed.push_back(wids[w]);

    }

    if( widsCrossed.size() > 0 ) return true;
    else return false;

  }



  //----------------------------------------------------------
  std::vector<double>  APAGeometryAlg::ThreeChanPos( uint32_t u, uint32_t v, uint32_t z )
  {

    // Say we've associated a U, V, and Z channel -- perhaps by associating hits
    // or cluster endpoints -- these don't necessarily all intersect, but they
    // are hopefully pretty close. Find the center of the 3 intersections.

    // get data needed along the way
    std::vector< geo::WireIDIntersection > UVIntersects;
    this->APAChannelsIntersect( u, v, UVIntersects );
    std::vector< double > UVzToZ(UVIntersects.size());
    geo::WireID Zwid = fGeom->ChannelToWire(z)[0];
    unsigned int cryo = Zwid.Cryostat;
    unsigned int tpc = Zwid.TPC;
    std::vector<geo::WireID> Uwids = fGeom->ChannelToWire(u);
    std::vector<geo::WireID> Vwids = fGeom->ChannelToWire(v);
    std::vector<geo::WireID> UwidsInTPC, VwidsInTPC;
    for(size_t i=0; i<Uwids.size(); i++) if( Uwids[i].TPC==tpc ) UwidsInTPC.push_back(Uwids[i]);
    for(size_t i=0; i<Vwids.size(); i++) if( Vwids[i].TPC==tpc ) VwidsInTPC.push_back(Vwids[i]);
    double Zcent[3] = {0.};
    fGeom->WireIDToWireGeo( Zwid ).GetCenter(Zcent);

    std::cout << "Zcent = " << Zcent[2] << ", UVintersects zpos = ";
    for(size_t uv=0; uv<UVIntersects.size(); uv++){
      std::cout << UVIntersects[uv].z << ", ";
    }
    std::cout << "\n";

    /////////////////////////////////////
    /////////////////////////////////////

    //std::cout << "U = " << u << " V = " << v << ", " << UVIntersects.size() << std::endl;

    if( UVIntersects.size() == 0 ){
      if( UwidsInTPC.size() > 1 || VwidsInTPC.size() > 1 ) 
	throw cet::exception("ThreeChanPos")    << "U/V channels don't intersect, bad return.\n";

      // Now assume there are only one of each u and v wireIDs on in this TPC
      mf::LogWarning("ThreeChanPos") << "No U/V intersect, exceptional channels. See if U or V intersects Z\n";
      std::vector<double> yzCenter(2,0.);
      geo::WireID Uwid = UwidsInTPC[0];
      geo::WireID Vwid = VwidsInTPC[0];
      geo::WireIDIntersection UZInt, VZInt;
      bool checkUZ = fGeom->WireIDsIntersect( Uwid, Zwid, UZInt );
      bool checkVZ = fGeom->WireIDsIntersect( Vwid, Zwid, VZInt );
      if( !checkUZ && !checkVZ )
	throw cet::exception("NoChanIntersect") << "No channels intersect, bad return.\n";
      if( checkUZ && !checkVZ ){ yzCenter[0] = UZInt.y; yzCenter[1] = UZInt.z; }
      if( checkVZ && !checkUZ ){ yzCenter[0] = VZInt.y; yzCenter[1] = VZInt.z; }
      if( checkUZ &&  checkVZ ){ yzCenter[0] = (VZInt.y+UZInt.y)/2; yzCenter[1] = (VZInt.z+UZInt.z)/2; }
      return yzCenter;
    }

    /////////////////////////////////////
    /////////////////////////////////////


    // In case the uv channels intersect twice on the same side, choose the best case.
    // Note: this will not happen for APAs with UV angle at about 36, but will for 45
    std::cout << "UVzToZ = ";
    for( size_t widI = 0; widI < UVIntersects.size(); widI++ ){
      UVzToZ[widI] = std::abs( UVIntersects[widI].z - Zcent[2] );
      std::cout << UVzToZ[widI] << ", ";
    }
    std::cout<<"\n";

    unsigned int bestWidI = 0;
    double minZdiff = fGeom->Cryostat(cryo).TPC(tpc).Length(); // start it out at maximum z
    for( unsigned int widI = 0; widI < UVIntersects.size(); widI++ ){

      //std::cout << "widI = " << widI << std::endl; 

      if( UVIntersects[widI].TPC == tpc  &&  UVzToZ[widI] < minZdiff ){ 
	minZdiff = UVzToZ[widI];
	bestWidI = widI;
	//std::cout << "bestWidI = " << bestWidI << std::endl;
      }
    }
    geo::WireIDIntersection ChosenUVInt = UVIntersects[bestWidI];
    
    // Now having the UV intersection, get the UZ and VZ
    double UVInt[3] = {0.};
    UVInt[1] = ChosenUVInt.y; UVInt[2] = ChosenUVInt.z;
    geo::WireID Uwid = this->NearestWireIDOnChan( UVInt, u, 0, tpc, cryo );
    geo::WireID Vwid = this->NearestWireIDOnChan( UVInt, v, 1, tpc, cryo );
    geo::WireIDIntersection UZInt, VZInt;
    bool checkUZ = fGeom->WireIDsIntersect( Uwid, Zwid, UZInt );
    bool checkVZ = fGeom->WireIDsIntersect( Vwid, Zwid, VZInt );

    std::cout << "UZint.z = " << UZInt.z << " (" << checkUZ << "), VZint.z = " << VZInt.z << " (" << checkVZ << ")\n";
    
    // find the center
    std::vector<double> yzCenter(2,0.);


    if( !checkUZ || !checkVZ ){
      //throw cet::exception("ThreeChanPos")  << "WireIDs were expected to intersect.\n";

      std::cout << "ChosenUVint.y = " << ChosenUVInt.y << "ChosenUVint.z = " << ChosenUVInt.z << std::endl;
      
      //temporary case
      yzCenter[0] = ChosenUVInt.y;
      yzCenter[1] = ChosenUVInt.z;


    } else {

    yzCenter[0] =  (ChosenUVInt.y + UZInt.y + VZInt.y)/3; 
    yzCenter[1] =  (ChosenUVInt.z + UZInt.z + VZInt.z)/3;

    }


    return yzCenter;

  } 




  //----------------------------------------------------------
  bool APAGeometryAlg::APAChannelsIntersect( uint32_t chan1, uint32_t chan2, 
					     std::vector< geo::WireIDIntersection >& IntersectVector )
  {

    
    // Get the WireIDs and view for each channel, make sure views are different
    geo::WireIDIntersection    widIntersect;
    std::vector< geo::WireID > wids1 = fGeom->ChannelToWire( chan1 );
    std::vector< geo::WireID > wids2 = fGeom->ChannelToWire( chan2 );
    geo::View_t                view1 = fGeom->View( chan1 );
    geo::View_t                view2 = fGeom->View( chan2 );
    if( view1 == view2 ){
      mf::LogWarning("APAChannelsIntersect") << "Comparing two channels in the same view, return false";
      return false;     }
    if( wids1[0].Cryostat         != wids2[0].Cryostat         ||
	this->ChannelToAPA(chan1) != this->ChannelToAPA(chan2)    ){
      throw cet::exception("APAChannelsIntersect") << "Comparing two channels in in different APAs: "
						   << "channel " << chan1 << " in Cryo " 
						   << wids1[0].Cryostat << ", APA " << this->ChannelToAPA(chan1)
						   << ", and channel " << chan2 << " in Cryo " 
						   << wids2[0].Cryostat << ", APA " << this->ChannelToAPA(chan2)
						   << "\n";
      return false;     }
 

    // Loop through wids1 and see if wids2 has any intersecting wires, 
    // given that the WireIDs are in the same TPC
    for( unsigned int i1 = 0; i1 < wids1.size() ; i1++){
      for( unsigned int i2 = 0; i2 < wids2.size() ; i2++){

	// make sure it is reasonable to intersect
	if( wids1[i1].Plane    == wids2[i2].Plane ||
	    wids1[i1].TPC      != wids2[i2].TPC   ||      // not reasonable for a *WireID*
	    wids1[i1].Cryostat != wids2[i2].Cryostat ) continue;

// 	std::cout << "Checking: \n WireID 1 = (" 
// 		  << wids1[i1].Cryostat << "," << wids1[i1].TPC << ","
// 		  << wids1[i1].Plane    << "," << wids1[i1].Wire
// 		  << ") \n WireID 2 = ("
// 		  << wids2[i2].Cryostat << "," << wids2[i2].TPC << ","
// 		  << wids2[i2].Plane    << "," << wids2[i2].Wire << ")" << std::endl;

	// Check if they even intersect; if they do, push back
	if( fGeom->WireIDsIntersect( wids1[i1], wids2[i2], widIntersect ) ){

//	  std::cout << "we have an intersect" << std::endl;

	  IntersectVector.push_back( widIntersect );
	}
      }
    }

    // Of all considered configurations, there are never more than 
    // 4 intersections per channel pair
    if( IntersectVector.size() > 4 ){
      mf::LogWarning("APAChannelsIntersect") << "Got " << IntersectVector.size() 
					     << " intersections for channels "
					     << chan1 << " and " << chan2 
					     << " - never expect more than 4, so far";  }


    // With increasing IntersectVector index, the WireID 
    // vector indices of the intersecting wireIDs increase.
    //   This matches the direction in which ChannelToWire
    //   builds its output WireID vector in the APA/35t Alg
    std::sort( IntersectVector.begin(), IntersectVector.end() );

    // return true if any intersection points were found
    if( IntersectVector.size() == 0 ) return false;
    else return true;

  }




} //end namespace apa
