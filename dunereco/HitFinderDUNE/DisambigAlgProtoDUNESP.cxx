////////////////////////////////////////////////////////////////////////
//
// DisambigAlgProtoDUNESP.cxx
//
// trj@fnal.gov
// tjyang@fnal.gov
//
// description
//
// An algorithm for disambiguation for the single-phase protoDUNE detector that is robust against missing the third wire plane.  
// Assumes hits cannot be found outside the outer APA's. Assumes that the volumes outside of the APA's are in fact TPC's.
// Assumes 36-degree DUNE geometry with no more than one intersection of an induction-plane wire with a collection-plane wire.
// Since the U and V angles are the same, assume that the a U wire intersects a V wire in at most one place as well.
// Assumes there are two induction planes. kU, kV, and one collection plane kZ.
//
//
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "DisambigAlgProtoDUNESP.h"
#include "larevt/Filters/ChannelFilter.h"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TH1D.h"

namespace dune{

  DisambigAlgProtoDUNESP::DisambigAlgProtoDUNESP(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset); 
  }

  //----------------------------------------------------------
  void DisambigAlgProtoDUNESP::reconfigure(fhicl::ParameterSet const& p)
  {

    fTimeCut = p.get<double>("TimeCut");
    fDistanceCut = p.get<double>("DistanceCut");
    fDcut2 = fDistanceCut * fDistanceCut;
  }


  //----------------------------------------------------------
  //----------------------------------------------------------
  void DisambigAlgProtoDUNESP::RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits   )
  {
    fDisambigHits.clear();

    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<geo::Geometry> geo;

    size_t napas = geo->NTPC()/2;

    // do only one APA at a time.

    for (size_t apa=0; apa<napas; apa++)
      {

	int tpc = longTPC(apa); // for purposes of evaluating time offsets for time matching.
	//tick offsets differ for the even and odd TPC's
        int cryostat = 0; // protoDUNE-SP has only one cryostat

	//  TH1D *histu = new TH1D("histu","histu",4000,0,4000);
	//  TH1D *histv = new TH1D("histv","histv",4000,0,4000);
	//  TH1D *histz = new TH1D("histz","histz",4000,0,4000);

	// copy hit ptrs to local storage

	std::vector<art::Ptr<recob::Hit> > hitsUV[2];  // index 0=U, 1=V
	std::vector<art::Ptr<recob::Hit> > hitsZ;

	//std::cout << "Dumping hits for this APA" << std::endl;

	for (size_t i = 0; i<OrigHits.size(); ++i)
	  {
	    unsigned int hitapa(0), hitcryo(0);
	    fAPAGeo.ChannelToAPA(OrigHits[i]->Channel(), hitapa, hitcryo);
	    if (apa == hitapa)
	      {
		//std::cout << *OrigHits[i] << std::endl;

		switch (OrigHits[i]->View())
		  {
		  case geo::kU:
		    hitsUV[0].push_back(OrigHits[i]);
		    //      if (OrigHits[i]->WireID().TPC==1) 
		    //	histu->Fill(OrigHits[i]->PeakTime()
		    //		    - detprop->GetXTicksOffset(0,
		    //					       tpc,
		    //					       cryostat)
		    //		    ,OrigHits[i]->Charge());

		    break;
		  case geo::kV:
		    hitsUV[1].push_back(OrigHits[i]);
		    //      if (OrigHits[i]->WireID().TPC==1) 
		    //	histv->Fill(OrigHits[i]->PeakTime()
		    //		    - detprop->GetXTicksOffset(1,
		    //					       tpc,
		    //					       cryostat)
		    //		    ,OrigHits[i]->Charge());
		    break;
		  case geo::kZ:
		    hitsZ.push_back(OrigHits[i]);
		    //      if (OrigHits[i]->WireID().TPC==1) 
		    //	histz->Fill(OrigHits[i]->PeakTime()
		    //		    - detprop->GetXTicksOffset(2,
		    //					       tpc,
		    //					       cryostat)
		    //		    ,OrigHits[i]->Charge());
		    break;
		  default:
		    throw cet::exception("DisambigAlgProtoDUNESP") <<": hit view unkonwn. \n";
		  }
	      }
	  }

	//  std::cout << " DisambigAlgProtoDUNESP timing means: " <<histu->GetMean()<<" "<<histv->GetMean()<<" "<<histz->GetMean()<<std::endl;
	//  delete histu;
	//  delete histv;
	//  delete histz;

	size_t zsize = hitsZ.size();

	// loop over U and V planes to disambiguate.
	// Identify hits matching in time in the Z and the other induction planes

	for (size_t uv=0; uv<2; uv++)  // uv is the current induction plane
	  {
	    size_t other=1-uv;  // the index of the other induction plane
	    size_t uvsize = hitsUV[uv].size();    // this induction plane's hits
	    size_t othersize = hitsUV[other].size();  // the other induction plane's hits

	    for (size_t iuv = 0; iuv < uvsize; ++iuv)
	      {
		//std::cout << "Attempting to Disambiguate hit: " << *hitsUV[uv][iuv] << std::endl;
		//std::cout << "Getting time offset: " << hitsUV[uv][iuv]->WireID().Plane << " " << tpc << " " << cryostat << std::endl;
		//std::cout << "time offset: " << detprop->GetXTicksOffset(hitsUV[uv][iuv]->WireID().Plane,tpc,cryostat) << std::endl;

		double tuv = hitsUV[uv][iuv]->PeakTime()
		  - detprop->GetXTicksOffset(hitsUV[uv][iuv]->WireID().Plane,tpc,cryostat);

		std::vector<size_t> zmatches;  // list of indices of time-matched hits
		std::vector<size_t> othermatches;
		for (size_t z = 0; z < zsize; z++)
		  {
		    double tz = hitsZ[z]->PeakTime()
		      - detprop->GetXTicksOffset(hitsZ[z]->WireID().Plane,tpc,cryostat);
		    //std::cout << "Looking for a Z match: " << tuv << " " << tz << " Offset: " <<
                    // detprop->GetXTicksOffset(hitsZ[z]->WireID().Plane,tpc,cryostat)
		    //      << " plane: " << hitsZ[z]->WireID().Plane << " TPC: " << tpc << " Cryo: " << cryostat
		    //	      << std::endl;

		    if (std::abs(tuv-tz)<fTimeCut)
		      {
			zmatches.push_back(z);
		      }
		  }
		for (size_t iother = 0; iother < othersize; iother++)
		  {
		    double tother = hitsUV[other][iother]->PeakTime()
		      - detprop->GetXTicksOffset(hitsUV[other][iother]->WireID().Plane,tpc,cryostat);
		    if (std::abs(tuv-tother)<fTimeCut)
		      {
			othermatches.push_back(iother);
		      }
		  }

		//std::cout << "Number of time-matched Z hits: " << zmatches.size() << std::endl;
		//std::cout << "Number of time-matched Other-Ind View hits: " << othermatches.size() << std::endl;

		if (zmatches.size() == 0 && othermatches.size() == 0) continue;   // hit has no matches -- possibly noise

		// find out how many of these matched-in-time hits are also matched in space.
		// take triplets over doublets.  Assume we are on the sensitive side of the APA.
	  
		std::vector< double > zmatchz;
		std::vector< double > zmatchy;
		std::vector< double > othermatchz;
		std::vector< double > othermatchy;

		std::vector<geo::WireID>  wires = geo->ChannelToWire(hitsUV[uv][iuv]->Channel());
		size_t wsize = wires.size();
		std::vector<size_t> ndoublets(wsize,0);
		std::vector<size_t> ntriplets(wsize,0);

		for (size_t w=0; w<wsize; w++)
		  {
		    ndoublets[w] = 0;
		    ntriplets[w] = 0;
		    geo::WireID uvwire = wires[w];
		    if ( notOuterWire(uvwire) )
		      {
			for (size_t z=0; z<zmatches.size(); z++)
			  {
			    geo::WireID zwire = geo->ChannelToWire(hitsZ[zmatches[z]]->Channel())[0];
			    if ( notOuterWire(zwire) )  // we really shouldn't have any hits on the outer z wires
			      {
			        geo::WireIDIntersection isect;
			        if (geo->WireIDsIntersect(zwire,uvwire,isect))
			          {
				    zmatchz.push_back(isect.z);
				    zmatchy.push_back(isect.y);
			          }
			      }
			  }
			// There is at most one intersection of the u channel with a v channel. Still have to loop over possibilities though.
			for (size_t iother=0; iother<othermatches.size(); iother++)
			  {
			    std::vector<geo::WireID>  otherwires = geo->ChannelToWire(hitsUV[other][othermatches[iother]]->Channel());
			    for (size_t otherw = 0; otherw < otherwires.size(); otherw++)
			      {
				geo::WireID otherwire = otherwires[otherw];
				if (notOuterWire(otherwire))
				  {
				    geo::WireIDIntersection isect;
				    if (geo->WireIDsIntersect(otherwire,uvwire,isect))
				      {
					othermatchz.push_back(isect.z);
					othermatchy.push_back(isect.y);
				      }
				  }
			      }
			  }

			// look for triplets among the matches, and then doublets
			for (size_t izmatch=0;izmatch<zmatchz.size();izmatch++)
			  {
			    bool found_triplet = false;
			    for (size_t iothermatch=0;iothermatch<othermatchz.size();iothermatch++)
			      {
				double disz = zmatchz[izmatch]-othermatchz[iothermatch]; 
				double disy = zmatchy[izmatch]-othermatchy[iothermatch]; 
				if ( disz*disz + disy*disy < fDcut2 )
				  {
				    ntriplets[w]++;
				    found_triplet = true;
				  }
			      }
			    if (!found_triplet) ndoublets[w] ++;
			  }
			for (size_t iothermatch=0;iothermatch<othermatchz.size();iothermatch++)
			  {
			    bool found_triplet = false;
			    for (size_t izmatch=0;izmatch<zmatchz.size();izmatch++)
			      {
				double disz = zmatchz[izmatch]-othermatchz[iothermatch]; 
				double disy = zmatchy[izmatch]-othermatchy[iothermatch]; 
				if ( disz*disz + disy*disy < fDcut2 )
				  {
				    found_triplet = true;
				  }
			      }
			    if (!found_triplet) ndoublets[w] ++;
			  }  // end summing up doublets and triplets
		      }  // end check if not an outer wire
		  } // end loop on wire candidates for U channels

		bool found_disambig = false;
		size_t bestwire = 0;
		for (size_t w=0; w<wsize; w++)
		  {
		    if (!found_disambig && (ntriplets[w]>0 || ndoublets[w] > 0))
		      {
			bestwire = w;
			found_disambig = true;
		      }
		    if (ntriplets[w] > ntriplets[bestwire]) 
		      {
			bestwire = w;
			found_disambig = true;
		      }
		    if ( (ntriplets[w] == ntriplets[bestwire]) && (ndoublets[w] > ndoublets[bestwire])) 
		      {
			bestwire = w;
			found_disambig = true;
		      }
		  }
		if (found_disambig)
		  {
		    fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[uv][iuv],wires[bestwire]));
		  }
		else
		  {
		    //mf::LogWarning("DisambigAlg35t")<<"Could not find disambiguated hit for  "<<*hitsUV[uv][iuv]<<"\n";
		    //std::cout <<"Could not find disambiguated hit for  "<<*hitsUV[uv][iuv]<<"\n";
		  }

	      } // end loop over induction hits in this plane
	  } // end loop over induction planes 
      } // end loop over APA's
  }

// method to tell us whether a particular wire is in an outer TPC.  Hardcoded outer TPC numbers.
// might be quicker to check the lowest two bits: 00 or 11 mean outer TPC.

  bool DisambigAlgProtoDUNESP::notOuterWire(geo::WireID wireid)
  {
    auto tpc = wireid.TPC;
    bool result;
    if (tpc == 0 || tpc == 3 || tpc == 4 || tpc == 7 || tpc == 8 || tpc == 11)
      {  
        result = false;
      }
    else
      {
        result = true;
      }
    return(result);
  }

  int DisambigAlgProtoDUNESP::longTPC(int apa)
  {
    int tpc=0;

    if (apa == 0) tpc=1;
    if (apa == 1) tpc=2;
    if (apa == 2) tpc=5;
    if (apa == 3) tpc=6;
    if (apa == 4) tpc=9;
    if (apa == 5) tpc=10;

    return(tpc);
  }
}

// todo -- look at hits nearby undisambiguated hits and assign similar disambiguation assignments
