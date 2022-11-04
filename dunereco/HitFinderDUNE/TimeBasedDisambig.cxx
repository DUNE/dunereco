////////////////////////////////////////////////////////////////////////
//
// TimeBasedDisambig.cxx
//
// trj@fnal.gov
// tjyang@fnal.gov
// iamjaejunkim@gmail.com
//
// description
//
// Based on time based hit matching algorithm
//
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "TimeBasedDisambig.h"
//#include "MCCheater/BackTracker.h"
#include "larevt/Filters/ChannelFilter.h"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TH1D.h"
#include "TF1.h"

using namespace std;

namespace dune{

TimeBasedDisambig::TimeBasedDisambig(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset); 
}

//----------------------------------------------------------
void TimeBasedDisambig::reconfigure(fhicl::ParameterSet const& p)
{
  
  fTimeCut = p.get<double>("TimeCut");
  fDistanceCut = p.get<double>("DistanceCut");
  fDistanceCutClu = p.get<double>("DistanceCutClu");
}


//----------------------------------------------------------
//----------------------------------------------------------
void TimeBasedDisambig::RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits )
//void TimeBasedDisambig::RunDisambig()
{
  //fDisambigHits.clear();

  //create geometry and backtracker servicehandle object
  art::ServiceHandle<geo::Geometry> geo;
  //art::ServiceHandle<cheat::BackTracker> bt;
  //define timeoffset and the maximum chargeratio allowed when matching induction plane hits
  //the timeoffset array is going to be replaced with defined values in the detectorproperties package later
  double timeoffset[3]={10.4,5.2,0.0};
  double rmax=2.5;

  //create vectors for induction plane and collection plane hits
  std::vector<art::Ptr<recob::Hit> > hitsUV;
  std::vector<art::Ptr<recob::Hit> > hitsZ;

  //create a vector where the disambiguated hits are going to be stored
  //std::vector<art::Ptr<recob::Hit> > DisambiguatedHits;
  //std::unique_ptr<std::vector<recob::Hit> > DisambiguatedHits(new std::vector<recob::Hit>);
  //std::vector< <recob::Hit> > DisambiguatedHits;
  std::vector< recob::Hit > DisambiguatedHits;

  //double HitsTotal=0;
  //double HitsCorrect=0;

  //loop over detector volumn in each component
  for (unsigned int Cstat=0; Cstat < geo->Ncryostats(); ++Cstat){
    for (unsigned int APA=0; APA < geo->Cryostat(Cstat).NTPC()/2; ++APA){
      hitsUV.clear();
      hitsZ.clear();
      //loop over all induction and collection plane hits and save the hits within each component
      for (size_t i = 0; i<OrigHits.size(); ++i){
	if (OrigHits[i]->WireID().Cryostat==Cstat && (OrigHits[i]->WireID().TPC==APA*2 || OrigHits[i]->WireID().TPC==APA*2+1)){	  
	  
	  //save induction and collection plane hits in two separate vectors
	  if (OrigHits[i]->View()!=geo::kZ){
	    hitsUV.push_back(OrigHits[i]);
	  } else {
	    hitsZ.push_back(OrigHits[i]);
	  }

	}
      }

      //define variables to save boundary information
      double origin[3]={0.0,0.0,0.0};
      double world[3]={0.0,0.0,0.0};
      //double xbound[2]={-2000,2000};
      //define xbound using peaktime for now but need to be replaced with coverted position in the future
      double xbound[2]={0,2000};
      double ybound[2]= {0.0,0.0};
      double zbound[2]= {0.0,0.0};
      geo->Cryostat(Cstat).TPC(APA*2).LocalToWorld(origin, world);
      ybound[0]=world[1]-geo->Cryostat(Cstat).TPC(APA*2).HalfHeight();
      ybound[1]=world[1]+geo->Cryostat(Cstat).TPC(APA*2).HalfHeight();
      zbound[0]=world[2]-0.5*geo->Cryostat(Cstat).TPC(APA*2).Length();
      zbound[1]=world[2]+0.5*geo->Cryostat(Cstat).TPC(APA*2).Length();

      //position-correction function needs to be applied on hits in certain time and longitudinal position ranges
      //thus need to create a struct of vectors where hit informations including their positions returned in the
      //disambiguation algorithm is going to be saved.
      HitPos h;
      std::vector<HitPos> InductionAll;
      InductionAll.clear();


      //loop over all induction plane hits
      for (unsigned int uv0=0; uv0 < hitsUV.size(); ++uv0){
	double PeakTimeMinusUV=hitsUV[uv0]->PeakTimePlusRMS(-1.)+timeoffset[1];
	double PeakTimeUV=hitsUV[uv0]->PeakTime()+timeoffset[1];
	if (hitsUV[uv0]->View()==geo::kU){
	  PeakTimeMinusUV+=timeoffset[1];
	  PeakTimeUV+=timeoffset[1];
	}

	//define a varible to save charge per num
	double ChargePerNum0=hitsUV[uv0]->Integral()/hitsUV[uv0]->Multiplicity();
	
	//loop over induction plane hits on a different view and among the induction plane hits whose charge per
	//num to that of the hit in loop is within 2.5 find an induction plane hit whose peaktime is closest to
	//that is in loop.
	double MeanPeakTimeUV=99999;
	unsigned int ChannelUV1=0;
	for (unsigned int uv1=0; uv1 < hitsUV.size(); ++uv1){
	  double PeakTimeMinusUV1=hitsUV[uv1]->PeakTimePlusRMS(-1.)+timeoffset[1];
	  double PeakTimeUV1=hitsUV[uv1]->PeakTime()+timeoffset[1];
	  if (hitsUV[uv1]->View()==geo::kU){
	    PeakTimeMinusUV1+=timeoffset[1];
	    PeakTimeUV1+=timeoffset[1];
	  }
	  double rChargeUV=hitsUV[uv0]->Integral()/hitsUV[uv1]->Integral();
	  double ChargePerNum1=hitsUV[uv1]->Integral()/hitsUV[uv1]->Multiplicity();
	  if (hitsUV[uv0]->Multiplicity()>0 && hitsUV[uv1]->Integral()>0) rChargeUV=ChargePerNum0/ChargePerNum1;

	  if (hitsUV[uv0]->View()!=hitsUV[uv1]->View() &&
	      fabs(PeakTimeMinusUV-PeakTimeMinusUV1)<MeanPeakTimeUV &&
	      fabs(PeakTimeUV-PeakTimeUV1)<20 &&
	      (rChargeUV<rmax && rChargeUV>(1/rmax)) ){
	    MeanPeakTimeUV=fabs(PeakTimeMinusUV-PeakTimeMinusUV1);
	    ChannelUV1=hitsUV[uv1]->Channel();
	  }
	}
	
	//define variables to save positions returned
	unsigned int ChannelUV0=hitsUV[uv0]->Channel();
       	double CandPosition[3]={0.0,0.0,0.0}; //candidate positions going to be returned in intersectionpoint function
	double FinalPosition[3]={0.0,0.0,0.0}; //final position among the candidate positions
	unsigned int FinalWireID=0;//final wireid for the induction plane hit in the loop
	double MinVertBetUVandZ=99999;
	//unsigned int FinalWire=0;

	//read wireid vectors of induction plane hits matched in the above loop
	std::vector<geo::WireID> wireiduv0 = geo->ChannelToWire(ChannelUV0);
	std::vector<geo::WireID> wireiduv1 = geo->ChannelToWire(ChannelUV1);

	//loop over the induction plane hit wireids, find candidate intersection positions using intersectionpoint
	//function and find a position whose longitudinal position is closest to the average position returned above
	//and save the positions and wireid in variables defined above
	for (unsigned int wireid0=0; wireid0 < wireiduv0.size(); ++wireid0){
	  double VertPos=0; //vertical position of collection plane wire
	  double MinPeakTime=99999; //variable to store peaktime information temporarily
	  //loop over collection plane hits and among the collection plane hits whose charge per num to that 
	  //of the hit in loop is within 2.5 find a collection plane hit that has closest time to the induction plane
	  //hit in the loop.
	  for (unsigned int z=0; z < hitsZ.size(); ++z){
	    double PeakTimeMinusZ=hitsZ[z]->PeakTimePlusRMS(-1.);
	    double PeakTimeZ=hitsZ[z]->PeakTime();
	    double ChargePerNumZ=hitsZ[z]->Integral()/hitsZ[z]->Multiplicity();
	    double rChargeZ=hitsUV[uv0]->Integral()/hitsZ[z]->Integral();
	    if (hitsUV[uv0]->Multiplicity()>0 && hitsZ[z]->Integral()>0) rChargeZ=ChargePerNum0/ChargePerNumZ;
	    std::vector<geo::WireID> wires = geo->ChannelToWire(hitsZ[z]->Channel());
            auto const vert = geo->Wire(wires[0]).GetCenter();
	    if (fabs(PeakTimeMinusUV-PeakTimeMinusZ)<MinPeakTime && fabs(PeakTimeUV-PeakTimeZ)<20 &&
		(rChargeZ<rmax && rChargeZ>(1/rmax)) &&
	  	wireiduv0[wireid0].TPC==wires[0].TPC){
	      MinPeakTime=fabs(PeakTimeMinusUV-PeakTimeMinusZ);
              VertPos=vert.Z();
	    } 
	  }

	  //loop over collection plane hits and among the collection plane hits whose charge per num to that
	  //of the hit in the loop is within 2.5 find collection plane hits whose time difference to the induction
	  //plane hit is within 10 and take the average position of the collection plane hits.
	  double VertPosMean=0; //average position of all collection plane hits whose peaktime are in the window
	  double RMS=0; //rms of the collection plane hit positions 
	  unsigned int CollectionHits=0; //total no of collection plane hits in the time window
	  for (unsigned int z=0; z < hitsZ.size(); ++z){
	    double PeakTimeMinusZ=hitsZ[z]->PeakTimePlusRMS(-1.);
	    //double PeakTimeZ=hitsZ[z]->PeakTime();
	    double ChargePerNumZ=hitsZ[z]->Integral()/hitsZ[z]->Multiplicity();
	    double rChargeZ=hitsUV[uv0]->Integral()/hitsZ[z]->Integral();
	    if (hitsUV[uv0]->Multiplicity()>0 && hitsZ[z]->Integral()>0) rChargeZ=ChargePerNum0/ChargePerNumZ;
	    std::vector<geo::WireID> wires = geo->ChannelToWire(hitsZ[z]->Channel());
            auto const vert = geo->Wire(wires[0]).GetCenter();
	    if (fabs(fabs(PeakTimeMinusUV-PeakTimeMinusZ)-MinPeakTime)<10 &&
		//fabs(PeakTimeUV-PeakTimeZ)<20 &&
		(rChargeZ<rmax && rChargeZ>(1/rmax)) &&
		wireiduv0[wireid0].TPC==wires[0].TPC ){
	      CollectionHits++;
              VertPosMean+=vert.Z();
              RMS+=(VertPos-vert.Z())*(VertPos-vert.Z());
	    } 
	  }
	  VertPosMean=VertPosMean/double(CollectionHits); //average position
	  RMS=RMS/double(CollectionHits); //rms of the collection plane hits
	  
	  
	  //loop over the wireid of the induction plane hit whose time is closest to the hit in the first loop
	  //and find intersection and find a wire whose returned position is closest to the position of the collection
	  //plane hit found in the loop above.
	  for (unsigned int wireid1=0; wireid1 < wireiduv1.size(); ++wireid1){
	    if (abs(int(wireid0)-int(wireid1))<3 && wireiduv0[wireid0].TPC==wireiduv1[wireid1].TPC){
	      geo->IntersectionPoint(wireiduv0[wireid0].Wire, wireiduv1[wireid1].Wire,
				     wireiduv0[wireid0].Plane, wireiduv1[wireid1].Plane,
				     wireiduv0[wireid0].Cryostat, wireiduv0[wireid0].TPC,
				     CandPosition[1], CandPosition[2]);
	      //double VertBetUVandZ=fabs(VertPosMean-CandPosition[2]);
	      double VertBetUVandZ=fabs(VertPos-CandPosition[2]);
	      if (CandPosition[1]>=ybound[0] && CandPosition[1]<=ybound[1] && CandPosition[2]>=zbound[0] && CandPosition[2]<=zbound[1]){
		if (VertBetUVandZ<MinVertBetUVandZ){
		  MinVertBetUVandZ=VertBetUVandZ;
		  FinalWireID=wireid0;
		  //FinalWire=wireiduv0[wireid0].Wire;
		  FinalPosition[1]=CandPosition[1];
		  FinalPosition[2]=CandPosition[2];
		}
	      }
	    }
	  }
	}

	//read simulated position of the induction plane hit using backtracker
	// trj: skip this as it is part of the analysis and not the reco alg
	//std::vector<double> SimPosition;
	//unsigned int SimTPC, SimCstat, SimWire, SimWireID;
	//SimPosition=bt->HitToXYZ(hitsUV[uv0]);
	//if (SimPosition[2]!=99999){
	// double xyzpos[3];
	//xyzpos[0]=SimPosition[0];
	//xyzpos[1]=SimPosition[1];
	//xyzpos[2]=SimPosition[2];
	// geo->PositionToTPC(xyzpos, SimTPC, SimCstat);
	//SimWire=int(0.5+geo->WireCoordinate(xyzpos[1], xyzpos[2], hitsUV[uv0]->WireID().Plane, SimTPC, SimCstat));
	// for (unsigned int wireid0=0; wireid0 < wireiduv0.size(); ++wireid0){
	//   if (abs(int(wireiduv0[wireid0].Wire)-int(SimWire))<2) SimWireID=wireid0;
	// }
	//}

	//save timeoffset corrected peaktime
	double PeakTime=hitsUV[uv0]->PeakTime();
	if (hitsUV[uv0]->View()==geo::kU) PeakTime+=timeoffset[0];
	if (hitsUV[uv0]->View()==geo::kV) PeakTime+=timeoffset[1];

	//save hit information in a struct of a vector including positions returned in the disambiguation module
	//and simulated hit information
	h.Channel=hitsUV[uv0]->Channel();
	h.StartTick=0;
	h.EndTick=0;
	h.PeakTime=hitsUV[uv0]->PeakTime();
	h.SigmaPeakTime=hitsUV[uv0]->SigmaPeakTime();
	h.RMS=hitsUV[uv0]->RMS();
	h.PeakAmplitude=hitsUV[uv0]->PeakAmplitude();
	h.SigmaPeakAmplitude=hitsUV[uv0]->SigmaPeakAmplitude();
	h.SummedADC=hitsUV[uv0]->SummedADC();
	h.Integral=hitsUV[uv0]->Integral();
	h.SigmaIntegral=hitsUV[uv0]->SigmaIntegral();
	h.Multiplicity=hitsUV[uv0]->Multiplicity();
	h.LocalIndex=hitsUV[uv0]->LocalIndex();
	h.GoodnessOfFit=hitsUV[uv0]->GoodnessOfFit();
	h.DegreesOfFreedom=hitsUV[uv0]->DegreesOfFreedom();
	h.View=hitsUV[uv0]->View();
	h.SignalType=hitsUV[uv0]->SignalType();
	h.WireID=wireiduv0[FinalWireID];
	h.FinalYPos=FinalPosition[1];
	h.FinalZPos=FinalPosition[2];
	h.TPC=hitsUV[uv0]->WireID().TPC;
	h.DisambigWireID=FinalWireID;
	//h.SimYPos=SimPosition[1];
	//h.SimZPos=SimPosition[2];
	//h.SimWireID=SimWireID;
	InductionAll.push_back(h);
	
      }//end of induction plane hit loop


      //the disambiguation algorithm above returns a vector whose returned intersection position of two
      //induction plane hits and returns a wireid of a wire whose position is closest to the returned
      //position.  the function below is to correct the position of hits whose side is correct disambiguated
      //but the position is off.  in the function, the detector volumne is going to be divided in xz position
      //space and yposition of induction plane hits in the position space is going to be saved in a histogram,
      //find a bin whose bincontents is highest, calculate the mean position as fitting the histogram around
      //the bin, find a bin whose distance to the mean position found is about 100 and correct the position
      //and wireid of hits in the bin.
      int BinPerPos=4;
      unsigned int CellsPerTPC=BinPerPos*BinPerPos;//cells in xz position space
      double TotalXbin=(xbound[1]-xbound[0])/double(BinPerPos);
      double TotalZbin=(zbound[1]-zbound[0])/double(BinPerPos);
      //define histograms to save y position returned in the disambiguation algorithm above for all cells in
      //xz position space
      std::vector<TH1F *> YPosXZ(CellsPerTPC,0);
      for (unsigned int cell=0; cell<YPosXZ.size(); cell++){
	YPosXZ[cell]=new TH1F(Form("cell%d",cell),Form("cell%d",cell), 40, ybound[0], ybound[1]);
      }

      //loop over all induction plane hits and save the returned yposition in a histogram
      for (unsigned int uv0=0; uv0 < InductionAll.size(); ++uv0){
	int XbinForHit=int( (InductionAll[uv0].PeakTime-xbound[0])/TotalXbin );
	int ZbinForHit=int( (InductionAll[uv0].FinalZPos-zbound[0])/TotalZbin );
	if (InductionAll[uv0].PeakTime>xbound[1]) XbinForHit=BinPerPos-1;
	if (InductionAll[uv0].FinalZPos<zbound[0]){
	  ZbinForHit=0;
	} else if (InductionAll[uv0].FinalZPos>zbound[1]){
	  ZbinForHit=BinPerPos-1;
	}
	int BinForHit=XbinForHit*BinPerPos+ZbinForHit;
	
	YPosXZ[BinForHit]->Fill(InductionAll[uv0].FinalYPos);
      }//end of induction plane hit loop

      //define variables to save minimum and maximum range of the bin whose content is highest,
      //a variable to indicate whether the position-correction function should be applied or not
      //and variables to save minimum and maximum range of the bin whose mean position to the
      //position of the highestcontent bin is 100cm away.
      std::vector<double> FirstPeakXmin(CellsPerTPC,0);
      std::vector<double> FirstPeakXmax(CellsPerTPC,0);
      std::vector<double> ShouldCorrect(CellsPerTPC,0);
      std::vector<double> SecondPeakXmin(CellsPerTPC,0);
      std::vector<double> SecondPeakXmax(CellsPerTPC,0);

      //loop over all histograms
      for (unsigned int cell=0; cell<CellsPerTPC; cell++){
	//define variables to save highest bincontent and the bin index
	double HitsInFirstPeak=0;
	int FirstPeakBin=0;
	for (unsigned int bin=0; bin<40; ++bin){
	  if (HitsInFirstPeak<YPosXZ[cell]->GetBinContent(bin+1)){
	    HitsInFirstPeak=YPosXZ[cell]->GetBinContent(bin+1);
	    FirstPeakBin=bin;
	  }
	}
	
	//define the minimum and maximum range around the highest content bin
	double PeakBinXmin=double(FirstPeakBin)*((ybound[1]-ybound[0])/40.0)+ybound[0]-10;
	double PeakBinXmax=double(FirstPeakBin)*((ybound[1]-ybound[0])/40.0)+ybound[0]+10;

	//define a guassian histogram to fit the yposition histogram in the range calculated above
	TF1 *YPosFitFirstPeak = new TF1("YPosFitFirstPeak","gaus", PeakBinXmin, PeakBinXmax);
	YPosFitFirstPeak->SetParameter(0, 100);
	YPosFitFirstPeak->SetParameter(1, PeakBinXmin+10);
	YPosFitFirstPeak->SetParameter(2, 0.1);
	//YPosFitFirstPeak->SetLineColor(4);
	YPosXZ[cell]->Fit("YPosFitFirstPeak", "Q", "", PeakBinXmin-10, PeakBinXmax+10);
	//save interal, mean and rms in variables
	double FirstPeakBinAll=YPosFitFirstPeak->GetParameter(0);
	double FirstPeakBinMean=YPosFitFirstPeak->GetParameter(1);
	double FirstPeakBinRMS=YPosFitFirstPeak->GetParameter(2);
	
	FirstPeakXmin[cell]=FirstPeakBinMean-fabs(FirstPeakBinRMS)*2;
	FirstPeakXmax[cell]=FirstPeakBinMean+fabs(FirstPeakBinRMS)*2;
	
	//define the position range where most incorrectly disambiguated hits are residing.
	//depending on whether the mean position returned in above guassian fit is above or below
	//the midpoint of the detector volume in ydirection, the position range could be 100cm below
	//or above the mean position returned.
	if (FirstPeakBinMean>(ybound[0]+ybound[1])*0.5){
	  SecondPeakXmin[cell]=FirstPeakXmin[cell]-100;
	  SecondPeakXmax[cell]=FirstPeakXmax[cell]-100;
	} else if (FirstPeakBinMean<(ybound[0]+ybound[1])*0.5){
	  SecondPeakXmin[cell]=FirstPeakXmin[cell]+100;
	  SecondPeakXmax[cell]=FirstPeakXmax[cell]+100;
	}
	
	//define a guassian histogram to fit the yposition histogram in the range where most incorrectly
	//disambiguated hits are residing
	TF1 *YPosFitSecondPeak = new TF1("YPosFitSecondPeak","gaus", SecondPeakXmin[cell], SecondPeakXmax[cell]);
	YPosFitSecondPeak->SetParameter(0, 100);
	YPosFitSecondPeak->SetParameter(1, SecondPeakXmin[cell]+10);
	YPosFitSecondPeak->SetParameter(2, 0.1);
	//YPosFitSecondPeak->SetLineColor(4);
	YPosXZ[cell]->Fit("YPosFitSecondPeak", "Q", "", SecondPeakXmin[cell]-20, SecondPeakXmax[cell]+20);
	//define variables to save integral, mean and rms of the guassian fit result.
	double SecondPeakBinAll=YPosFitSecondPeak->GetParameter(0);
	double SecondPeakBinMean=YPosFitSecondPeak->GetParameter(1);
	double SecondPeakBinRMS=YPosFitSecondPeak->GetParameter(2);
	
	SecondPeakXmin[cell]=SecondPeakBinMean-fabs(SecondPeakBinRMS)*4;
	SecondPeakXmax[cell]=SecondPeakBinMean+fabs(SecondPeakBinRMS)*4;
	
	//we assumed that the position-correction function works when more than half of the induction plane
	//hits in the xz position cell is correctly disambiguated, meaning that most hits in the bin range whose
	//bincontents is highest is correctly disambiguated.  in this case, the hits in the highest bincontent
	//ranges are in the firstpeak range and the incorrectly disambiguated hits in the secondpeak range.  we
	//only correct the position of the hits in the secondpeak range if the hits in the firstpeak is more than 5. 
	ShouldCorrect[cell]=0;
	if (FirstPeakBinAll>SecondPeakBinAll &&
	    HitsInFirstPeak>5 && 
	    SecondPeakBinAll>0 && SecondPeakBinAll!=100 &&
	    SecondPeakXmin[cell]>ybound[0] && SecondPeakXmax[cell]<ybound[1]){
	  ShouldCorrect[cell]=1;
	}
	
	//delete histograms used
	YPosFitFirstPeak->Delete();
	YPosFitSecondPeak->Delete();
	YPosXZ[cell]->Delete();
		
      }//end of histogram loop


      //now we know the hits whose position and wireid to be corrected so loop over the induction plane hits
      //again, find the index of the hit in the xz position space and correct the position and wireid if the position
      //returned in the disambiguation algorithm above is in the range where most incorrectly disambiguated hits could be residing.
      for (unsigned int uv0=0; uv0 < InductionAll.size(); ++uv0){
	int XbinForHit=int( (InductionAll[uv0].PeakTime-xbound[0])/TotalXbin );
	int ZbinForHit=int( (InductionAll[uv0].FinalZPos-zbound[0])/TotalZbin );
	if (InductionAll[uv0].PeakTime>xbound[1]) XbinForHit=BinPerPos-1;
	if (InductionAll[uv0].FinalZPos<zbound[0]){
	  ZbinForHit=0;
	} else if (InductionAll[uv0].FinalZPos>zbound[1]){
	  ZbinForHit=BinPerPos-1;
	}
	int BinForHit=XbinForHit*BinPerPos+ZbinForHit;//bin index for the hit


	//shift the position and wireidreturned in the disambiguation algorithm
	std::vector<geo::WireID> wireiduv0 = geo->ChannelToWire(InductionAll[uv0].Channel);	
	if (ShouldCorrect[BinForHit]!=0){
	  if (InductionAll[uv0].FinalYPos>SecondPeakXmin[BinForHit] &&
	      InductionAll[uv0].FinalYPos<SecondPeakXmax[BinForHit]){
	    InductionAll[uv0].FinalYPos+=InductionAll[uv0].FinalYPos+
	    (FirstPeakXmin[BinForHit]+FirstPeakXmax[BinForHit])*0.5-(SecondPeakXmin[BinForHit]+SecondPeakXmax[BinForHit])*0.5;

	    if (InductionAll[uv0].DisambigWireID>1){
	      InductionAll[uv0].DisambigWireID=InductionAll[uv0].DisambigWireID-2;
	    } else if (InductionAll[uv0].DisambigWireID<2 && InductionAll[uv0].DisambigWireID+2<wireiduv0.size()){
	      InductionAll[uv0].DisambigWireID=InductionAll[uv0].DisambigWireID+2;
	    }
	  }
	}

	//save the hit information in a hit object
	recob::Hit hit(InductionAll[uv0].Channel,
		       InductionAll[uv0].StartTick,
		       InductionAll[uv0].EndTick,
		       InductionAll[uv0].PeakTime,
		       InductionAll[uv0].SigmaPeakTime,
		       InductionAll[uv0].RMS,
		       InductionAll[uv0].PeakAmplitude,
		       InductionAll[uv0].SigmaPeakAmplitude,
		       InductionAll[uv0].SummedADC,
		       InductionAll[uv0].Integral,
		       InductionAll[uv0].SigmaIntegral,
		       InductionAll[uv0].Multiplicity,
		       InductionAll[uv0].LocalIndex,
		       InductionAll[uv0].GoodnessOfFit,
		       InductionAll[uv0].DegreesOfFreedom,
		       InductionAll[uv0].View,
		       InductionAll[uv0].SignalType,
		       wireiduv0[InductionAll[uv0].DisambigWireID]);
	
	DisambiguatedHits.push_back(hit);

	//if (InductionAll[uv0].SimYPos!=99999){
	//HitsTotal+=InductionAll[uv0].Integral;
	//if (InductionAll[uv0].DisambigWireID==InductionAll[uv0].SimWireID) HitsCorrect+=InductionAll[uv0].Integral;
	//}
	
      }//end of induction plane hit loop

    }
  }

  //cout << "testing dune disambig" << "efficiency:" << double(HitsCorrect)*100.0/double(HitsTotal) << endl;

  
}

} //end namespace apa
