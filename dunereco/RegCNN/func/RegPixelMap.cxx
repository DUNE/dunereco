////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap.cxx
/// \brief   RegPixelMap for RegCNN modifed from PixelMap.cxx
/// \author  Ilsoo Seong - iseong@uci.edu
///
/// \modified Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <iomanip>
#include "dune/RegCNN/func/RegPixelMap.h"

namespace cnn
{

  RegPixelMap::RegPixelMap(unsigned int nWire, unsigned int nWRes,
          unsigned int nTdc, unsigned int nTRes, const RegCNNBoundary& bound, const bool& prongOnly):
  fNWire(nWire),
  fNWRes(nWRes),
  fNTdc(nTdc),
  fNTRes(nTRes),
  fInPM(0),
  fTPC(0),
  fdist(100000),
  fPE(nWire*nTdc),
  fPEX(nWire*nTdc),
  fPEY(nWire*nTdc),
  fPEZ(nWire*nTdc),
  fPur(nWire*nTdc),
  fPurX(nWire*nTdc),
  fPurY(nWire*nTdc),
  fPurZ(nWire*nTdc),
  fLab(nWire*nTdc),
  fLabX(nWire*nTdc),
  fLabY(nWire*nTdc),
  fLabZ(nWire*nTdc),
  fBound(bound),
  fProngOnly(prongOnly),
  fProngTagX(nWire*nTdc),
  fProngTagY(nWire*nTdc),
  fProngTagZ(nWire*nTdc)
  {
      std::cout<<"here :"<<fNWire<<", "<<fNTdc<<", "<<fNWRes<<", "<<fNTRes<<std::endl;
  }

  void RegPixelMap::FillInputVector(float* input) const
  {
    unsigned int i = 0;

    for(const auto& pe:fPE){
      input[i] = pe;
      ++i;
    }

  }


  void RegPixelMap::Add(const int& wire, const int& tdc, const unsigned int& view, const double& pe, const unsigned int& tpc, int hit_prong_tag)
  {
    // keep these for now although we only use fPE
    const HitType label = kEmptyHit;
    const double purity=0.0;
    if(fBound.IsWithin(wire, tdc, view)){
      fInPM = 1; // any hit within the boundary
      GetTPC(wire, tdc, view, tpc);
      fPE[GlobalToIndex(wire, tdc, view)] += (float)pe;
      fLab[GlobalToIndex(wire,tdc, view)] = label;
      fPur[GlobalToIndex(wire,tdc, view)] = purity;
      if(view==0){
	    fPEX[GlobalToIndexSingle(wire,tdc, view)] += (float)pe;
        fProngTagX[GlobalToIndexSingle(wire, tdc, view)] = hit_prong_tag;
	    fLabX[GlobalToIndexSingle(wire,tdc, view)] = label;
	    fPurX[GlobalToIndexSingle(wire,tdc, view)] = purity;
        // FIXIT
        //if (pe!=0) {
        //    std::cout<<view<<" | ";
        //    std::cout<<wire<<", "<<tdc<<", "<<GlobalToIndexSingle(wire,tdc,view)<<", "<<pe<<std::endl;
        //}
      }
      if(view==1){
	    fPEY[GlobalToIndexSingle(wire,tdc, view)] += (float)pe;
        fProngTagY[GlobalToIndexSingle(wire, tdc, view)] = hit_prong_tag;
	    fLabY[GlobalToIndexSingle(wire,tdc, view)] = label;
	    fPurY[GlobalToIndexSingle(wire,tdc, view)] = purity;
      }
     if(view==2){
	    fPEZ[GlobalToIndexSingle(wire,tdc, view)] += (float)pe;
        fProngTagZ[GlobalToIndexSingle(wire, tdc, view)] = hit_prong_tag;
	    fLabZ[GlobalToIndexSingle(wire,tdc, view)] = label;
	    fPurZ[GlobalToIndexSingle(wire,tdc, view)] = purity;
        //std::cout << "Q = " << pe << " " << fPEZ[GlobalToIndexSingle(wire,tdc, view)] << std::endl;
      }
   }
  }

  void RegPixelMap::Finish() {
      // fProngOnly=True means only the primary prong is selected to creat pixel maps 
      //            and that prong needs to be either a muon or antimuon (FIXIT)?
      // Caveat 1: the prong tag of that pixel is determined by the track id of the track 
      //         associlated with the last hit that associated with that pixel. It means
      //         if this pixel-associated hits belong to different prongs, we could 
      //         either add deposited energy from other prongs (if other prongs' hit is
      //         in the front), or throw the energy from primary prong away (if there is
      //         one or more hits from other prongs in the last). A better way could be 
      //         determining the prong tag by looping the hits, instead of the pixels. 
      //         That means we should create a new pixel map only with the spacepoints 
      //         associated with the primary prong, instead of using the prong tag
      //         (little effect on the results, ignored for now)
      if (fProngOnly) {
          std::cout<<"Do Prong Only selection ......"<<std::endl;
          for (unsigned int i_p= 0; i_p< fPE.size(); ++i_p) {
              if (fProngTagX[i_p] != 0)
                  fPEX[i_p] = 0;
              if (fProngTagY[i_p] != 0)
                  fPEY[i_p] = 0;
              if (fProngTagZ[i_p] != 0)
                  fPEZ[i_p] = 0;
          } // end of i_p
      } // end of fProngOnly
  }

  unsigned int  RegPixelMap::GlobalToIndex(const int& wire,
                                        const int& tdc,
                                        const unsigned int& view)
  {

    //int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view))/2;
    int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view)+(int)fNWRes)/2;
    int meanTDC  = (fBound.LastTDC(view)+fBound.FirstTDC(view)+(int)fNTRes/2)/2;

    //unsigned int internalWire = (unsigned int)( (wire-meanWire) + fNWire/2 );
    unsigned int internalWire = (unsigned int)(round((float)((unsigned int)(wire-meanWire)+fNWire*fNWRes/2)/(float)fNWRes) );
    //unsigned int internalTdc  = (unsigned int)( round((float)(tdc-meanTDC)/(float)fNTRes + fNTdc/2) );
    unsigned int internalTdc  = (unsigned int)( round((float)((unsigned int)(tdc-meanTDC)+fNTdc*fNTRes/2)/(float)fNTRes) );

    unsigned int index = internalWire * fNTdc + internalTdc % fNTdc;
   
    //if (internalTdc <1.1){
    //std::cout << "====> " << meanWire <<" " << meanTDC << " " << internalWire << " " << internalTdc << " " << index << std::endl;
    //std::cout << "    => " << wire << " " << tdc << " " << wire-meanWire << " " << round((float)(tdc-meanTDC)/(float)fNTRes) << std::endl;
    //}
    assert(index < fPE.size());
    return index;
  }

  unsigned int  RegPixelMap::LocalToIndex(const unsigned int& wire,
                                       const unsigned int& tdc) const
  {
    unsigned int index = wire * fNTdc + tdc % fNTdc;

    assert(index < fPE.size());
    return index;
  }

  unsigned int  RegPixelMap::GlobalToIndexSingle(const int& wire,
                                              const int& tdc,
                                              const unsigned int& view)

  {

    //int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view))/2;
    int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view)+(int)fNWRes/2)/2;
    int meanTDC  = (fBound.LastTDC(view)+fBound.FirstTDC(view)+(int)fNTRes/2)/2;

    //unsigned int internalWire = (unsigned int)( (wire-meanWire) + fNWire/2 );
    //unsigned int internalWire = (unsigned int)( round((float)((unsigned int)(wire-meanWire)+fNWire*fNWRes/2)/(float)fNWRes) );
    unsigned int internalWire = (unsigned int)(round((float)(wire - meanWire + fNWire*fNWRes/2)/(float)fNWRes));

    //unsigned int internalTdc  = (unsigned int)( round( (float)(tdc-meanTDC)/(float)fNTRes + fNTdc/2) );
    //unsigned int internalTdc  = (unsigned int)( round((float)((unsigned int)(tdc-meanTDC)+fNTdc*fNTRes/2)/(float)fNTRes) );
    unsigned int internalTdc  = (unsigned int)(round((float)(tdc - meanTDC + fNTdc*fNTRes/2)/(float)fNTRes));

    unsigned int index = internalWire * fNTdc + internalTdc % fNTdc;

    assert(index < fPEX.size());

    return index;
  }

  void  RegPixelMap::GetTPC(const int& wire,
                             const int& tdc,
                             const unsigned int& view,
                             const unsigned int& tpc)

  {

    //int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view))/2;
    int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view)+(int)fNWRes)/2;
    int meanTDC  = (fBound.LastTDC(view)+fBound.FirstTDC(view)+(int)fNTRes/2)/2;
    double dist = pow((wire-meanWire)*0.5,2)+pow((tdc-meanTDC)*0.08,2);
    dist = sqrt(dist);
    if (dist < fdist){
        fdist = dist;
        fTPC = tpc;
    }
  }


  void RegPixelMap::Print()
  {

    // Start by doing even wires
    for(unsigned int iTdc = 0; iTdc < fNTdc; ++iTdc)
    {
      for(unsigned int iWire = 0; iWire < fNWire; iWire += 2)
      {
        unsigned int index = LocalToIndex(iWire, iTdc);
        if( fPE[index] > 0)
        {
          std::cout << "*";
        }
        else
        {
          std::cout << " ";
        }

      }
      std::cout << std::endl;
    }
    // Then do odd wires
    for(unsigned int iTdc = 0; iTdc < fNTdc; ++iTdc)
    {
      for(unsigned int iWire = 1; iWire < fNWire; iWire += 2)
      {
        unsigned int index = LocalToIndex(iWire, iTdc);
        if( fPE[index] > 0)
        {
          std::cout << "*";
        }
        else
        {
          std::cout << " ";
        }

      }
      std::cout << std::endl;
    }

  }

  TH2F* RegPixelMap::ToTH2() const
  {

    // Create a histogram, use twice as many tdcs to distinguish views
    TH2F* hist = new TH2F("RegPixelMap", ";Wire;Tdc", fNWire, 0, fNWire,
                                                      fNTdc*3, 0, fNTdc*3);

    for(unsigned int iWire = 0; iWire < fNWire; ++iWire)
    {
      for(unsigned int iTdc = 0; iTdc < fNTdc; ++iTdc)
      {
        // Add 1 to in each bin to skip underflow
        hist->SetBinContent(iWire+1, iTdc + fNTdc*(iWire%3) + 1,
                            fPE[LocalToIndex(iWire, iTdc)]);

      }
    }
    return hist;
  }

  TH2F* RegPixelMap::ToLabTH2() const
  {

    // Create a histogram, use twice as many tdcs to distinguish views
    TH2F* hist = new TH2F("RegPixelMap", ";Wire;Tdc", fNWire, 0, fNWire,
                                                      fNTdc*3, 0, fNTdc*3);

    for(unsigned int iWire = 0; iWire < fNWire; ++iWire)
    {
      for(unsigned int iTdc = 0; iTdc < fNTdc; ++iTdc)
      {
        // Add 1 to in each bin to skip underflow
        hist->SetBinContent(iWire+1, iTdc + fNTdc*(iWire%3) + 1,
                            (double)fLab[LocalToIndex(iWire, iTdc)]);

      }
    }
    return hist;
  }

  TH2F* RegPixelMap::SingleViewToTH2(const unsigned int& view) const
  {

    // Create a histogram
    TH2F* hist = new TH2F("RegPixelMap", ";Wire;Tdc", fNWire, 0, fNWire,
                          			      fNTdc, 0, fNTdc);

    for(unsigned int iWire = 0; iWire < fNWire; ++iWire)
    {
      for(unsigned int iTdc = 0; iTdc < fNTdc; ++iTdc)
      {
        // Add 1 to in each bin to skip underflow
	if(view==0){
        hist->SetBinContent(iWire+1, iTdc + 1,
                            fPEX[LocalToIndex(iWire, iTdc)]);
        // FIXIT
        //if(fPEX[LocalToIndex(iWire, iTdc)]!=0) {
        //    std::cout<<iWire<<", "<<iTdc<<", "<<fPEX[LocalToIndex(iWire, iTdc)]/500<<std::endl;
        //}
	}
        if(view==1){
        hist->SetBinContent(iWire+1, iTdc + 1,
                            fPEY[LocalToIndex(iWire, iTdc)]);
	}
        if(view==2){
          hist->SetBinContent(iWire+1, iTdc + 1,
                              fPEZ[LocalToIndex(iWire, iTdc)]);
	}
      }
    }
    return hist;
  }

  std::ostream& operator<<(std::ostream& os, const RegPixelMap& m)
  {
    os << "RegPixelMap with " << m.NPixel() << " pixels, "
                           << m.NWire() << " wires"
                << " by "  << m.NTdc()  << " tdcs" ;
    return os;
  }
}
