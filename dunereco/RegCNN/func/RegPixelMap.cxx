////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap.cxx
/// \brief   RegPixelMap for RegCNN modifed from PixelMap.cxx
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <iomanip>
#include "dune/RegCNN/func/RegPixelMap.h"

namespace cnn
{

  RegPixelMap::RegPixelMap(unsigned int nWire, unsigned int nTdc, unsigned int nTRes,
                     const RegCNNBoundary& bound):
  fNWire(nWire),
  fNTdc(nTdc),
  fNTRes(nTRes),
  fInPM(0),
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
  fBound(bound)
  {}

  void RegPixelMap::FillInputVector(float* input) const
  {
    unsigned int i = 0;

    for(const auto& pe:fPE){
      input[i] = pe;
      ++i;
    }

  }


  void RegPixelMap::Add(const int& wire, const int& tdc, const unsigned int& view, const double& pe)
  {
    // keep these for now although we only use fPE
    const HitType label = kEmptyHit;
    const double purity=0.0;
    if(fBound.IsWithin(wire, tdc, view)){
      fInPM = 1; // any hit within the boundary
      fPE[GlobalToIndex(wire,tdc, view)] += (float)pe;
      fLab[GlobalToIndex(wire,tdc, view)] = label;
      fPur[GlobalToIndexSingle(wire,tdc, view)] = purity;
      if(view==0){
	fPEX[GlobalToIndexSingle(wire,tdc, view)] += (float)pe;
	fLabX[GlobalToIndexSingle(wire,tdc, view)] = label;
	fPurX[GlobalToIndexSingle(wire,tdc, view)] = purity;
      }
      if(view==1){
	fPEY[GlobalToIndexSingle(wire,tdc, view)] += (float)pe;
	fLabY[GlobalToIndexSingle(wire,tdc, view)] = label;
	fPurY[GlobalToIndexSingle(wire,tdc, view)] = purity;

      }
     if(view==2){
	fPEZ[GlobalToIndexSingle(wire,tdc, view)] += (float)pe;
	fLabZ[GlobalToIndexSingle(wire,tdc, view)] = label;
	fPurZ[GlobalToIndexSingle(wire,tdc, view)] = purity;
      //std::cout << "Q = " << pe << " " << fPEZ[GlobalToIndexSingle(wire,tdc, view)] << std::endl;
      }
   }
  }

  unsigned int  RegPixelMap::GlobalToIndex(const int& wire,
                                        const int& tdc,
                                        const unsigned int& view)
  {

    int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view))/2;
    int meanTDC  = (fBound.LastTDC(view)+fBound.FirstTDC(view)+(int)fNTRes/2)/2;

    unsigned int internalWire = (unsigned int)( (wire-meanWire) + fNWire/2 );
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

    int meanWire = (fBound.LastWire(view)+fBound.FirstWire(view))/2;
    int meanTDC  = (fBound.LastTDC(view)+fBound.FirstTDC(view)+(int)fNTRes/2)/2;

    unsigned int internalWire = (unsigned int)( (wire-meanWire) + fNWire/2 );
    //unsigned int internalTdc  = (unsigned int)( round( (float)(tdc-meanTDC)/(float)fNTRes + fNTdc/2) );
    unsigned int internalTdc  = (unsigned int)( round((float)((unsigned int)(tdc-meanTDC)+fNTdc*fNTRes/2)/(float)fNTRes) );

    unsigned int index = internalWire * fNTdc + internalTdc % fNTdc;

    assert(index < fPEX.size());

    return index;
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
