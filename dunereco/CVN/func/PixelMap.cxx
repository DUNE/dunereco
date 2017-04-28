////////////////////////////////////////////////////////////////////////
/// \file    PixelMap.h
/// \brief   PixelMap for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include "dune/CVN/func/PixelMap.h"

namespace cvn
{

  PixelMap::PixelMap(unsigned int nWire, unsigned int nTdc,
                     const Boundary& bound):
  fNWire(nWire),
  fNTdc(nTdc),
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

  void PixelMap::FillInputVector(float* input) const
  {
    unsigned int i = 0;

    for(const auto& pe:fPE){
      input[i] = pe;
      ++i;
    }

  }


  void PixelMap::Add(const unsigned int& wire, const double& tdc, const unsigned int& view, const double& pe)
  {
    const HitType label = kEmptyHit;
    const double purity=0.0;
    if(fBound.IsWithin(wire, tdc, view)){
      fPE[GlobalToIndex(wire,tdc, view)] += pe;
      fLab[GlobalToIndex(wire,tdc, view)] = label;
      fPur[GlobalToIndexSingle(wire,tdc, view)] = purity;
      if(view==0){
	fPEX[GlobalToIndexSingle(wire,tdc, view)] += pe;//Why +=?
	fLabX[GlobalToIndexSingle(wire,tdc, view)] = label;
	fPurX[GlobalToIndexSingle(wire,tdc, view)] = purity;
      }
      if(view==1){
	fPEY[GlobalToIndexSingle(wire,tdc, view)] += pe;
	fLabY[GlobalToIndexSingle(wire,tdc, view)] = label;
	fPurY[GlobalToIndexSingle(wire,tdc, view)] = purity;
      }
     if(view==2){
	fPEZ[GlobalToIndexSingle(wire,tdc, view)] += pe;
	fLabZ[GlobalToIndexSingle(wire,tdc, view)] = label;
	fPurZ[GlobalToIndexSingle(wire,tdc, view)] = purity;
      }
   }
  }

  unsigned int  PixelMap::GlobalToIndex(const unsigned int& wire,
                                        const double& tdc,
                                        const unsigned int& view)
  {

    unsigned int internalWire =  wire - fBound.FirstWire(view);

    double upperTL=fBound.LastTDC(view);
    double lowerTL=fBound.FirstTDC(view);
    double timestep=(upperTL-lowerTL)/double(fNTdc);
    double roundChannel=round((tdc-lowerTL)/timestep);

    unsigned int internalTdc  =  roundChannel;

    unsigned int index = internalWire * fNTdc + internalTdc % fNTdc;

    assert(index < fPE.size());

    return index;
  }

  unsigned int  PixelMap::LocalToIndex(const unsigned int& wire,
                                       const unsigned int& tdc) const
  {
    unsigned int index = wire * fNTdc + tdc % fNTdc;

    assert(index < fPE.size());
    return index;
  }

  unsigned int  PixelMap::GlobalToIndexSingle(const unsigned int& wire,
                                              const double& tdc,
                                              const unsigned int& view)

  {

    unsigned int internalWire =  wire - fBound.FirstWire(view);

    double upperTL=fBound.LastTDC(view);
    double lowerTL=fBound.FirstTDC(view);
    double timestep=(upperTL-lowerTL)/double(fNTdc);
    double roundChannel=round((tdc-lowerTL)/timestep);

    unsigned int internalTdc  =  roundChannel;

    unsigned int index = internalWire * fNTdc + internalTdc % fNTdc;

    assert(index < fPEX.size());

    return index;
  }

  void PixelMap::Print()
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

  TH2F* PixelMap::ToTH2() const
  {

    // Create a histogram, use twice as many tdcs to distinguish views
    TH2F* hist = new TH2F("PixelMap", ";Wire;Tdc", fNWire, 0, fNWire,
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

  TH2F* PixelMap::ToLabTH2() const
  {

    // Create a histogram, use twice as many tdcs to distinguish views
    TH2F* hist = new TH2F("PixelMap", ";Wire;Tdc", fNWire, 0, fNWire,
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

  TH2F* PixelMap::SingleViewToTH2(const unsigned int& view) const
  {

    // Create a histogram
    TH2F* hist = new TH2F("PixelMap", ";Wire;Tdc", fNWire, 0,
                          fNWire,
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

  std::ostream& operator<<(std::ostream& os, const PixelMap& m)
  {
    os << "PixelMap with " << m.NPixel() << " pixels, "
                           << m.NWire() << " wires"
                << " by "  << m.NTdc()  << " tdcs" ;
    return os;
  }
}
