////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap3D.cxx
/// \brief   RegPixelMap3D for RegCNN modifed from PixelMap.cxx
/// \author  Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <iomanip>
#include "dune/RegCNN/func/RegPixelMap3D.h"

namespace cnn
{

  RegPixelMap3D::RegPixelMap3D(const RegCNNBoundary3D& bound, const bool& cropped, const bool& prongOnly):
      fBound(bound),
      fCropped(cropped),
      fProngOnly(prongOnly),
      fInPM(0),
      fPE(fBound.NBins(0)*fBound.NBins(1)*fBound.NBins(2)),
      fPECropped(32*32*32),                   // Fixed to 32x32x32
      fProngTag(fBound.NBins(0)*fBound.NBins(1)*fBound.NBins(2))
  {
      x_axis.Set(fBound.NBins(0), fBound.StartPos(0), fBound.StopPos(0));
      y_axis.Set(fBound.NBins(1), fBound.StartPos(1), fBound.StopPos(1));
      z_axis.Set(fBound.NBins(2), fBound.StartPos(2), fBound.StopPos(2));
  }

  void RegPixelMap3D::AddHit(float rel_x, float rel_y, float rel_z, float charge, int hit_prong_tag) 
  {
      if (fBound.IsWithin(rel_x, rel_y, rel_z)) {
          fInPM = 1;
          int xbin = x_axis.FindBin(rel_x);
          int ybin = y_axis.FindBin(rel_y);
          int zbin = z_axis.FindBin(rel_z);
          fPE[LocalToIndex(xbin-1, ybin-1, zbin-1)] += charge;
          fProngTag[LocalToIndex(xbin-1, ybin-1, zbin-1)] = hit_prong_tag;
      }
  }

  void RegPixelMap3D::Finish() {
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
              if (fProngTag[i_p] != 0)
                  fPE[i_p] = 0;
          } // end of i_p
      } // end of fProngOnly

      if (fCropped) {
          std::cout<<"Crop pixel size to 32x32x32 ......"<<std::endl;
          int cropped_xbin_low = 50 - 16;    //int cropped_xbin_high = 50 + 16;
          int cropped_ybin_low = 50 - 16;    //int cropped_ybin_high = 50 + 16;
          int cropped_zbin_low = 0;          //int cropped_zbin_high = 32;
          for (int i_x= 0; i_x< 32; ++i_x) {
              for (int i_y= 0; i_y< 32; ++i_y) {
                  for (int i_z= 0; i_z< 32; ++i_z) {
                      int cropped_index = i_x*32*32 + i_y*32 + i_z%32;
                      int ii_x  = i_x + cropped_xbin_low;
                      int ii_y  = i_y + cropped_ybin_low;
                      int ii_z  = i_z + cropped_zbin_low;
                      int index =  LocalToIndex(ii_x, ii_y, ii_z);
                      fPECropped[cropped_index] = fPE[index];
                  } // end of i_z
              } // end of i_y
          } // end of i_x
      } // end of fCropped
  }

  unsigned int RegPixelMap3D::LocalToIndex(const unsigned int& bin_x, 
          const unsigned int& bin_y, const unsigned int& bin_z) const 
  {
      unsigned int index = bin_x*fBound.NBins(1)*fBound.NBins(2) + bin_y*fBound.NBins(2) + bin_z%fBound.NBins(2);
      assert(index < fPE.size());
      return index;
  }

  TH3F* RegPixelMap3D::ToTH3() const {
      // Create a histogram
      TH3F* hist = new TH3F("RegPixelMap3D", "X:Y:Z", fBound.NBins(0), fBound.StartPos(0), fBound.StopPos(0),
              fBound.NBins(1), fBound.StartPos(1), fBound.StopPos(1),
              fBound.NBins(2), fBound.StartPos(2), fBound.StopPos(2));
      for (int ix= 0; ix< fBound.NBins(0); ++ix) {
          for (int iy= 0; iy< fBound.NBins(1); ++iy) {
              for (int iz= 0; iz< fBound.NBins(2); ++iz) {
                  hist->SetBinContent(ix+1, iy+1, iz+1, fPE[LocalToIndex(ix, iy, iz)]);
              }
          }
      }
      return hist;
  }

  TH3F* RegPixelMap3D::ToCroppedTH3() const {
      // Create a histogram
      TH3F* hist = new TH3F("RegCroppedPixelMap3D", "X:Y:Z", 32, 0, 32*x_axis.GetBinWidth(0),
              32, 0, 32*y_axis.GetBinWidth(0),
              32, 0, 32*z_axis.GetBinWidth(0));
      for (int ix= 0; ix< 32; ++ix) {
          for (int iy= 0; iy< 32; ++iy) {
              for (int iz= 0; iz< 32; ++iz) {
                  int cropped_index = ix*32*32 + iy*32 + iz%32;
                  hist->SetBinContent(ix+1, iy+1, iz+1, fPECropped[cropped_index]);
              }
          }
      }
      return hist;
  }

  std::ostream& operator<<(std::ostream& os, const RegPixelMap3D& m)
  {
    os << "RegPixelMap3D...... ";
    return os;
  }
}
