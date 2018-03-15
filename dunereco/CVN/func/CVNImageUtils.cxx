#include <vector>
#include <iostream>

#include "dune/CVN/func/CVNImageUtils.h"

cvn::CVNImageUtils::CVNImageUtils(){
  // Set a default image size
  SetImageSize(500,500,3);
  SetPixelMapSize(2880,500);
  // Defualt is to reverse the y view
  fViewReverse = {false,true,false};

  fUseLogScale = false;
}

cvn::CVNImageUtils::CVNImageUtils(unsigned int nWires, unsigned int nTDCs, unsigned int nViews){
  SetImageSize(nWires,nTDCs,nViews);
  SetPixelMapSize(2880,500);
  fUseLogScale = false;
}

cvn::CVNImageUtils::~CVNImageUtils(){

}

unsigned char cvn::CVNImageUtils::ConvertChargeToChar(float charge){

  float peCorrChunk;
  float truncateCorr;
  float centreScale = 0.7;
  if(fUseLogScale){
    float scaleFrac=(log(charge)/log(1000));
    truncateCorr= ceil(centreScale*scaleFrac*255.0);
      }
  else{
    peCorrChunk = (1000.) / 255.0;
    truncateCorr = ceil((charge)/(peCorrChunk));
  }
  if (truncateCorr > 255) return (unsigned char)255;
  else return (unsigned char)truncateCorr;

}

void cvn::CVNImageUtils::SetImageSize(unsigned int nWires, unsigned int nTDCs, unsigned int nViews){
  fNWires = nWires;
  fNTDCs = nTDCs;
  fNViews = nViews;
}

void cvn::CVNImageUtils::SetViewReversal(bool reverseX, bool reverseY, bool reverseZ){
  fViewReverse = {reverseX,reverseY,reverseZ};
}

void cvn::CVNImageUtils::SetLogScale(bool setLog){
  fUseLogScale = setLog;
}

void cvn::CVNImageUtils::SetPixelMapSize(unsigned int nWires, unsigned int nTDCs){
  fPixelMapWires = nWires;
  fPixelMapTDCs = nTDCs;
}

void cvn::CVNImageUtils::ConvertPixelMapToPixelArray(PixelMap &pm, std::vector<unsigned char> &pix){

  SetPixelMapSize(pm.fNWire,pm.fNTdc);
  // Strip out the charge vectors and use these
  ConvertChargeVectorsToPixelArray(pm.fPEX,pm.fPEY,pm.fPEZ,pix);

}


void cvn::CVNImageUtils::ConvertChargeVectorsToPixelArray(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                                          std::vector<float> &v2pe, std::vector<unsigned char> &pix){

  // Reverse requested views
  if(fViewReverse[0]) ReverseView(v0pe);
  if(fViewReverse[1]) ReverseView(v1pe);
  if(fViewReverse[2]) ReverseView(v2pe);
 
  // Calculate the integrated charge in each plane.
  std::vector< std::vector<float> > chargeVec;

  // Loop over v views, w wires and t tdcs 
  for (unsigned int view = 0; view < fNViews; ++view)
  {
    std::vector<float> tempChargeVec;

    for (unsigned int wire = 0; wire < fPixelMapWires; ++wire)
    {
      float totCharge = 0;

      for (unsigned int time = 0; time < fPixelMapTDCs; ++time)
      {
        float val = 0.;
        unsigned int element = time + fPixelMapTDCs * wire;
        if(view == 0 ){
          val = v0pe[element];
        }
        if(view == 1 ){
          val = v1pe[element];
        }
        if(view == 2 ){
          val = v2pe[element];
        }
        totCharge += val;
      }
      tempChargeVec.push_back(totCharge);
    }
    chargeVec.push_back(tempChargeVec);
  }

  // Now we want start and end planes for each view
  std::vector<unsigned int> imageStartWire(3,0);
  std::vector<unsigned int> imageEndWire(3,0);

  for(unsigned int view = 0; view < chargeVec.size(); ++view){
    for(unsigned int wire = 0; wire < chargeVec[view].size(); ++wire){

      // If we have got to fNWires from the end, the start needs to be this wire
      if(chargeVec[view].size() - wire == fNWires){
        imageStartWire[view] = wire;
        imageEndWire[view] = wire + fNWires - 1;
        break;
      }

      // For a given plane, look to see if the next 20 planes are empty. If not, this can be out start plane.
      int nEmpty = 0;
      for(unsigned int nextWire = wire + 1; nextWire <= wire + 20; ++nextWire){
        if(chargeVec[view][nextWire] == 0.0) ++nEmpty;
      }
      if(nEmpty < 5){
        imageStartWire[view] = wire;
        imageEndWire[view] = wire + fNWires - 1;
        break;
      }
    }
  }

  // Actually write the values to the pixel array
  for (unsigned int view = 0; view < fNViews; ++view)
  {
    for (unsigned int wire = imageStartWire[view]; wire < imageEndWire[view]; ++wire)
    {
      for (unsigned int time = 0; time < fNTDCs; ++time)
      {
        unsigned int i = time + fNTDCs * ((wire - imageStartWire[view]) + fNWires * view);
        float val =0.;
        unsigned int element = time + fNTDCs * wire;
        if(view == 0 ){
          val = v0pe[element];
        }
        if(view == 1 ){
          val = v1pe[element];
        }
        if(view == 2 ){
          val = v2pe[element];
        }
        pix[i] = ConvertChargeToChar(val);
      }
    }
  }

  return;

}

void cvn::CVNImageUtils::ReverseView(std::vector<float> &peVec){

  std::vector<float> vecCopy(peVec.size(),0.); 

  for (unsigned int w = 0; w < fPixelMapWires; ++w)
  {
    // Get our new plane number
    unsigned int newPlane = fPixelMapWires - w - 1;

    for (unsigned int t = 0; t < fPixelMapTDCs; ++t)
    {
      float val = peVec[t + fPixelMapTDCs * w];
      vecCopy[t + fPixelMapTDCs * newPlane] = val;
    }
  }

  // Copy the values back into the original vector
  for(unsigned int e = 0; e < peVec.size(); ++e){
    float val = vecCopy[e];
    peVec[e] = val;
  }

}

