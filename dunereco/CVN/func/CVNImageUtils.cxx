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

void cvn::CVNImageUtils::SetViewReversal(std::vector<bool> reverseViews){
  if(reverseViews.size() != 3){
    std::cout << "Expected three views for view reversals... using defaults." << std::endl;
  }
  else{
    SetViewReversal(reverseViews[0],reverseViews[1],reverseViews[2]);
  }
  return;
}

void cvn::CVNImageUtils::SetLogScale(bool setLog){
  fUseLogScale = setLog;
}

void cvn::CVNImageUtils::SetPixelMapSize(unsigned int nWires, unsigned int nTDCs){
  fPixelMapWires = nWires;
  fPixelMapTDCs = nTDCs;
}

void cvn::CVNImageUtils::ConvertPixelMapToPixelArray(const PixelMap &pm, std::vector<unsigned char> &pix){

  SetPixelMapSize(pm.fNWire,pm.fNTdc);

  // Strip out the charge vectors and use these
  std::vector<float> v0pe = pm.fPEX;
  std::vector<float> v1pe = pm.fPEY;
  std::vector<float> v2pe = pm.fPEZ;

  ConvertChargeVectorsToPixelArray(v0pe,v1pe,v2pe,pix);

}


void cvn::CVNImageUtils::ConvertChargeVectorsToPixelArray(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                                          std::vector<float> &v2pe, std::vector<unsigned char> &pix){

  // Get the vectors
  cvn::ViewVector view0;
  cvn::ViewVector view1;
  cvn::ViewVector view2;
  ConvertChargeVectorsToViewVectors(v0pe, v1pe, v2pe, view0, view1, view2); 

  // Actually write the values to the pixel array
  for (unsigned int view = 0; view < fNViews; ++view){
    for (unsigned int wire = 0; wire < fNWires; ++wire){
      for (unsigned int time = 0; time < fNTDCs; ++time){

        unsigned char val = 0;
        // Get the index for the pixel map
        if(view == 0 ){ val = view0[wire][time]; }
        if(view == 1 ){ val = view1[wire][time]; }
        if(view == 2 ){ val = view2[wire][time]; }

        // Get the index for the final image
        unsigned int i = time + fNTDCs * (wire + fNWires * view);
        pix[i] = val;
  
      }
    }
  }

  return;

}

void cvn::CVNImageUtils::ConvertPixelMapToImageVector(const cvn::PixelMap &pm, cvn::ImageVector &imageVec){

  SetPixelMapSize(pm.fNWire,pm.fNTdc);

  // Strip out the charge vectors and use these
  std::vector<float> v0pe = pm.fPEX;
  std::vector<float> v1pe = pm.fPEY;
  std::vector<float> v2pe = pm.fPEZ;

  ConvertChargeVectorsToImageVector(v0pe, v1pe, v2pe, imageVec);
}

void cvn::CVNImageUtils::ConvertPixelMapToImageVectorF(const cvn::PixelMap &pm, cvn::ImageVectorF &imageVec){

  SetPixelMapSize(pm.fNWire,pm.fNTdc);

  // Strip out the charge vectors and use these
  std::vector<float> v0pe = pm.fPEX;
  std::vector<float> v1pe = pm.fPEY;
  std::vector<float> v2pe = pm.fPEZ;

  ConvertChargeVectorsToImageVectorF(v0pe, v1pe, v2pe, imageVec);
}

void cvn::CVNImageUtils::ConvertChargeVectorsToImageVector(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                                           std::vector<float> &v2pe, cvn::ImageVector &imageVec){

  cvn::ViewVector view0;
  cvn::ViewVector view1;
  cvn::ViewVector view2;

  ConvertChargeVectorsToViewVectors(v0pe, v1pe, v2pe, view0, view1, view2);

  cvn::ImageVector newImage = BuildImageVector(view0,view1,view2);

  imageVec = newImage;
}

void cvn::CVNImageUtils::ConvertChargeVectorsToImageVectorF(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                                           std::vector<float> &v2pe, cvn::ImageVectorF &imageVec){

  cvn::ViewVector view0;
  cvn::ViewVector view1;
  cvn::ViewVector view2;

  ConvertChargeVectorsToViewVectors(v0pe, v1pe, v2pe, view0, view1, view2);

  // Convert the ViewVector to ViewVectorF
  cvn::ViewVectorF floatView0 = ConvertViewVecToViewVecF(view0);
  cvn::ViewVectorF floatView1 = ConvertViewVecToViewVecF(view1);
  cvn::ViewVectorF floatView2 = ConvertViewVecToViewVecF(view2);

  cvn::ImageVectorF newImage = BuildImageVectorF(floatView0,floatView1,floatView2);

  imageVec = newImage;
}


void cvn::CVNImageUtils::ConvertChargeVectorsToViewVectors(std::vector<float> &v0pe, std::vector<float> &v1pe, std::vector<float> &v2pe,
                                       cvn::ViewVector& view0, cvn::ViewVector& view1, cvn::ViewVector& view2){

  // Reverse requested views
  if(fViewReverse[0]) ReverseView(v0pe);
  if(fViewReverse[1]) ReverseView(v1pe);
  if(fViewReverse[2]) ReverseView(v2pe);

  // Get the integrated charge for each wire
  std::vector< std::vector<float> > wireCharges;
  for (unsigned int view = 0; view < fNViews; ++view){

    std::vector<float> tempChargeVec;
    for (unsigned int wire = 0; wire < fPixelMapWires; ++wire){

      float totCharge = 0;
      for (unsigned int time = 0; time < fPixelMapTDCs; ++time){
        float val = 0.;
        unsigned int element = time + fPixelMapTDCs * wire;
        if(view == 0 ){ val = v0pe[element]; }
        if(view == 1 ){ val = v1pe[element]; }
        if(view == 2 ){ val = v2pe[element]; }
        totCharge += val;
      }
      tempChargeVec.push_back(totCharge);
    }
    wireCharges.push_back(tempChargeVec);
  }

  // Get the integrated charge for each tdc
  std::vector< std::vector<float> > tdcCharges;
  for (unsigned int view = 0; view < fNViews; ++view){

    std::vector<float> tempChargeVec;
    for (unsigned int time = 0; time < fPixelMapTDCs; ++time){


      float totCharge = 0;
      for (unsigned int wire = 0; wire < fPixelMapWires; ++wire){

        float val = 0.;
        unsigned int element = time + fPixelMapTDCs * wire;
        if(view == 0 ){ val = v0pe[element]; }
        if(view == 1 ){ val = v1pe[element]; }
        if(view == 2 ){ val = v2pe[element]; }
        totCharge += val;
      }
      tempChargeVec.push_back(totCharge);
    }
    tdcCharges.push_back(tempChargeVec);
  }

  // The output image consists of a rectangular region of the pixel map
  // We want to find the start and end wires for each view
  std::vector<unsigned int> imageStartWire(3,0);
  std::vector<unsigned int> imageEndWire(3,0);
  // And the same for TDCs
  std::vector<unsigned int> imageStartTDC(3,0);
  std::vector<unsigned int> imageEndTDC(3,0);

  for(unsigned int view = 0; view < wireCharges.size(); ++view){
    GetMinMaxWires(wireCharges[view],imageStartWire[view],imageEndWire[view]);
    GetMinMaxTDCs(tdcCharges[view],imageStartTDC[view],imageEndTDC[view]);

//    std::cout << " Wires: " << imageStartWire[view] << ", " << imageEndWire[view] << " :: TDCs: "
//                            << imageStartTDC[view]  << ", " << imageEndTDC[view]  << std::endl;
  }

  // Write the values to the three vectors
  for (unsigned int view = 0; view < fNViews; ++view){
    cvn::ViewVector viewChargeVec;
    for (unsigned int wire = imageStartWire[view]; wire <= imageEndWire[view]; ++wire){
      std::vector<unsigned char> wireTDCVec;
      for (unsigned int time = imageStartTDC[view]; time <= imageEndTDC[view]; ++time){

        // Get the index for the pixel map
        unsigned int element = time + fPixelMapTDCs * wire;

        // We have to convert to char and then convert back to a float
        unsigned char val = 0;
        if(view == 0){ val = ConvertChargeToChar(v0pe[element]); }
        if(view == 1){ val = ConvertChargeToChar(v1pe[element]); }
        if(view == 2){ val = ConvertChargeToChar(v2pe[element]); }
        wireTDCVec.push_back(val);
      }
      viewChargeVec.push_back(wireTDCVec);
    }
    if(view == 0) view0 = viewChargeVec;
    if(view == 1) view1 = viewChargeVec;
    if(view == 2) view2 = viewChargeVec;
  }

  return;

}

void cvn::CVNImageUtils::ConvertPixelArrayToImageVectorF(const std::vector<unsigned char> &pixelArray, cvn::ImageVectorF &imageVec){

  // The pixel arrays is built with indices i = tdc + nTDCs(wire + nWires*view)

  cvn::ViewVectorF view0; 
  cvn::ViewVectorF view1; 
  cvn::ViewVectorF view2;
 
  for(unsigned int v = 0; v < fNViews; ++v){
    for(unsigned int w = 0; w < fNWires; ++w){
      std::vector<float> wireVec;
      for(unsigned int t = 0; t < fNTDCs; ++t){
        unsigned int index = t + fNTDCs*(w + fNWires*v);
        wireVec.push_back(pixelArray[index]);
      }
      if(v==0) view0.push_back(wireVec);
      if(v==1) view1.push_back(wireVec);
      if(v==2) view2.push_back(wireVec);
    }
  }

  imageVec = BuildImageVectorF(view0,view1,view2);

}

void cvn::CVNImageUtils::GetMinMaxWires(std::vector<float> &wireCharges, unsigned int &minWire, unsigned int &maxWire){

  minWire = 0;
  maxWire = fNWires;

  for(unsigned int wire = 0; wire < wireCharges.size(); ++wire){

    // If we have got to fNWires from the end, the start needs to be this wire
    if(wireCharges.size() - wire == fNWires){
      minWire = wire;
      maxWire = wire + fNWires - 1;
      return;
    }

    // For a given plane, look to see if the next 20 planes are empty. If not, this can be out start plane.
    int nEmpty = 0;
    for(unsigned int nextWire = wire + 1; nextWire <= wire + 20; ++nextWire){
      if(wireCharges[nextWire] == 0.0) ++nEmpty;
    }
    if(nEmpty < 5){
      minWire = wire;
      maxWire = wire + fNWires - 1;
      return;
    }
  }

}

void cvn::CVNImageUtils::GetMinMaxTDCs(std::vector<float> &tdcCharges, unsigned int &minTDC, unsigned int &maxTDC){

  minTDC = 0;
  maxTDC = fNTDCs;

  for(unsigned int tdc = 0; tdc < tdcCharges.size(); ++tdc){

    // If we have got to fNWires from the end, the start needs to be this wire
    if(tdcCharges.size() - tdc == fNTDCs){
      minTDC = tdc;
      maxTDC = tdc + fNTDCs - 1;
      return;
    } 

    // For a given plane, look to see if the next 20 planes are empty. If not, this can be out start plane.
    int nEmpty = 0;
    for(unsigned int nextTDC = tdc + 1; nextTDC <= tdc + 20; ++nextTDC){
      if(tdcCharges[nextTDC] == 0.0) ++nEmpty;
    } 
    if(nEmpty < 5){
      minTDC = tdc;
      maxTDC = tdc + fNTDCs - 1;
      return;
    } 
  }

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

cvn::ViewVectorF cvn::CVNImageUtils::ConvertViewVecToViewVecF(cvn::ViewVector view){

  cvn::ViewVectorF newVec; 
  for(size_t w = 0; w < view.size(); ++w){
    std::vector<float> thisWire;
    for(size_t t = 0; t < view[w].size(); ++t){
      float chargeSC = static_cast<float>(view[w][t]);
      thisWire.push_back(chargeSC);
    }
    newVec.push_back(thisWire);
  }
  return newVec;
}

cvn::ImageVectorF cvn::CVNImageUtils::ConvertImageVecToImageVecF(cvn::ImageVector image){

  cvn::ImageVectorF newImage; 
  for(size_t w = 0; w < image.size(); ++w){
    cvn::ViewVectorF thisWire;
    for(size_t t = 0; t < image[w].size(); ++t){
      std::vector<float> thisTime;
      for(size_t v = 0; v < image[w][t].size(); ++v){
        float chargeSC = static_cast<float>(image[w][t][v]);
        thisTime.push_back(chargeSC);
      }
      thisWire.push_back(thisTime);
    }
    newImage.push_back(thisWire);
  }
  return newImage;
}

cvn::ImageVector cvn::CVNImageUtils::BuildImageVector(cvn::ViewVector v0, cvn::ViewVector v1, cvn::ViewVector v2){

  // Tensorflow wants things in the arrangement <wires, TDCs, views>
  cvn::ImageVector image;
  for(unsigned int w = 0; w < v0.size(); ++w){
    std::vector<std::vector<unsigned char> > wireVec;
    for(unsigned int t = 0; t < v0[0].size(); ++t){
      std::vector<unsigned char> timeVec;
      timeVec.push_back(v0[w][t]);
      timeVec.push_back(v1[w][t]);
      timeVec.push_back(v2[w][t]);
      wireVec.push_back(timeVec);
    } // Loop over tdcs
    image.push_back(wireVec);
  } // Loop over wires
  
  return image;

}

cvn::ImageVectorF cvn::CVNImageUtils::BuildImageVectorF(cvn::ViewVectorF v0, cvn::ViewVectorF v1, cvn::ViewVectorF v2){

  // Tensorflow wants things in the arrangement <wires, TDCs, views>
  cvn::ImageVectorF image;
  for(unsigned int w = 0; w < v0.size(); ++w){
    std::vector<std::vector<float> > wireVec;
    for(unsigned int t = 0; t < v0[0].size(); ++t){
      std::vector<float> timeVec;
      timeVec.push_back(v0[w][t]);
      timeVec.push_back(v1[w][t]);
      timeVec.push_back(v2[w][t]);
      wireVec.push_back(timeVec);
    } // Loop over tdcs
    image.push_back(wireVec);
  } // Loop over wires
  
  

  return image;
}


