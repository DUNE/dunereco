////////////////////////////////////////////////////////////////////////
/// \file    TFNetHandler.cxx
/// \brief   TFNetHandler for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead   - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include "cetlib/getenv.h"

#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dune/CVN/art/TFNetHandler.h"
#include "dune/CVN/func/CVNImageUtils.h"

#include "TH2D.h"
#include "TCanvas.h"

namespace cvn
{

  TFNetHandler::TFNetHandler(const fhicl::ParameterSet& pset):
    fLibPath(cet::getenv(pset.get<std::string>("LibPath", ""))),
    fTFProtoBuf  (fLibPath+"/"+pset.get<std::string>("TFProtoBuf")),
    fUseLogChargeScale(pset.get<bool>("ChargeLogScale")),
    fImageWires(pset.get<unsigned int>("NImageWires")),
    fImageTDCs(pset.get<unsigned int>("NImageTDCs")),
    fReverseViews(pset.get<std::vector<bool> >("ReverseViews"))
  {

    // Construct the TF Graph object. The empty vector {} is used since the protobuf
    // file gives the names of the output layer nodes
    mf::LogInfo("TFNetHandler") << "Loading network: " << fTFProtoBuf << std::endl;
    fTFGraph = tf::Graph::create(fTFProtoBuf.c_str(),{});
    if(!fTFGraph){
      art::Exception(art::errors::Unknown) << "Tensorflow model not found or incorrect";
    }

  }
  
  std::vector<float> TFNetHandler::Predict(const PixelMap& pm)
  {
   
    CVNImageUtils imageUtils;

    // Configure the image utility  
    imageUtils.SetViewReversal(fReverseViews);
    imageUtils.SetImageSize(fImageWires,fImageTDCs,3);
    imageUtils.SetLogScale(fUseLogChargeScale);

    ImageVectorF thisImage;
    imageUtils.ConvertPixelMapToImageVectorF(pm,thisImage);
    std::vector<ImageVectorF> vecForTF;

    vecForTF.push_back(thisImage);

    auto cvnResults = fTFGraph->run(vecForTF);

//    std::cout << "Number of CVN result vectors " << cvnResults.size() << " with " << cvnResults[0].size() << " categories" << std::endl;

    std::cout << "Classifier summary: ";
    for(auto const v : cvnResults[0]){
      std::cout << v << ", ";
    }
    std::cout << std::endl;

    // Leigh - test the new framework for easier access to multiple output architectures
    std::vector<std::vector<std::vector<float>>> newCVNResult = fTFGraph->runMultiOutput(vecForTF);
    std::vector<float> newResult = newCVNResult[0][0];
    std::cout << "New method classifier summary: ";
    for(auto const v : newResult){
      std::cout << v << ", ";
    }
    std::cout << std::endl;

    return cvnResults[0];
  }
 
  // The standard output has 13 elements, this function sums the convenient ones 
  std::vector<float> TFNetHandler::PredictFlavour(const PixelMap& pm){

    std::vector<float> fullResults = this->Predict(pm);

    std::vector<float> flavourResults;

    // First element is CC numu
    float sumNumu  = fullResults[0] + fullResults[1] + fullResults[2] + fullResults[3];
    // Then CC nue
    float sumNue   = fullResults[4] + fullResults[5] + fullResults[6] + fullResults[7];
    // Then CC nutau
    float sumNutau = fullResults[8] + fullResults[9] + fullResults[10] + fullResults[11];
    // End with NC
    float sumNC    = fullResults[12];

    flavourResults.push_back(sumNumu);
    flavourResults.push_back(sumNue);
    flavourResults.push_back(sumNutau);
    flavourResults.push_back(sumNC);

    return flavourResults;
  }

}

