////////////////////////////////////////////////////////////////////////
/// \file    TFNetHandler.cxx
/// \brief   TFNetHandler for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead   - leigh.howard.whitehead@cern.ch
///          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
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
    fTFGraph = tf::Graph::create(fTFProtoBuf.c_str(),{},pset.get<int>("NInputs"),pset.get<int>("NOutputs"));
    if(!fTFGraph){
      art::Exception(art::errors::Unknown) << "Tensorflow model not found or incorrect";
    }

  }
 
  // Check the network outputs
  bool check(const std::vector< std::vector< float > > & outputs)
  {
    if (outputs.size() == 1) return true;
    size_t aux = 0;
    for (size_t o = 0; o < outputs.size(); ++o)
    {   
        size_t aux2 = 0;

        for (size_t i = 0; i < outputs[o].size(); ++i)
            if (outputs[o][i] == 0.0 || outputs[o][i] == 1.0)
                aux2++;
        if (aux2 == outputs[o].size()) aux++;
    }
    return aux == outputs.size() ? false : true;        
  }

  // Fill outputs with value -3
  void fillEmpty(std::vector< std::vector< float > > & outputs)
  {
    for (size_t o = 0; o < outputs.size(); ++o)
    {
        for (size_t i = 0; i < outputs[o].size(); ++i)
            outputs[o][i] = -3.0;
    }
    return;
  }

  std::vector< std::vector<float> > TFNetHandler::Predict(const PixelMap& pm)
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

    std::vector< std::vector< std::vector< float > > > cvnResults; // shape(samples, #outputs, output_size)
    bool status = false;

    int counter = 0;

    do{ // do until it gets a correct result
        // std::cout << "Number of CVN result vectors " << cvnResults.size() << " with " << cvnResults[0].size() << " categories" << std::endl;
        cvnResults = fTFGraph->run(vecForTF);
        status = check(cvnResults[0]);
        //std::cout << "Status: " << status << std::endl;
        counter++;
        if(counter==10){
            std::cout << "Error, CVN never outputing a correct result. Filling result with zeros.";
            std::cout << std::endl;
            fillEmpty(cvnResults[0]);
            break;
        }
    }while(status == false);

    std::cout << "Classifier summary: ";
    std::cout << std::endl;
    int output_index = 0;
    for(auto const & output : cvnResults[0])
    {
      std::cout << "Output " << output_index++ << ": ";
      for(auto const v : output)
          std::cout << v << ", ";
      std::cout << std::endl;
    }
    std::cout << std::endl;

    // Leigh - test the new framework for easier access to multiple output architectures
//    std::vector<std::vector<std::vector<float>>> newCVNResult = fTFGraph->runMultiOutput(vecForTF);
//    std::vector<float> newResult = newCVNResult[0][0];
//    std::cout << "New method classifier summary: ";
//    for(auto const v : newResult){
//      std::cout << v << ", ";
//    }
//    std::cout << std::endl;

    return cvnResults[0];
  }

  /* 
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
  */

}

