////////////////////////////////////////////////////////////////////////
/// \file    CaffeNetHandler.h
/// \brief   CaffeNetHandler for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead   - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include "cetlib/getenv.h"

#include "dune/CVN/art/CaffeNetHandler.h"

#include "dune/CVN/func/CVNImageUtils.h"

namespace cvn
{

  /// This function just returns the thing you give it, but in the meanwhile
  /// it sets the glog level to what you pass in.
  int LogLevel(int level)
  {
    FLAGS_minloglevel = level;
    return level;
  }

// Delete this constructer... it isn't very constructive (sorry)
//  CaffeNetHandler::CaffeNetHandler( const std::string& libPath,
//                                    const std::string& deployProto,
//                                    const std::string& modelFile,
//				    const std::string& featureMap):
//  fLibPath(cet::getenv(libPath)+"/duneCVNNetwork/"),
//  fDeployProto(fLibPath + deployProto),
//  fModelFile  (fLibPath + modelFile),
//  fFeatureMap     (featureMap),
//  fLogLevel(LogLevel(google::WARNING)),
//  fNet(fDeployProto, caffe::TEST)
//  {
//    fNet.CopyTrainedLayersFrom(fModelFile);
//  }

  CaffeNetHandler::CaffeNetHandler(const fhicl::ParameterSet& pset):
    fLibPath(cet::getenv(pset.get<std::string>("LibPath", ""))+"/duneCVNNetwork/"),
    fDeployProto(fLibPath+pset.get<std::string>("DeployProto")),
    fModelFile  (fLibPath+pset.get<std::string>("ModelFile")),
    fFeatureMap (pset.get<std::string>("FeatureMap")),
    fLogLevel(LogLevel(google::WARNING)),
    fUseLogChargeScale(pset.get<bool>("ChargeLogScale")),
    fImageWires(pset.get<unsigned int>("NImageWires")),
    fImageTDCs(pset.get<unsigned int>("NImageTDCs")),
    fReverseViews(pset.get<std::vector<bool> >("ReverseViews")),
    fNet(fDeployProto, caffe::TEST)
  {
    fNet.CopyTrainedLayersFrom(fModelFile);
  }
  
  int CaffeNetHandler::NOutput() const
  {
    return fNet.output_blobs()[1]->shape(1);
  }
  
  int CaffeNetHandler::NFeatures() const
  {
    if(fFeatureMap==""){
      return 1;
    }
    return fNet.blob_by_name(fFeatureMap)->shape(1);
  }
  
  std::pair<const float*, const float*> CaffeNetHandler::Predict(const PixelMap& pm)
  {
    std::vector<caffe::Datum> datums;
    datums.push_back(PixelMapToDatum(pm));
    boost::dynamic_pointer_cast< caffe::MemoryDataLayer<float> >
      (fNet.layers()[0])->AddDatumVector(datums);
    std::vector<caffe::Blob<float>*> result = fNet.ForwardPrefilled();
    
    if(fFeatureMap==""){
      float features[]={-1};
      return std::make_pair<const float*, const float*>(result[1]->cpu_data(), features);
    }
    
    //Pull out the features stored in the final layer of network before the MLP layer that produces the final classifier:
    const boost::shared_ptr<caffe::Blob<float>> features = fNet.blob_by_name(fFeatureMap);
    
    return std::make_pair<const float*, const float*>(result[1]->cpu_data(), features->cpu_data());
  }
  
  caffe::Datum CaffeNetHandler::PixelMapToDatum(const PixelMap& pm)
  {

    unsigned int views = 3;

    // Use the CVNImageUtil class to produce the size of image we need
    CVNImageUtils imageUtils(fImageWires,fImageTDCs,views);
    // Make the pixel array that we will write into the datum
    std::vector<unsigned char> pixelArray(views*fImageWires*fImageTDCs,0);
    imageUtils.SetViewReversal(fReverseViews);
    imageUtils.ConvertPixelMapToPixelArray(pm,pixelArray);

    caffe::Datum datum;
    datum.set_channels(views);
    datum.set_height(fImageWires);
    datum.set_width(fImageTDCs);

    datum.set_data(pixelArray.data(), views*fImageWires*fImageTDCs);
    return datum;
    
   }


}
