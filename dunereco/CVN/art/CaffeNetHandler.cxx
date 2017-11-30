////////////////////////////////////////////////////////////////////////
/// \file    CaffeNetHandler.h
/// \brief   CaffeNetHandler for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include "cetlib/getenv.h"

#include "dune/CVN/art/CaffeNetHandler.h"

namespace cvn
{

  /// This function just returns the thing you give it, but in the meanwhile
  /// it sets the glog level to what you pass in.
  int LogLevel(int level)
  {
    FLAGS_minloglevel = level;
    return level;
  }


  CaffeNetHandler::CaffeNetHandler( const std::string& libPath,
                                    const std::string& deployProto,
                                    const std::string& modelFile,
				    const std::string& featureMap):
  fLibPath(cet::getenv(libPath)+"/duneCVNNetwork/"),
  fDeployProto(fLibPath + deployProto),
  fModelFile  (fLibPath + modelFile),
  fFeatureMap     (featureMap),
  fLogLevel(LogLevel(google::WARNING)),
  fNet(fDeployProto, caffe::TEST)
  {
    fNet.CopyTrainedLayersFrom(fModelFile);
  }

  CaffeNetHandler::CaffeNetHandler(const fhicl::ParameterSet& pset):
    CaffeNetHandler(
                    pset.get<std::string>("LibPath", ""),
                    pset.get<std::string>("DeployProto"),
                    pset.get<std::string>("ModelFile"),
                    pset.get<std::string>("FeatureMap")
                    )
  {}
  
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
  
  char ConvertToChar(float n)
  {
    float peCorrChunk;
    // This value is 2.0 / (PE/GeV)
    peCorrChunk = 1000.0/255.0;
    
    float truncateCorr = ceil(n/peCorrChunk);
    if (truncateCorr > 255){ 
      truncateCorr= 255.0;
    }
    return (char)truncateCorr;
  }
  
  caffe::Datum PixelMapToDatum(const PixelMap& pm)
  {


    // caffe::Datum datum;
    // char* pixels = NULL;
    // int channels(3), planes(0), cells(0);
    // datum.set_channels(channels);
    // planes = pm.fNWire;
    // cells  = pm.fNTdc;
    // datum.set_height(planes);
    // datum.set_width(cells);
    // pixels = new char[channels*planes*cells];
    // for (int iChan = 0; iChan < channels; ++iChan)
    // {
    //   for (int iPlane = 0; iPlane < planes; ++iPlane)
    //   {
    //     for (int iCell = 0; iCell < cells; ++iCell)
    //     {
    //       int i = iCell + cells*(iPlane + planes*iChan);
    //       float val =0.;
    //       if(iChan == 0 ){
    //         val=pm.fPEX.at(iCell + cells*iPlane);
    //       }
    //       if(iChan == 1 ){
    //         val=pm.fPEY.at(iCell + cells*iPlane);
    //       }
    //       if(iChan == 2 ){
    //         val=pm.fPEZ.at(iCell + cells*iPlane);
    //       }
    //       char pix = ConvertToChar(val);
    //       pixels[i] = pix;
    //     }
    //   }
    // }
    // datum.set_data(pixels, channels*planes*cells);
    // delete[] pixels;
    // return datum;
  
    caffe::Datum datum;
    char* pixels = NULL;
    int channels(3), planes(0), cells(0);
    datum.set_channels(channels);
    planes = 500; //pm.fNWire; TODO Make fcl parameter
    cells  = pm.fNTdc;
    datum.set_height(planes);
    datum.set_width(cells);
    pixels = new char[channels*planes*cells];
    for (int iChan = 0; iChan < channels; ++iChan)
    {
      for (int iPlane = 0; iPlane < planes; ++iPlane)
      {
        for (int iCell = 0; iCell < cells; ++iCell)
        {
          int i = iCell + cells*(iPlane + planes*iChan);
          float val =0.;
          if(iChan == 0 ){
            val=pm.fPEX.at(iCell + cells*iPlane);
          }
          if(iChan == 2 ){
            val=pm.fPEZ.at(iCell + cells*iPlane);
          }
          char pix = ConvertToChar(val);
          pixels[i] = pix;
        }
      }
    }

    //work out where to flip the second induction view, so they point the same way
    int maxPlaneWithCharge=0;
    for (int iChan = 0; iChan < channels; ++iChan)
    {
      for (int iPlane = 0; iPlane < planes; ++iPlane)
        {
          for (int iCell = 0; iCell < cells; ++iCell)
            {
              float val =0.;
              if(iChan == 1 ){
                val=pm.fPEY.at(iCell + cells*iPlane);
                if(val>0){
                  maxPlaneWithCharge=iPlane;
                }
              }
            }
        }
    }
    int newMapEdge=planes;
    if(maxPlaneWithCharge>planes){
      newMapEdge=maxPlaneWithCharge+1;
    }
    
    //fill the second induction view
    for (int iChan = 0; iChan < channels; ++iChan)
      {
        for (int iPlane = newMapEdge; iPlane > (newMapEdge-planes); --iPlane)
          {
            for (int iCell = 0; iCell < cells; ++iCell)
              {
                int planeCount=planes-iPlane;
                int i = iCell + cells*(planeCount + planes*iChan);
                float val =0.;
                if(iChan == 1 ){
                  val=pm.fPEY.at(iCell + cells*iPlane);
                  char pix = ConvertToChar(val);
                  pixels[i] = pix;
                }
              }
          }
      }

    datum.set_data(pixels, channels*planes*cells);
    delete[] pixels;
    return datum;
   }


}
