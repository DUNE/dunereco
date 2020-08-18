////////////////////////////////////////////////////////////////////////
/// \file    RegCNNNumuHandler.cxx
/// \brief   RegCNNNumuHandler for numu energy estimation
/// \author  Wenjie Wu - wenjieww@uci.edu 
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include "cetlib/getenv.h"

#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dune/RegCNN/art/RegCNNNumuHandler.h"

#include "TH2D.h"
#include "TCanvas.h"

namespace cnn
{
  RegCNNNumuHandler::RegCNNNumuHandler(const fhicl::ParameterSet& pset):
  fTFHandlerContained (pset.get<fhicl::ParameterSet> ("TFNetHandlerContained")),
  fTFHandlerExiting   (pset.get<fhicl::ParameterSet> ("TFNetHandlerExiting"))
  {
  }
 
  std::vector<float> RegCNNNumuHandler::Predict(const RegPixelMap& pm, bool fLongestTrackContained)
  {
    std::vector<float> cnnResults;
    if (fLongestTrackContained) {
      cnnResults = fTFHandlerContained.Predict(pm);
    } else {
      cnnResults = fTFHandlerExiting.Predict(pm);
    }

    //std::cout << "Number of CNN result vectors " << cnnResults.size() << " with " << cnnResults[0].size() << " categories" << std::endl;

    //std::cout << "summary: ";
    //for(auto const v : cnnResults[0]){
    //  std::cout << v << ", ";
    //}
    //std::cout << std::endl;

    return cnnResults;
  }
} // end namespace cnn

////////////////////////////////////////////////////////////////////////
