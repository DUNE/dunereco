////////////////////////////////////////////////////////////////////////
/// \file    RegCNNVtxNetHandler.h
/// \brief   Reg CNN Vertex Handler
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_VTXHANDLER_H
#define REGCNN_VTXHANDLER_H

#include <vector>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "dune/RegCNN/art/RegPixelMapProducer.h"
#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/art/TFRegNetHandler.h"

namespace cnn
{
  /// RegCNNVtxHandler, basic output of CNN neural net
  class RegCNNVtxHandler
  {
  public:
    RegCNNVtxHandler(const fhicl::ParameterSet& pset);

    std::vector<float> GetVertex(detinfo::DetectorClocksData const& clockData,
                                 detinfo::DetectorPropertiesData const& detProp,
                                 art::Event& evt, const RegPixelMap &pixelmap);

  private:
    void FindGlobalVertices(const RegPixelMap &pm, std::vector<float> &outputs);
    void getCM(const RegPixelMap& pm, float* cm_list);

    std::string fHitsModuleLabel;
    cnn::TFRegNetHandler fTFHandler1st;
    cnn::TFRegNetHandler fTFHandler2nd;
    unsigned int fTdcWidth;        
    unsigned int fWireLength;     
    unsigned int fTimeResolution;  
    unsigned int fWireResolution;  
    unsigned int fGlobalWireMethod;
//    bool fProngOnly;
    RegPixelMapProducer fProducer;

  };

}

#endif  // REGCNN_VTXHANDLER_H
