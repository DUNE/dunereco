////////////////////////////////////////////////////////////////////////
/// \file    RegCNNVtxHandler.cxx
/// \brief   RegCNNVtxHandler 
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "canvas/Persistency/Common/FindManyP.h" 
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "dune/RegCNN/art/RegCNNVtxHandler.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cnn
{
  RegCNNVtxHandler::RegCNNVtxHandler(const fhicl::ParameterSet& pset):
    fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
    fTFHandler1st     (pset.get<fhicl::ParameterSet> ("TFNetHandler1st")),
    fTFHandler2nd     (pset.get<fhicl::ParameterSet> ("TFNetHandler2nd")),
    fTdcWidth         (pset.get<unsigned short>     ("TdcWidth")),
    fWireLength       (pset.get<unsigned short>     ("WireLength")),
    fTimeResolution   (pset.get<unsigned short>     ("TimeResolution")),
    fWireResolution   (pset.get<unsigned short>     ("WireResolution")),
    fGlobalWireMethod (pset.get<int>                ("GlobalWireMethod")),
//    fProngOnly        (pset.get<bool>               ("ProngOnly")),
    fProducer      (fWireLength, fWireResolution, fTdcWidth, fTimeResolution, fGlobalWireMethod, 0, 1)
  {
  }
  std::vector<float> RegCNNVtxHandler::GetVertex(detinfo::DetectorClocksData const& clockData,
                                                 detinfo::DetectorPropertiesData const& detProp,
                                                 art::Event &evt, const RegPixelMap &pixelmap){

      std::vector< art::Ptr< recob::Hit > > hitlist;
      auto hitListHandle = evt.getHandle< std::vector< recob::Hit > >(fHitsModuleLabel);
      if (hitListHandle)
        art::fill_ptr_vector(hitlist, hitListHandle);
      art::FindManyP<recob::Wire> fmwire(hitListHandle, evt, fHitsModuleLabel);

      /// Load in the pixel map
      // get vertex
      std::vector <float> Result(3, -99999);
      if (pixelmap.fInPM){
              std::vector<float> networkOutput(6, -99999);
              networkOutput = fTFHandler1st.Predict(pixelmap);
              FindGlobalVertices(pixelmap, networkOutput);
              std::vector <float> center_of_mass(7,0);
	      for (int ii = 0; ii < 6; ii++){
		      center_of_mass[ii] = networkOutput[ii];
	      }
	      networkOutput[1] += 7; networkOutput[5] += 7;
              networkOutput[3] -= 9;

      	      RegPixelMap pm;
              pm = fProducer.CreateMap(clockData, detProp, hitlist, fmwire, networkOutput);
      	      if (pm.fInPM){
		      center_of_mass[6] = (float)(pm.fTPC%4); // add TPC info
              	      Result = fTFHandler2nd.Predict(pm, center_of_mass);
	      }
      } // end of pixelmap
      return Result;
  } // end og GetVertex


  void RegCNNVtxHandler::FindGlobalVertices(const RegPixelMap& pm, std::vector<float> &outputs)
  {
    for (int ii = 0; ii < 3; ii++){
        float vtx_wire = outputs[ii*2+1];
        float vtx_tick = outputs[ii*2];
        float dwire = pm.fNWire/2-vtx_wire;
        float dtick = (pm.fNTdc/2-vtx_tick)*pm.fNTRes;
        float cm_wire = pm.fBound.fFirstWire[ii]+pm.fNWire/2-dwire;
        float cm_tdc  = pm.fBound.fFirstTDC[ii]+pm.fNTdc*pm.fNTRes/2-dtick;
        outputs[ii*2+1] = cm_wire;
        outputs[ii*2] = cm_tdc;
    }
  }
  void RegCNNVtxHandler::getCM(const RegPixelMap& pm, float* cm_list)
  {
    for (int ii = 0; ii < 3; ii++){
        float mean_wire = pm.fBound.fFirstWire[ii]+pm.fNWire/2;
        float mean_tdc  = pm.fBound.fFirstTDC[ii]+pm.fNTdc*pm.fNTRes/2;
        cm_list[2*ii] = mean_tdc;
        cm_list[2*ii+1] = mean_wire;
    }
  }



}
