//////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       EventImageData
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan (Dorota.Stefan@cern.ch), May 2016
//              Stripped from SPMultiTpcDump by Leigh Whitehead (leigh.howard.whitehead@cern.ch)
//
// Produces "global image" of events in multi-tpc detector, could be FD workspace or ProtoDUNE (SP).
// Projections of ADC are saved to 1-byte / pixel (TH1C) format, int-size (TH2I) projection of PDG
// truth and float-size (TH2F) projection of true charge deposits can be saved as well. We use this
// to dump deconvoluted ADC images for preparation of various classifiers and CNN based pattern reco.
//
// Moved from larreco/RecoAlg/ImagePatternAlgs since it became specific to DUNE geometries.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include "dune/CVN/adcutils/EventImageData.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

// C++ Includes
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH2C.h" // ADC map
#include "TH2I.h" // PDG+vertex info map
#include "TH2F.h" // deposit map


nnet::EventImageData::EventImageData(size_t w, size_t d, bool saveDep) :
    fVtxX(-9999), fVtxY(-9999),
    fProjX(-9999), fProjY(-9999),
    fSaveDep(saveDep)
{
    fAdc.resize(w, std::vector<float>(d, 0));
    if (saveDep) { fDeposit.resize(w, std::vector<float>(d, 0)); }
    fPdg.resize(w, std::vector<int>(d, 0));
}


void nnet::EventImageData::setProjXY(const TrainingDataAlg & dataAlg, float x, float y, size_t gw, bool flipw, size_t gd, bool flipd)
{
  if (flipw) { fProjX = gw + dataAlg.NWires() - x - 1; }
  else { fProjX = gw + x; }

  if (flipd) { fProjY = gd + dataAlg.NScaledDrifts() - y/dataAlg.DriftWindow() - 1; }
  else { fProjY = gd + y/dataAlg.DriftWindow(); }
}

void nnet::EventImageData::addTpc(const TrainingDataAlg & dataAlg, size_t gw, bool flipw, size_t gd, bool flipd)
{
  float zero = dataAlg.ZeroLevel();
  for (size_t w = 0; w < dataAlg.NWires(); ++w)
  {
      size_t drift_size = dataAlg.NScaledDrifts();
      std::vector<float> & dstAdc = fAdc[gw + w];
      std::vector<float> & dstDep = fDeposit[gw + w];
      std::vector<int> & dstPdg = fPdg[gw + w];
      const float* srcAdc = 0;
      const float* srcDep = 0;
      const int* srcPdg = 0;
      if (flipw)
      {
          srcAdc = dataAlg.wireData(dataAlg.NWires() - w - 1).data();
          srcDep = dataAlg.wireEdep(dataAlg.NWires() - w - 1).data();
          srcPdg = dataAlg.wirePdg(dataAlg.NWires() - w - 1).data();
      }
      else
      {
          srcAdc = dataAlg.wireData(w).data();
          srcDep = dataAlg.wireEdep(w).data();
          srcPdg = dataAlg.wirePdg(w).data();
      }
      if (flipd)
      {
          for (size_t d = 0; d < drift_size; ++d) { dstAdc[gd + d] += srcAdc[drift_size - d - 1] - zero; }
          if (fSaveDep) { for (size_t d = 0; d < drift_size; ++d) { dstDep[gd + d] += srcDep[drift_size - d - 1]; } }
          for (size_t d = 0; d < drift_size; ++d)
          {
              int code = srcPdg[drift_size - d - 1];
              int best_pdg = code & nnet::TrainingDataAlg::kPdgMask;
              int vtx_flags = (dstPdg[gd + d] | code) & nnet::TrainingDataAlg::kVtxMask;
              dstPdg[gd + d] = vtx_flags | best_pdg; // now just overwrite pdg and keep all vtx flags
              
              if (code & nnet::TrainingDataAlg::kNuPri) { fVtxX = gw + w; fVtxY = gd + d; } // dstAdc[gd + d] = 127; }
          }
      }
      else
      {
          for (size_t d = 0; d < drift_size; ++d) { dstAdc[gd + d] += srcAdc[d] - zero; }
          if (fSaveDep) { for (size_t d = 0; d < drift_size; ++d) { dstDep[gd + d] += srcDep[d]; } }
          for (size_t d = 0; d < drift_size; ++d)
          {
              int code = srcPdg[d];
              int best_pdg = code & nnet::TrainingDataAlg::kPdgMask;
              int vtx_flags = (dstPdg[gd + d] | code) & nnet::TrainingDataAlg::kVtxMask;
              dstPdg[gd + d] = vtx_flags | best_pdg; // now just overwrite pdg and keep all vtx flags

              if (code & nnet::TrainingDataAlg::kNuPri) { fVtxX = gw + w; fVtxY = gd + d; } // dstAdc[gd + d] = 127; }
          }
      }
  }
}

bool nnet::EventImageData::findCrop(size_t max_area_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const
{
  size_t max_cut = max_area_cut / 4;
  float adcThr = 10;

  w0 = 0;
  size_t cut = 0;
  while (w0 < fAdc.size())
  {
      for (auto const d : fAdc[w0]) { if (d > adcThr) cut++; }
      if (cut < max_cut) w0++;
      else break;
  }
  w1 = fAdc.size() - 1;
  cut = 0;
  while (w1 > w0)
  {
      for (auto const d : fAdc[w1]) { if (d > adcThr) cut++; }
      if (cut < max_cut) w1--;
      else break;
  }
  w1++;

  d0 = 0;
  cut = 0;
  while (d0 < fAdc.front().size())
  {
      for (size_t i = w0; i < w1; ++i) { if (fAdc[i][d0] > adcThr) cut++; }
      if (cut < max_cut) d0++;
      else break;
  }
  d1 = fAdc.front().size() - 1;
  cut = 0;
  while (d1 > d0)
  {
      for (size_t i = w0; i < w1; ++i) { if (fAdc[i][d1] > adcThr) cut++; }
      if (cut < max_cut) d1--;
      else break;
  }
  d1++;

  unsigned int margin = 32;
  if ((w1 - w0 > 8) && (d1 - d0 > 8))
  {
      if (w0 < margin) w0 = 0;
      else w0 -= margin;

      if (w1 > fAdc.size() - margin) w1 = fAdc.size();
      else w1 += margin;
      
      if (d0 < margin) d0 = 0;
      else d0 -= margin;
      
      if (d1 > fAdc.front().size() - margin) d1 = fAdc.front().size();
      else d1 += margin;
      
      return true;
  }
  else return false;
}

