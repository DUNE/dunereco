/**
 * \file BDTModel.cxx
 *
 * \brief Implementation of compact binary BDT model loader.
 *
 * \author Pierre Granger (pierre.granger@cern.ch)
 * \date 2026
 */

#include "BDTModel.h"
#include "cetlib/search_path.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <fstream>
#include <cmath>
#include <cstring>

namespace flashmatch {

bool BDTModel::Load(const std::string& filename)
{
  fLoaded = false;
  fTrees.clear();
  fNodes.clear();

  std::string full_path = filename;

  // Check if file exists directly
  std::ifstream test(full_path.c_str(), std::ios::binary);
  if (!test.is_open()) {
    // Try resolving via FW_SEARCH_PATH
    std::string resolved;
    cet::search_path sp("FW_SEARCH_PATH");
    if (sp.find_file(filename, resolved)) {
      full_path = resolved;
    } else {
      // Try resolving via FHICL_FILE_PATH
      cet::search_path fcl_sp("FHICL_FILE_PATH");
      if (fcl_sp.find_file(filename, resolved)) {
        full_path = resolved;
      }
    }
  } else {
    test.close();
  }

  std::ifstream in(full_path, std::ios::binary);
  if (!in.is_open()) {
    mf::LogError("BDTModel") << "Unable to open binary model file: " << filename
                             << " (attempted full path: " << full_path << ")";
    return false;
  }

  char magic[4];
  in.read(magic, 4);
  if (std::strncmp(magic, "BDT1", 4) != 0) {
    mf::LogError("BDTModel") << "Invalid magic header in " << full_path
                             << " (expected BDT1)";
    return false;
  }

  uint32_t version = 0;
  in.read(reinterpret_cast<char*>(&version), sizeof(uint32_t));
  if (version != 1) {
    mf::LogError("BDTModel") << "Unsupported model version: " << version;
    return false;
  }

  uint32_t total_nodes = 0;
  in.read(reinterpret_cast<char*>(&fNumTrees), sizeof(uint32_t));
  in.read(reinterpret_cast<char*>(&fNumFeatures), sizeof(uint32_t));
  in.read(reinterpret_cast<char*>(&total_nodes), sizeof(uint32_t));
  in.read(reinterpret_cast<char*>(&fBaseMargin), sizeof(float));

  fTrees.resize(fNumTrees);
  for (size_t t = 0; t < fNumTrees; ++t) {
    in.read(reinterpret_cast<char*>(&fTrees[t].node_offset), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&fTrees[t].num_nodes), sizeof(uint16_t));
  }

  fNodes.resize(total_nodes);
  for (size_t i = 0; i < total_nodes; ++i) {
    in.read(reinterpret_cast<char*>(&fNodes[i].feat), sizeof(int16_t));
    in.read(reinterpret_cast<char*>(&fNodes[i].threshold), sizeof(float));
    in.read(reinterpret_cast<char*>(&fNodes[i].left), sizeof(int16_t));
    in.read(reinterpret_cast<char*>(&fNodes[i].right), sizeof(int16_t));
  }

  if (!in.good()) {
    mf::LogError("BDTModel") << "Failed to read complete model data from " << full_path;
    return false;
  }

  fLoaded = true;
  mf::LogInfo("BDTModel") << "Loaded BDT model from " << full_path
                          << " (trees: " << fNumTrees
                          << ", features: " << fNumFeatures
                          << ", nodes: " << total_nodes
                          << ", base_margin: " << fBaseMargin << ")";
  return true;
}

float BDTModel::PredictRaw(const float* x) const
{
  if (!fLoaded) return 0.0f;
  float total = fBaseMargin;
  const auto* node_data = fNodes.data();
  for (size_t t = 0; t < fNumTrees; ++t) {
    const auto& tree = fTrees[t];
    int node_idx = 0;
    while (true) {
      const auto& node = node_data[tree.node_offset + node_idx];
      if (node.left == -1) {
        total += node.threshold;
        break;
      }
      if (x[node.feat] < node.threshold) {
        node_idx = node.left;
      } else {
        node_idx = node.right;
      }
    }
  }
  return total;
}

float BDTModel::PredictProba(const float* x) const
{
  float raw = PredictRaw(x);
  return 1.0f / (1.0f + std::exp(-raw));
}

} // namespace flashmatch
