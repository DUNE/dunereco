#ifndef DUNERECO_FLASHMATCH_BDTMODEL_H
#define DUNERECO_FLASHMATCH_BDTMODEL_H

/**
 * \file BDTModel.h
 *
 * \brief Compact binary XGBoost decision tree model loader and evaluator.
 *
 * Decouples model weights from algorithmic C++ code, enabling zero-recompilation
 * model swapping and sub-millisecond initialization from lightweight .bin assets.
 *
 * \author Pierre Granger (pierre.granger@cern.ch)
 * \date 2026
 */

#include <vector>
#include <string>
#include <cstdint>
#include <cstddef>

namespace flashmatch {

struct BDTNode {
  int16_t feat = -1;      ///< Feature index (0..12, or -1 if leaf)
  float threshold = 0.0f; ///< Split threshold if decision node; leaf value if leaf
  int16_t left = -1;      ///< Left child index relative to tree start (-1 if leaf)
  int16_t right = -1;     ///< Right child index relative to tree start (-1 if leaf)
};

struct BDTTree {
  uint32_t node_offset = 0; ///< Offset in global node array
  uint16_t num_nodes = 0;   ///< Number of nodes in this tree
};

class BDTModel {
public:
  BDTModel() = default;

  /**
   * @brief Load model from binary format (.bin) file.
   *
   * Searches the provided path, FW_SEARCH_PATH, and FHICL_FILE_PATH.
   * @param filename Path or filename of binary model
   * @return true if successfully loaded, false otherwise
   */
  bool Load(const std::string& filename);

  bool IsLoaded() const { return fLoaded; }
  size_t NumTrees() const { return fNumTrees; }
  size_t NumFeatures() const { return fNumFeatures; }
  size_t TotalNodes() const { return fNodes.size(); }
  float BaseMargin() const { return fBaseMargin; }

  /// Predict raw log-odds: sum_t f_t(x) + base_margin
  float PredictRaw(const float* x) const;

  /// Predict classification probability: sigmoid(PredictRaw(x)) in [0, 1]
  float PredictProba(const float* x) const;

private:
  bool fLoaded = false;
  float fBaseMargin = 0.0f;
  uint32_t fNumFeatures = 0;
  uint32_t fNumTrees = 0;
  std::vector<BDTTree> fTrees;
  std::vector<BDTNode> fNodes;
};

} // namespace flashmatch

#endif // DUNERECO_FLASHMATCH_BDTMODEL_H
