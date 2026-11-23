#ifndef DUNERECO_FLASHMATCH_FLASHMATCHBDTALG_H
#define DUNERECO_FLASHMATCH_FLASHMATCHBDTALG_H

/**
 * \file FlashMatchBDTAlg.h
 *
 * \brief Pure C++ algorithmic engine for BDT-based charge-light matching in DUNE.
 *
 * Evaluates candidate flashes against a reconstructed charge position using
 * the 13-feature XGBoost BDT model loaded from an external compact binary asset,
 * selecting the most compatible flash and determining interaction T0.
 *
 * \author Pierre Granger (pierre.granger@cern.ch)
 * \date 2026
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstddef>

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "BDTModel.h"

namespace flashmatch {

struct MatchResult {
  bool matched = false;
  int flash_idx = -1;       ///< Index in original flash collection
  double flash_time = 0.0;  ///< Reconstructed flash time [us]
  double bdt_score = 0.0;   ///< BDT classification score in [0, 1]
  double pe = 0.0;          ///< Flash total PE
  double dy = 0.0;          ///< Centroid delta Y [cm]
  double dz = 0.0;          ///< Centroid delta Z [cm]
  double chi2 = 0.0;        ///< Spatial chi2
  int n_candidates = 0;     ///< Number of candidate flashes passing PE threshold
};

class FlashMatchBDTAlg {
public:
  struct Config {
    std::string weight_file = "FlashMatch/bdt_augmented.bin"; ///< Path or filename of binary model weights
    double pe_threshold = 100.0;                              ///< Minimum flash PE to consider as candidate
    double spatial_sigma = 62.64;                             ///< Spatial resolution scale sigma [cm]
    double min_bdt_score = 0.30;                              ///< Minimum BDT score to declare a positive match
    bool fallback_to_charge = true;                           ///< Allow fallback to charge centroid if no vertex
  };

  FlashMatchBDTAlg();
  explicit FlashMatchBDTAlg(const fhicl::ParameterSet& pset);
  explicit FlashMatchBDTAlg(const Config& config);

  void Configure(const fhicl::ParameterSet& pset);
  void Configure(const Config& config);

  /**
   * @brief Match a collection of flashes to a reconstructed charge position.
   *
   * @param flashes Input collection of recob::OpFlash
   * @param charge_y Y coordinate of reconstructed vertex or centroid [cm]
   * @param charge_z Z coordinate of reconstructed vertex or centroid [cm]
   * @param has_charge Whether the charge position is valid
   * @return MatchResult Output matching summary and best flash index
   */
  MatchResult Match(const std::vector<recob::OpFlash>& flashes,
                    double charge_y,
                    double charge_z,
                    bool has_charge = true) const;

  const Config& GetConfig() const { return fConfig; }
  const BDTModel& GetModel() const { return fModel; }

private:
  Config fConfig;
  BDTModel fModel;
};

} // namespace flashmatch

#endif // DUNERECO_FLASHMATCH_FLASHMATCHBDTALG_H
