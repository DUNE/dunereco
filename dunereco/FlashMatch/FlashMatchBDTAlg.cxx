/**
 * \file FlashMatchBDTAlg.cxx
 *
 * \brief Pure C++ algorithmic engine for BDT-based charge-light matching in DUNE.
 *
 * \author Pierre Granger (pierre.granger@cern.ch)
 * \date 2026
 */

#include "FlashMatchBDTAlg.h"
#include <numeric>
#include <limits>
#include <cmath>

namespace flashmatch {

FlashMatchBDTAlg::FlashMatchBDTAlg()
  : fConfig{}
{}

FlashMatchBDTAlg::FlashMatchBDTAlg(const fhicl::ParameterSet& pset)
{
  Configure(pset);
}

FlashMatchBDTAlg::FlashMatchBDTAlg(const Config& config)
  : fConfig(config)
{
  Configure(config);
}

void FlashMatchBDTAlg::Configure(const fhicl::ParameterSet& pset)
{
  fConfig.weight_file = pset.get<std::string>("WeightFile", "FlashMatch/bdt_augmented.bin");
  fConfig.pe_threshold = pset.get<double>("PEThreshold", 100.0);
  fConfig.spatial_sigma = pset.get<double>("SpatialSigma", 62.64);
  fConfig.min_bdt_score = pset.get<double>("MinBDTScore", 0.30);
  fConfig.fallback_to_charge = pset.get<bool>("FallbackToCharge", true);

  fModel.Load(fConfig.weight_file);
}

void FlashMatchBDTAlg::Configure(const Config& config)
{
  fConfig = config;
  fModel.Load(fConfig.weight_file);
}

MatchResult FlashMatchBDTAlg::Match(const std::vector<recob::OpFlash>& flashes,
                                    double charge_y,
                                    double charge_z,
                                    bool has_charge) const
{
  MatchResult result;

  if (!fModel.IsLoaded() || !has_charge || !std::isfinite(charge_y) || !std::isfinite(charge_z)) {
    return result;
  }

  // Identify candidates passing PE threshold
  struct Candidate {
    int orig_idx;
    double pe;
    float rank;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(flashes.size());

  double pe_max = 0.0;
  for (size_t i = 0; i < flashes.size(); ++i) {
    const auto& f = flashes[i];
    if (f.TotalPE() >= fConfig.pe_threshold) {
      candidates.push_back({static_cast<int>(i), f.TotalPE(), 0.0f});
      if (f.TotalPE() > pe_max) {
        pe_max = f.TotalPE();
      }
    }
  }

  const int n_cand = static_cast<int>(candidates.size());
  result.n_candidates = n_cand;

  if (n_cand == 0) {
    return result;
  }

  // Calculate PE ranks (0 = highest PE, 1 = 2nd highest, ...)
  std::vector<size_t> order(n_cand);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&candidates](size_t a, size_t b) {
    return candidates[a].pe > candidates[b].pe;
  });

  for (size_t r = 0; r < order.size(); ++r) {
    candidates[order[r]].rank = static_cast<float>(r);
  }

  float best_score = -1.0f;
  int best_cand_idx = -1;
  double best_dy = 0.0;
  double best_dz = 0.0;
  double best_chi2 = 0.0;

  static constexpr size_t NUM_FEATURES = 13;
  float x[NUM_FEATURES];

  for (int c = 0; c < n_cand; ++c) {
    const auto& cand = candidates[c];
    const auto& f = flashes[cand.orig_idx];

    const double pe = cand.pe;
    const double dy = f.YCenter() - charge_y;
    const double dz = f.ZCenter() - charge_z;
    const double abs_dy = std::abs(dy);
    const double abs_dz = std::abs(dz);
    const double dy_norm = dy / fConfig.spatial_sigma;
    const double dz_norm = dz / fConfig.spatial_sigma;
    const double chi2 = dy_norm * dy_norm + dz_norm * dz_norm;
    const double chi = std::sqrt(chi2);
    const double yw = f.YWidth();
    const double zw = f.ZWidth();
    const double width_sum = yw + zw;
    const double lpe = std::log10(std::max(pe, 1.0));
    const double pe_ratio = pe / std::max(pe_max, 1.0);
    const double log10_pe_ratio = std::log10(std::max(pe_ratio, 1e-6));

    x[0] = static_cast<float>(lpe);
    x[1] = static_cast<float>(dy);
    x[2] = static_cast<float>(dz);
    x[3] = static_cast<float>(abs_dy);
    x[4] = static_cast<float>(abs_dz);
    x[5] = static_cast<float>(chi2);
    x[6] = static_cast<float>(chi);
    x[7] = static_cast<float>(yw);
    x[8] = static_cast<float>(zw);
    x[9] = static_cast<float>(width_sum);
    x[10] = cand.rank;
    x[11] = static_cast<float>(log10_pe_ratio);
    x[12] = static_cast<float>(n_cand);

    // Safeguard against NaN/inf (mimics np.nan_to_num)
    for (size_t fi = 0; fi < NUM_FEATURES; ++fi) {
      if (!std::isfinite(x[fi])) {
        x[fi] = 0.0f;
      }
    }

    const float score = fModel.PredictProba(x);

    if (score > best_score) {
      best_score = score;
      best_cand_idx = c;
      best_dy = dy;
      best_dz = dz;
      best_chi2 = chi2;
    }
  }

  if (best_cand_idx >= 0) {
    const auto& best_cand = candidates[best_cand_idx];
    const auto& best_flash = flashes[best_cand.orig_idx];

    result.flash_idx = best_cand.orig_idx;
    result.flash_time = best_flash.Time();
    result.bdt_score = static_cast<double>(best_score);
    result.pe = best_cand.pe;
    result.dy = best_dy;
    result.dz = best_dz;
    result.chi2 = best_chi2;

    if (best_score >= fConfig.min_bdt_score) {
      result.matched = true;
    }
  }

  return result;
}

} // namespace flashmatch
