#ifndef FLASH_MATCHER_HPP
#define FLASH_MATCHER_HPP

#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include "LowECluster.h"
#include "TGraphErrors.h"
#include <TH2.h>
#include <TTree.h>
#include <TMVA/TSpline1.h>
#include <canvas/Persistency/Common/Ptr.h>
#include <lardataobj/RecoBase/OpFlash.h>


// FUNCTIONS ------------------------------------------------------------------
// Compute the reconstructed energy from the charge, drift time and calibration parameters
// of the LikelihoodComputer class. The first three parameters belong to the class,
// while the last two are the charge of a cluster and the time difference between the candidate flash.
inline float give_me_Ereco(float calib_c, float calib_slope, float corr_lambda,
                           float dt, float charge){
  float q_corr = charge * exp(dt * corr_lambda);
  float E_reco = (q_corr - calib_c) / calib_slope;

  return E_reco;
}

// Returns the index of the TPC voxel corresponding to the input coordinates.
// The returned index is the one to be used to access the visibility map.
template <typename T>
int GetTPCIndex(T vertex_coor[3], TH1D* hgrid[3], std::vector<int>& total_to_tpc){
  int bin_x = hgrid[0]->FindBin(vertex_coor[0]);
  int bin_y = hgrid[1]->FindBin(vertex_coor[1]);
  int bin_z = hgrid[2]->FindBin(vertex_coor[2]);

  int tpc_index = total_to_tpc[bin_z + bin_y*hgrid[2]->GetNbinsX()
    + bin_x*hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()];

  return tpc_index; // Inside TPC
}


// From the hgridX histograms created by the PhotonVisibilityExpoort_module.cc, it returns
// a vector where the z-th + y-th * hgrid[2]->GetNbinsX() + x-th * hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()
// entry value is the number of the photoVisMap entry corresponding to the x,y,z coordinates.
template <typename T>
std::vector<int> GetCryoToTPCMap(TH1D* hgrid[3], T tpc_min[3], T tpc_max[3]){
  size_t n_total_entries = hgrid[0]->GetNbinsX() * hgrid[1]->GetNbinsX() * hgrid[2]->GetNbinsX();
  std::vector<int> total_to_tpc(n_total_entries, -1);
  int tpc_entry = 0;
  for (int ix=0; ix<hgrid[0]->GetNbinsX(); ix++){
    for (int iy=0; iy<hgrid[1]->GetNbinsX(); iy++){
      for (int iz=0; iz<hgrid[2]->GetNbinsX(); iz++){
        double x = hgrid[0]->GetBinCenter(ix);
        double y = hgrid[1]->GetBinCenter(iy);
        double z = hgrid[2]->GetBinCenter(iz);

        if (x > tpc_min[0] && x < tpc_max[0] && 
          y > tpc_min[1] && y < tpc_max[1] && 
          z > tpc_min[2] && z < tpc_max[2]) {
          total_to_tpc[iz + iy*hgrid[2]->GetNbinsX() + ix*hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()] = tpc_entry;
          tpc_entry++;
        }

      }
    }
  }
  return total_to_tpc;
}



class LikelihoodComputer{

public:

  float E_reco;
  float x_reco;
  float drift_velocity;
  TGraphErrors* g_logms;
  TGraphErrors* g_sigmas;
  TF1* f_lognormal;
  TF1* f_logms_trend;
  TF1* f_sigmas_trend;
  double trend_thr;
  TEfficiency* he_hit_prob;

  float GetNegativeLogLikelihoodMatch(const solar::LowECluster &tpc_cluster,
                                      const recob::OpFlash &pds_cluster);
  // Default constructor
  LikelihoodComputer() = default;
  // LikelihoodComputer constructor
  LikelihoodComputer(TString visibility_file_name,
                     std::string geom_identifier,
                     float LY_times_PDE,
                     TEfficiency* he_hit_prob,
                     TF1* f_lognormal,
                     TF1* f_logms_trend,
                     TF1* f_sigmas_trend,
                     float drift_velocity,
                     TGraphErrors* g_logms,
                     TGraphErrors* g_sigmas,
                     double trend_thr,
                     float calib_c, 
                     float calib_slope, 
                     float corr_lambda)
    : drift_velocity(drift_velocity), 
    g_logms(g_logms),
    g_sigmas(g_sigmas),
    f_lognormal(f_lognormal),
    f_logms_trend(f_logms_trend),
    f_sigmas_trend(f_sigmas_trend),
    trend_thr(trend_thr),
    he_hit_prob(he_hit_prob),
    visibility_file_name(visibility_file_name),
    geom_identifier(geom_identifier),
    LY_times_PDE(LY_times_PDE),
    calib_c(calib_c), 
    calib_slope(calib_slope), 
    corr_lambda(corr_lambda) {
    setprivatemembers();
  } // LikelihoodComputer constructor

  LikelihoodComputer CreateLikelihoodComputer(const double& fElectronScintYield, const std::string& fVisibilityFilename, double driftvelocity)
  {
    std::cout << fElectronScintYield << fVisibilityFilename << driftvelocity << std::endl;
    std::cout << "\n\n----------------------------\n FIX THE CREATION OF THE LIKELIHOOD COMPURER \n dunereco/LowEUtils/LikelihoodComputer.f \n.-----------------------------\n\n" << std::endl;
    LikelihoodComputer LC = LikelihoodComputer();
    return LC;
  }

private:
  // take in input
  TString visibility_file_name;
  std::string geom_identifier;
  size_t n_opdet;
  float LY_times_PDE;
  std::vector<std::vector<float>> opDet_visMapDirect;
  float calib_c, calib_slope, corr_lambda;

  // set inside the class
  TH1D* hgrid[3] = {nullptr};
  std::vector<int> cryo_to_tpc;
  float tpc_min[3];
  float tpc_max[3];
  float voxel_width;
  float anode_x;

  TMVA::TSpline1* g_he = nullptr;
  float xprob_max = 0.;

  void setprivatemembers();
}; // LikelihoodComputer

#endif // FLASH_MATCHER_HPP
