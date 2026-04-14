#ifndef FLASH_MATCHER_HPP
#define FLASH_MATCHER_HPP

#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include "LowECluster.h"
#include "TGraphErrors.h"
#include "TLeaf.h"
#include <TH2.h>
#include <TTree.h>
#include <TMVA/TSpline1.h>
#include <canvas/Persistency/Common/Ptr.h>
#include <filesystem>
#include <lardataobj/RecoBase/OpFlash.h>

// using namespace TMVA;

// FUNCTIONS ------------------------------------------------------------------
inline float give_me_Ereco(float calib_c, float calib_slope, float corr_lambda,
                           float dt, float charge){
  
  float q_corr = charge * exp(dt * corr_lambda);
  float E_reco = (q_corr - calib_c) / calib_slope;
  
  return E_reco;
}

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

    float GetLogLikelihoodMatch(const solar::LowECluster &tpc_cluster,
                                const recob::OpFlash &pds_cluster)
    {
      double clusterTime = tpc_cluster.getAverageTime();
      double flashTime = pds_cluster.Time();
      float delta_time = float(clusterTime - flashTime);
      E_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, delta_time, tpc_cluster.getTotalCharge());
      float reco_pe, exp_ph;
      x_reco = delta_time*drift_velocity;
      float vertex_coor[3] = {x_reco, tpc_cluster.getY(), tpc_cluster.getZ()};
      vertex_coor[0] = std::max(vertex_coor[0], tpc_min[0]+15); vertex_coor[0] = std::min(vertex_coor[0], tpc_max[0]-15); // HARD CODE
      vertex_coor[1] = std::max(vertex_coor[1], tpc_min[1]+15); vertex_coor[1] = std::min(vertex_coor[1], tpc_max[1]-15);
      vertex_coor[2] = std::max(vertex_coor[2], tpc_min[2]+15); vertex_coor[2] = std::min(vertex_coor[2], tpc_max[2]-15);
      int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);

      if (n_opdet != pds_cluster.PEs().size()){
        throw std::runtime_error("LikelihoodComputer: number of optical detectors in visibility file (" + std::to_string(n_opdet) + ") does not match the number of optical detectors in the OpFlash (" + std::to_string(pds_cluster.PEs().size()) + ").");
      }

      float log_likelihood = 0.;
      for (size_t idx_opdet=0; idx_opdet<pds_cluster.PEs().size(); idx_opdet++){
        float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];

        exp_ph = E_reco*LY_times_PDE*voxel_vis;
        if(exp_ph==0) exp_ph = E_reco*LY_times_PDE*1.e-15;
        float P_hit_mu = (exp_ph<xprob_max) ? g_he->Eval(exp_ph) : 1.; // Avoid weird extrapolation where we
        //                                                             // have no points in the efficiency graph,
        //                                                             // and just assume that the hit probability is 1 for very high expected PE.
        if (P_hit_mu <= 0.) P_hit_mu = 1.e-4;
        if (P_hit_mu >= 1.) P_hit_mu = 1. - 1.e-4;

        reco_pe = pds_cluster.PEs().at(idx_opdet);
        float term = 0.;

        if (reco_pe>0.0){
          if (reco_pe > trend_thr){
            f_lognormal->SetParameters(f_logms_trend->Eval(exp_ph), f_sigmas_trend->Eval(exp_ph));
            term = P_hit_mu*f_lognormal->Eval(reco_pe);
            log_likelihood += log(term);
          } else {
            f_lognormal->SetParameters(g_logms->Eval(exp_ph), g_sigmas->Eval(exp_ph));
            term = P_hit_mu*f_lognormal->Eval(reco_pe);
            log_likelihood += log(term);
          }
        } else {
          term = 1 - P_hit_mu;
          log_likelihood += log(term);
        }
      }

      return log_likelihood;
    } // GetLogLikelihoodMatch

    // Default constructor
    LikelihoodComputer() = default;
    // LikelihoodComputer constructor
    LikelihoodComputer(TString visibility_file_name,
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
    size_t n_opdet;
    float LY_times_PDE;
    std::vector<std::vector<float>> opDet_visMapDirect;
    float calib_c, calib_slope, corr_lambda;

    // set inside the class
    TH1D* hgrid[3] = {nullptr};
    std::vector<int> cryo_to_tpc;
    float tpc_min[3];
    float tpc_max[3];

    TMVA::TSpline1* g_he = nullptr;
    float xprob_max = 0.;


    void setprivatemembers(){
      TFile* visibility_file = TFile::Open(visibility_file_name, "READ");

      hgrid[0] = (TH1D*)visibility_file->Get("photovisAr/hgrid0")->Clone("hclone_hgrid0");
      hgrid[1] = (TH1D*)visibility_file->Get("photovisAr/hgrid1")->Clone("hclone_hgrid1");
      hgrid[2] = (TH1D*)visibility_file->Get("photovisAr/hgrid2")->Clone("hclone_hgrid2");

      TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
      Double_t coor_dim[3] = {0.};
      tDimensions->SetBranchAddress("dimension", coor_dim);
      tDimensions->GetEntry(0);
      tpc_min[0] = float(coor_dim[0]); tpc_min[1] = float(coor_dim[1]); tpc_min[2] = float(coor_dim[2]); // x,y,z
      tDimensions->GetEntry(1);
      tpc_max[0] = float(coor_dim[0]); tpc_max[1] = float(coor_dim[1]); tpc_max[2] = float(coor_dim[2]); // x,y,z
      std::cout << "TPC min: " << tpc_min[0] << ", " << tpc_min[1] << ", " << tpc_min[2] << std::endl;
      cryo_to_tpc = GetCryoToTPCMap(hgrid, tpc_min, tpc_max);

      TTree* photoVisMap = (TTree*)visibility_file->Get("photovisAr/photoVisMap");
      const size_t n_entriesmap = photoVisMap->GetEntries();
      TBranch* b_opDet_visDirect = photoVisMap->GetBranch("opDet_visDirect");
      n_opdet = b_opDet_visDirect->GetLeaf("opDet_visDirect")->GetLen();
      std::vector<float> opDet_visDirect(n_opdet, 0.);
      photoVisMap->SetBranchAddress("opDet_visDirect", opDet_visDirect.data());
      std::cout << "resizing" << std::endl;
      std::cerr << "resizing" << std::endl;
      opDet_visMapDirect.resize(n_entriesmap, std::vector<float>(n_opdet, 0.));
      std::cout << "resiz" << std::endl;
      std::cerr << "resiz" << std::endl;

      for (size_t VisMapEntry = 0; VisMapEntry < n_entriesmap; VisMapEntry++){
        photoVisMap->GetEntry(VisMapEntry);
        opDet_visMapDirect[VisMapEntry] = opDet_visDirect;
      }
      std::cout << "riz" << std::endl;
      std::cerr << "riz" << std::endl;

      TGraph* g_he_graph = (TGraph*)he_hit_prob->CreateGraph();
      xprob_max = g_he_graph->GetX()[g_he_graph->GetN()-1];
      g_he = new TMVA::TSpline1("spline", g_he_graph);

      visibility_file->Close();
      return;
    } // 

}; // LikelihoodComputer

#endif // FLASH_MATCHER_HPP
