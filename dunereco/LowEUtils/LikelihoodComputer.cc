#include "LikelihoodComputer.h"
#include "TLeaf.h"

void LikelihoodComputer::setprivatemembers()
{
  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");

  hgrid[0] = (TH1D*)visibility_file->Get("photovisAr/hgrid0")->Clone("hclone_hgrid0");
  hgrid[1] = (TH1D*)visibility_file->Get("photovisAr/hgrid1")->Clone("hclone_hgrid1");
  hgrid[2] = (TH1D*)visibility_file->Get("photovisAr/hgrid2")->Clone("hclone_hgrid2");
  voxel_width = 1.5*hgrid[0]->GetBinWidth(1); // HARD CODED: 1.5 factor to avoid edge effects in the visibility map

  TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
  Double_t coor_dim[3] = {0.};
  tDimensions->SetBranchAddress("dimension", coor_dim);
  tDimensions->GetEntry(0);
  tpc_min[0] = float(coor_dim[0]); tpc_min[1] = float(coor_dim[1]); tpc_min[2] = float(coor_dim[2]); // x,y,z
  tDimensions->GetEntry(1);
  tpc_max[0] = float(coor_dim[0]); tpc_max[1] = float(coor_dim[1]); tpc_max[2] = float(coor_dim[2]); // x,y,z
  std::cout << "TPC min: " << tpc_min[0] << ", " << tpc_min[1] << ", " << tpc_min[2] << std::endl;

  if (geom_identifier == "dune10kt")        anode_x = tpc_max[0]+tpc_min[0];
  else if (geom_identifier == "dunevd10kt") anode_x = tpc_max[0];
  else throw std::runtime_error("LikelihoodComputer: unknown geom_identifier " + geom_identifier + ". Supported values are: dune10kt, dunevd10kt.");
  std::cout << "Anode x = " << anode_x << std::endl;
  cryo_to_tpc = GetCryoToTPCMap(hgrid, tpc_min, tpc_max);

  TTree* photoVisMap = (TTree*)visibility_file->Get("photovisAr/photoVisMap");
  const size_t n_entriesmap = photoVisMap->GetEntries();
  TBranch* b_opDet_visDirect = photoVisMap->GetBranch("opDet_visDirect");
  n_opdet = b_opDet_visDirect->GetLeaf("opDet_visDirect")->GetLen();
  std::vector<float> opDet_visDirect(n_opdet, 0.);
  photoVisMap->SetBranchAddress("opDet_visDirect", opDet_visDirect.data());
  opDet_visMapDirect.resize(n_entriesmap, std::vector<float>(n_opdet, 0.));

  for (size_t VisMapEntry = 0; VisMapEntry < n_entriesmap; VisMapEntry++){
    photoVisMap->GetEntry(VisMapEntry);
    opDet_visMapDirect[VisMapEntry] = opDet_visDirect;
  }

  g_logms->Sort();  g_logms->SetBit(TGraph::kIsSortedX);
  g_sigmas->Sort(); g_sigmas->SetBit(TGraph::kIsSortedX);

  TGraph* g_he_graph = (TGraph*)he_hit_prob->CreateGraph();
  xprob_max = g_he_graph->GetX()[g_he_graph->GetN()-1];
  g_he = new TMVA::TSpline1("spline", g_he_graph);

  visibility_file->Close();
  return;
} // setprivatemembers

float LikelihoodComputer::GetNegativeLogLikelihoodMatch(const solar::LowECluster &tpc_cluster,
                                                        const recob::OpFlash &pds_cluster)
{
  float delta_time = tpc_cluster.getAverageTime() - float(pds_cluster.Time());
  E_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, delta_time, tpc_cluster.getTotalCharge());
  float reco_pe, exp_ph;
  x_reco = anode_x-delta_time*drift_velocity;
  float vertex_coor[3] = {x_reco, tpc_cluster.getY(), tpc_cluster.getZ()};
  vertex_coor[0] = std::max(vertex_coor[0], tpc_min[0]+voxel_width); vertex_coor[0] = std::min(vertex_coor[0], tpc_max[0]-voxel_width);
  vertex_coor[1] = std::max(vertex_coor[1], tpc_min[1]+voxel_width); vertex_coor[1] = std::min(vertex_coor[1], tpc_max[1]-voxel_width);
  vertex_coor[2] = std::max(vertex_coor[2], tpc_min[2]+voxel_width); vertex_coor[2] = std::min(vertex_coor[2], tpc_max[2]-voxel_width);
  int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);

  if (n_opdet != pds_cluster.PEs().size()){
    throw std::runtime_error("LikelihoodComputer: number of optical detectors in visibility file (" + std::to_string(n_opdet) + ") does not match the number of optical detectors in the OpFlash (" + std::to_string(pds_cluster.PEs().size()) + ").");
  }

  float n_hit = std::count_if(pds_cluster.PEs().begin(), pds_cluster.PEs().end(), [](float pe){ return pe > 0; });
  float weight = 1./(n_hit*n_hit);
  float term = 0.;
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

    if (reco_pe>0.0){
      if (exp_ph > trend_thr){
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

  return -log_likelihood*weight;
} // GetNegativeLogLikelihoodMatch
