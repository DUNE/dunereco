/**
 * @file PointResTree_module.cc
 * @author James Shen (jieran.shen@duke.edu)
 * @brief Retrieve and calculate relevant infromation from the art reco file.
 * Save all information in a tree.
 * @version 2.1
 * @date 2021-11-29
 * Specifically, this module does the following:
 *      - Retrieve relevant MC information for the neutrino and the electron.
 *      - Resolve track directional ambiguity via daughter flipping.
 *      - Reconstruct electron energy using charge information from the planes,
 * drift corrected by photon flashes.
 *
 */

#include <RtypesCore.h>
#ifndef PointResTree_H
#define PointResTree_H 1
#define DEBUG_MSG 0
// Framework includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes

// #include "dune/AnaUtils/DUNEAnaHitUtils.h"
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "Math/Functor.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPolyLine3D.h"
#include "TPrincipal.h"
#include "TTree.h"
#include "TVector3.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

using XYZVector = TVector3;
using track_loc =
    std::pair<std::string, size_t>; // <ModuleLabel, idx> of a track.
namespace dune {
/**
 * @brief Faster version of ParticleInventoryService, used to the generator
 * associated with a reco track.
 *
 */
struct genFinder {
private:
  typedef std::pair<int, std::string> track_id_to_string;
  std::vector<track_id_to_string> track_id_map;
  std::set<std::string> generator_names;
  bool isSorted = false;

public:
  void sort_now() {
    std::sort(this->track_id_map.begin(), this->track_id_map.end(),
              [](const auto &a, const auto &b) { return (a.first < b.first); });
    isSorted = true;
  }

  void add(const int &track_id, const std::string &gname) {
    this->track_id_map.push_back(std::make_pair(track_id, gname));
    generator_names.emplace(gname);
    isSorted = false;
  }

  bool has_gen(std::string gname) {
    return static_cast<bool>(generator_names.count(gname));
  };

  std::string get_gen(int tid) {
    if (!isSorted)
      this->sort_now();
    return std::lower_bound(
               track_id_map.begin(), track_id_map.end(), tid,
               [](const auto &a, const auto &b) { return (a.first < b); })
        ->second;
  };
};

class PointResTree : public art::EDAnalyzer {
public:
  PointResTree(fhicl::ParameterSet const &);
  /**
   * @brief Read in fcl parameters, configure them to proper class variables
   */
  void reconfigure(fhicl::ParameterSet const &);
  void beginJob() override;
  void analyze(art::Event const &) override;
  void endJob() override;

private:
  // NOT IN TREE:
  genFinder *gf;
  track_loc primary_trk_loc; ///< Keeping track of primary track's id in its
                             ///< producer.
  std::vector<recob::Hit> hits_U, hits_V, hits_Z;
  bool has_sp_vertex;

  // operation flags
  bool fCollectTracks;
  // art module labels
  std::string fGeneratorLabel;
  std::string fSimulationLabel;
  std::string fSimChannelModuleLabel;
  std::vector<std::string> fTrackModuleLabel;
  std::string fHitsModuleLabel;
  std::string fOpFlashLabel;
  std::string fHitFdModuleLabel;
  std::string fSpacePointModuleLabel;
  std::string fHitToSpacePointLabel;

  // histograms
  TH1F *vertex_reco_distance; ///< Distance between the SP vertex and truth
                              ///< position

  // Tree & branches
  // Note: Distance are in cm, time is in ms.
  TTree *tr;
  XYZVector truth_nu_dir, truth_e_dir; ///< Truth directions
  Double_t truth_nu_en, truth_e_en;    ///< Truth energies
  Double_t truth_en_deposition,
      truth_charge_deposition; ///< G4 deposition information
  Int_t nu_pdg;                ///< Flavor of the neutrino
  XYZVector truth_e_position;  ///< Truth vertex

  Int_t NParticles, NHits, NTrks; ///< Talley of important products.
  Double_t charge_U, charge_V,
      charge_Z; ///< Total charge on each plane (ADC).
  std::vector<Double_t> charge_Z_dist;
  ///< Vector of collection plane charge over different distance cut from
  ///< veretx.
  Double_t charge_corrected; ///< Charge after correcting for drift.
  Double_t cut_charge_loss;  ///< Portion of the charge cut by the distance cut.
  Double_t drift_time;       ///< Drift time (reconstructed from photons).

  XYZVector reco_e_dir;      ///< Reconstructed electron direction.
  Double_t reco_e_en;        ///< Reconstructed electron energy.
  XYZVector reco_e_position; ///< Reconstructed vertex.

  Double_t reco_distance; ///< Distance between reconstructed and truth vertex.
  Int_t primary_trk_id;   ///< ID of the primary tracks in the parallel array.

  std::vector<Double_t> trk_length; ///< Length of all vallid tracks.
  std::vector<Double_t> trk_charge; ///< Charge of all vallid tracks.
  std::vector<Double_t> trk_start_x, trk_start_y,
      trk_start_z; ///< Starting Position of all vallid tracks.
  std::vector<Double_t> trk_end_x, trk_end_y,
      trk_end_z; ///< Ending Position of all vallid tracks.
  std::vector<Double_t> trk_start_dir_x, trk_start_dir_y, trk_start_dir_z;
  ///< Starting Direction of all vallid tracks.
  std::vector<Double_t> trk_end_dir_x, trk_end_dir_y, trk_end_dir_z;
  ///< Starting Direction of all vallid tracks.

  // helper functions
  /**
   * @brief Rest all tree variables to bogus values.
   */
  void reset_variables();

  /**
   * @brief Write Neutrino and Electron truths to the tree.
   *
   * @param event art event.
   */
  void writeMCTruths(art::Event const &event);

  /**
   * @brief Given accumulated charges on the U, V, and Z planes, does drift
   * correction and reco energy of the electron.
   * @param event art event.
   */
  void calculate_electron_energy(art::Event const &event);

  /**
   * @brief Reads the track direction information, calculate the direction of
   * the primary track. Most importantly, this function carries out daughter
   * flipping.
   */
  void calculate_electron_direction();

  /**
   * @brief Determine if a hit should be included in the energy sum.
   * If there are 3D reco information (SpacePoints) associated with the hit,
   * use that as the track location. If there are no spacepoints, carry out a
   * 2D reconstruction using drift time and the center of the wire segment.
   * @param event art event.
   * @param hit   hit currently working on.
   * @param spacepoints   spacepoints associated with the hit.
   * @return Bool_t True if should be cut, false if shoudl be included.
   */
  Bool_t distance_cut(art::Event const &event, recob::Hit const &hit,
                      std::vector<art::Ptr<recob::SpacePoint>> spacepoints);

  /**
   * @brief Determine the vertex of the SN event using reconstructed
   * spacepoints Goes through all the spacepoints, select the top 10
   * spacepoints ranked by associated z-view hits' amplitude. It then tries to
   * select one of the ten, hopefully in the vinicinity of the truth vertex.
   * @param event art event.
   * @return XYZVector result, positional vector of the reco vertex.
   */
  XYZVector find_vertex_sp(art::Event const &event);
  /* Function to write hit as a distance of truth distance. */
  void write_multi_distance(art::Event const &, recob::Hit const &,
                            std::vector<art::Ptr<recob::SpacePoint>>);

}; // class PointResTree

} // namespace dune

#endif

namespace dune {
DEFINE_ART_MODULE(PointResTree)
}

dune::PointResTree::PointResTree(fhicl::ParameterSet const &parameterSet)
    : EDAnalyzer(parameterSet) {
  this->reconfigure(parameterSet);
}

void dune::PointResTree::reconfigure(fhicl::ParameterSet const &parameterSet) {
  fGeneratorLabel = parameterSet.get<std::string>(
      "GeneratorLabel"); // label for MCTruth block
  fSimulationLabel = parameterSet.get<std::string>(
      "SimulationLabel"); // label for MCParticle block, probably largeant
  fSimChannelModuleLabel =
      parameterSet.get<std::string>("SimChannelModuleLabel");
  fTrackModuleLabel =
      parameterSet.get<std::vector<std::string>>("TrackModuleLabel");
  fHitsModuleLabel = parameterSet.get<std::string>("HitsModuleLabel");
  fOpFlashLabel = parameterSet.get<std::string>("OpFlashLabel");
  fHitToSpacePointLabel = parameterSet.get<std::string>("HitToSpacePointLabel");
  fCollectTracks = parameterSet.get<bool>("CollectTracks");
  fSpacePointModuleLabel =
      parameterSet.get<std::string>("SpacePointModuleLabel");
}

void dune::PointResTree::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  vertex_reco_distance = tfs->make<TH1F>(
      "vertex_reco_distance",
      "Distance between reconstructed vertex (vis SP) and true vertex", 1000, 0,
      1000);

  tr = tfs->make<TTree>("tr", "tr");

  tr->Branch("truth_e_dir", &truth_e_dir);
  tr->Branch("truth_nu_dir", &truth_nu_dir);
  tr->Branch("truth_nu_en", &truth_nu_en);
  tr->Branch("truth_e_en", &truth_e_en);

  tr->Branch("truth_en_deposition", &truth_en_deposition);
  tr->Branch("truth_charge_deposition", &truth_charge_deposition);
  tr->Branch("truth_e_position", &truth_e_position);

  tr->Branch("nu_pdg", &nu_pdg);

  tr->Branch("NParticles", &NParticles);
  tr->Branch("NHits", &NHits);
  tr->Branch("NTrks", &NTrks);
  tr->Branch("charge_U", &charge_U);
  tr->Branch("charge_V", &charge_V);
  tr->Branch("charge_Z", &charge_Z);
  tr->Branch("charge_Z_dist", &charge_Z_dist);
  tr->Branch("charge_corrected", &charge_corrected);
  tr->Branch("cut_charge_loss", &cut_charge_loss);
  tr->Branch("drift_time", &drift_time);

  tr->Branch("reco_e_dir", &reco_e_dir);
  tr->Branch("reco_e_en", &reco_e_en);
  tr->Branch("reco_e_position", &reco_e_position);

  tr->Branch("reco_distance", &reco_distance);
  tr->Branch("primary_trk_id", &primary_trk_id);

  tr->Branch("trk_length", &trk_length);
  tr->Branch("trk_charge", &trk_charge);
  tr->Branch("trk_start_x", &trk_start_x);
  tr->Branch("trk_start_y", &trk_start_y);
  tr->Branch("trk_start_z", &trk_start_z);

  tr->Branch("trk_end_x", &trk_end_x);
  tr->Branch("trk_end_y", &trk_end_y);
  tr->Branch("trk_end_z", &trk_end_z);

  tr->Branch("trk_start_dir_x", &trk_start_dir_x);
  tr->Branch("trk_start_dir_y", &trk_start_dir_y);
  tr->Branch("trk_start_dir_z", &trk_start_dir_z);

  tr->Branch("trk_end_dir_x", &trk_end_dir_x);
  tr->Branch("trk_end_dir_y", &trk_end_dir_y);
  tr->Branch("trk_end_dir_z", &trk_end_dir_z);

} // beginJob()

void dune::PointResTree::endJob() {
  std::cout << "Ending Job" << std::endl;
} // endJob()

void dune::PointResTree::analyze(art::Event const &event) {
  // ana begin
  reset_variables();
  gf = new genFinder();
  // ... Create a map of track IDs to generator labels
  // Get a list of generator names.
  std::vector<art::Handle<std::vector<simb::MCTruth>>> mcHandles;
  event.getManyByType(mcHandles);
  std::vector<std::pair<int, std::string>> track_id_to_label;
  for (auto const &mcHandle : mcHandles) {
    const std::string &sModuleLabel = mcHandle.provenance()->moduleLabel();
    art::FindManyP<simb::MCParticle> findMCParts(mcHandle, event,
                                                 fSimulationLabel);
    std::vector<art::Ptr<simb::MCParticle>> mcParts = findMCParts.at(0);
    for (const art::Ptr<simb::MCParticle> ptr : mcParts) {
      int track_id = ptr->TrackId();
      gf->add(track_id, sModuleLabel);
    }
  }
  if (DEBUG_MSG)
    std::cout << "Getting Truths..." << std::endl;
  // Get MC Truths
  // FIXME: Unify truth gathering MCTruth block.
  writeMCTruths(event);
  // if (fSimulationLabel == "marley") writeMCTruths_marley(event);
  // if (fSimulationLabel == "largeant") writeMCTruths_largeant(event);

  // use spacepoint info to reconstruct an approximate vertex. This vertex is
  // then used to find correct tracks by scanning tracks close to it.
  auto sp_vertex = find_vertex_sp(event);
  auto sp_distance_th = 120; // if track is further away from sp_vertex than
                             // this, it gets thrown out

  if (fCollectTracks) {
    // get track info
    // determine primary track: IF pandora, choose longest track
    if (DEBUG_MSG)
      std::cout << "Getting Tracks..." << std::endl;
    Double_t max_trk_weight = 0;
    int trk_count = 0;
    for (std::string track_label : fTrackModuleLabel) {
      auto trk_list =
          event.getValidHandle<std::vector<recob::Track>>(track_label);
      art::FindManyP<recob::Hit> hitsFromTracks(trk_list, event, track_label);
      trk_count += trk_list->size();
      for (int i = 0; i < (int)trk_list->size();
           i++) { // using index loop to get track idx
        recob::Track const &trk = trk_list->at(i);
        XYZVector trk_start(trk.Start().X(), trk.Start().Y(), trk.Start().Z());
        XYZVector trk_end(trk.End().X(), trk.End().Y(), trk.End().Z());
        // throw away bad tracks
        if (has_sp_vertex && (trk_start - sp_vertex).Mag() > sp_distance_th &&
            (trk_end - sp_vertex).Mag() > sp_distance_th)
          continue;

        // loop through all charges associated with this track to
        // accumulate charge
        auto hits = hitsFromTracks.at(i);
        float charge = 0;
        for (auto hit : hits) {
          charge += hit->Integral();
        }
        trk_charge.push_back(charge);

        Double_t trk_weight =
            charge *
            (track_label == "pandoraShower" ? 10.0 : 1.0); // prioritize showers
        if (trk_weight > max_trk_weight) {
          max_trk_weight = trk_weight;
          primary_trk_id = trk_length.size(); // most recent entry
          primary_trk_loc = track_loc(track_label, i);
        }
        // write track info

        // track length debugger
        if (trk.Length() > 100) {
          double length = 0;
          size_t idx = 0;

          auto iLast = trk.Trajectory().LastValidPoint();
          while (idx != iLast) {
            auto current_loc = trk.Trajectory().LocationAtPoint(idx);
            auto next_idx = trk.Trajectory().NextValidPoint(idx);
            auto next_loc = trk.Trajectory().LocationAtPoint(next_idx);
            length += (next_loc - current_loc).R();
            idx = next_idx;
          } // while
          std::cout << "Trk.Length(): " << trk.Length()
                    << "\tSum of distances: " << length << std::endl;
        }
        trk_length.push_back(trk.Length());
        trk_start_x.push_back(trk.Start().X());
        trk_start_y.push_back(trk.Start().Y());
        trk_start_z.push_back(trk.Start().Z());
        trk_end_x.push_back(trk.End().X());
        trk_end_y.push_back(trk.End().Y());
        trk_end_z.push_back(trk.End().Z());

        trk_start_dir_x.push_back(trk.StartDirection().X());
        trk_start_dir_y.push_back(trk.StartDirection().Y());
        trk_start_dir_z.push_back(trk.StartDirection().Z());
        trk_end_dir_x.push_back(trk.EndDirection().X());
        trk_end_dir_y.push_back(trk.EndDirection().Y());
        trk_end_dir_z.push_back(trk.EndDirection().Z());

      } // end loop through trks
    }   // loop through track labels
    NTrks = trk_length.size();
    if (DEBUG_MSG) {
      std::cout << "Number of Tracks cut: " << trk_count - NTrks << std::endl;
      std::cout << "NTrks: " << NTrks << std::endl;
      std::cout << "SP Vertex Loc: " << sp_vertex.X() << " " << sp_vertex.Y()
                << " " << sp_vertex.Z() << std::endl;
      if (NTrks == 0)
        std::cout << "No Tracks FOUND!!" << std::endl;
      else {
        std::cout << "Primary Track Start: " << trk_start_x.at(primary_trk_id)
                  << " " << trk_start_y.at(primary_trk_id) << " "
                  << trk_start_z.at(primary_trk_id) << std::endl;
        std::cout << "Primary Track End: " << trk_end_x[primary_trk_id] << " "
                  << trk_end_y[primary_trk_id] << " "
                  << trk_end_z[primary_trk_id] << std::endl;
      }
    }

    // std::cout <<
    // std::distance(trk_charge.begin(),std::max_element(trk_charge.begin(),
    // trk_charge.end()))<<std::endl; std::cout << primary_trk_id <<
    // std::endl;
    calculate_electron_direction();
  }
  // auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService
  // const>()->DataFor(event);
  // get hit info
  if (DEBUG_MSG)
    std::cout << "Getting Hits..." << std::endl;
  auto hit_handle =
      event.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
  std::vector<art::Ptr<recob::Hit>> hit_list;
  art::fill_ptr_vector(hit_list, hit_handle);
  auto const detectorPropertiesData =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(
          event);
  art::ServiceHandle<cheat::BackTrackerService> bt;
  auto hit_to_sp = art::FindManyP<recob::SpacePoint>(hit_handle, event,
                                                     fHitToSpacePointLabel);
  float NHits_Z = 0;
  NHits = hit_list.size();
  double uncut_charge = 0;
  // double max_amplitude = 0.;
  // int brightest_hit_idx = -1;
  for (int i = 0; i < NHits; i++) {
    auto const hit = hit_list.at(i);
    if (hit->View() == geo::kU)
      charge_U += hit->Integral();
    if (hit->View() == geo::kV)
      charge_V += hit->Integral();
    if (hit->View() == geo::kZ)
      uncut_charge += hit->Integral();
    // for collection plane, apply distance cut
    if (hit->View() == geo::kZ) {
      auto spacepoints = hit_to_sp.at(i);
      NHits_Z++;
      write_multi_distance(event, *hit, spacepoints);
      if (distance_cut(event, *hit, spacepoints)) {
        continue;
      } else
        charge_Z += hit->Integral();
      charge_Z += hit->Integral();
    }
  } // end loop through hits

  if (DEBUG_MSG)
    std::cout << "Done looping through hits" << std::endl;
  cut_charge_loss = charge_Z / uncut_charge;
  // std::cout<<"Spacepoint reco hits / total hits: " << NspHits_Z/NHits_Z <<
  // std::endl; std::cout<<NHits_Z/NHits<< std::endl;
  // hit_distance->Fill(NspHits_Z/NHits_Z);
  if (NTrks != 0)
    reco_distance = (truth_e_position - reco_e_position).Mag();

  if (DEBUG_MSG)
    std::cout << "Calculate Electron Energy" << std::endl;
  calculate_electron_energy(event);

  delete gf;
  tr->Fill();
} // analyze

//=================================================================================================
// Private Helper methods

void dune::PointResTree::reset_variables() {
  truth_nu_dir.SetXYZ(-2, -2, -2);
  truth_e_dir.SetXYZ(-2, -2, -2);
  reco_e_dir.SetXYZ(-2, -2, -2);
  reco_e_position.SetXYZ(-2, -2, -2);

  reco_e_en = truth_nu_en = truth_e_en = -10;
  nu_pdg = 0;
  has_sp_vertex = false;
  primary_trk_loc = track_loc("NA", -1);

  NParticles = NHits = NTrks = 0;
  charge_U = charge_V = charge_Z = charge_corrected = 0;
  charge_Z_dist.resize(20, 0);
  std::fill(charge_Z_dist.begin(), charge_Z_dist.end(), 0);
  cut_charge_loss = 0;
  drift_time = -10;
  reco_distance = -10;
  primary_trk_id = -1;
  truth_en_deposition = 0;
  truth_charge_deposition = 0;

  trk_length.clear();
  trk_charge.clear();
  trk_start_x.clear();
  trk_start_y.clear();
  trk_start_z.clear();
  trk_end_x.clear();
  trk_end_y.clear();
  trk_end_z.clear();
  trk_start_dir_x.clear();
  trk_start_dir_y.clear();
  trk_start_dir_z.clear();
  trk_end_dir_x.clear();
  trk_end_dir_y.clear();
  trk_end_dir_z.clear();

} // reset_variables

void dune::PointResTree::writeMCTruths(art::Event const &event) {
  auto const &truth_blocks =
      event.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
  if (truth_blocks->size() != 1)
    std::cout << "WARNING: " << truth_blocks->size()
              << " MCTruth blocks found, expect 1" << std::endl;
  Bool_t wrote_nu = kFALSE;
  Bool_t wrote_e = kFALSE;
  for (auto truth_block : *truth_blocks) {
    for (int i = 0; i < truth_block.NParticles(); i++) {
      simb::MCParticle part = truth_block.GetParticle(i);
      if (abs(part.PdgCode()) == 12 || abs(part.PdgCode()) == 14 ||
          abs(part.PdgCode()) == 16) { // neutrino
        if (wrote_nu) {
          std::cout << "WARNING: More than one neutrino is found" << std::endl;
        }
        nu_pdg = part.PdgCode();
        truth_nu_dir.SetXYZ(part.Px(), part.Py(), part.Pz());
        truth_nu_dir.SetMag(1.0);
        truth_nu_en = part.E() * 1000;
        wrote_nu = kTRUE;
        continue;
      }

      if (part.PdgCode() == 11) { // electron
        if (wrote_e) {
          std::cout << "WARNING: More than one electron is found" << std::endl;
        }
        truth_e_dir.SetXYZ(part.Px(), part.Py(), part.Pz());
        truth_e_dir.SetMag(1);
        truth_e_en = (part.E() - part.Mass()) * 1000;
        truth_e_position = part.Position().Vect();
        wrote_e = kTRUE;
        continue;
      }
    }
  } // end go through truth
  if (!wrote_nu)
    std::cout << "WARNING: No neutrino was found!" << std::endl;
  if (!wrote_e)
    std::cout << "WARNING: No electron was found!" << std::endl;

  auto const &particle_list =
      event.getValidHandle<std::vector<simb::MCParticle>>(fSimulationLabel);
  NParticles = particle_list->size();
}

void dune::PointResTree::calculate_electron_energy(art::Event const &event) {
  // use collection plane charge for reconstruction
  Double_t charge_raw = charge_Z;
  // drift correction
  if (NTrks == 0 ||
      primary_trk_id == -1) { // no tracks found, do not do drift correction
    charge_corrected = charge_raw;
    return;
  }

  auto const detectorClockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(
          event);
  auto opflashHandle =
      event.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashLabel);
  std::vector<double> totalPEs;
  std::vector<float> flashtimes;
  for (unsigned int i = 0; i < opflashHandle->size(); i++) {
    auto flash = opflashHandle->at(i);
    totalPEs.push_back(flash.TotalPE());
    flashtimes.push_back(flash.Time()); // units us
  }

  std::vector<double> hitPeakTimes;
  drift_time = 0.0;
  // Create FindManyP object to find hits associated with Tracks
  auto primary_trk_label = primary_trk_loc.first;
  auto primary_trk_idx = primary_trk_loc.second;
  auto tracklistHandle =
      event.getValidHandle<std::vector<recob::Track>>(primary_trk_label);
  art::FindManyP<recob::Hit> hitsFromTracks(tracklistHandle, event,
                                            primary_trk_label);
  // hits in primary track:
  if (DEBUG_MSG)
    std::cout << "primary_trk_label is " << primary_trk_label
              << " hitsFromTracks size: " << hitsFromTracks.size() << std::endl;
  auto hits_primary = hitsFromTracks.at(primary_trk_idx); // pointers to hits
  // for(auto hit : hits_primary){
  //     std::cout<<detectorClockData.TPCTick2Time(hit->PeakTime()) <<
  //     std::endl;
  // }
  if (DEBUG_MSG)
    std::cout << "Done indexing hitsFromTracks" << std::endl;
  Double_t hitTime =
      detectorClockData.TPCTick2Time(hits_primary[0]->PeakTime());
  Double_t current_drift_time = 0;
  while (totalPEs.size() != 0) {
    int max_flash_idx = std::distance(
        totalPEs.begin(), std::max_element(totalPEs.begin(), totalPEs.end()));
    current_drift_time = hitTime - flashtimes[max_flash_idx];
    // std::cout<<current_drift_time<<std::endl;
    if (current_drift_time > 0 &&
        current_drift_time < 2400.0) { // drift time found
      drift_time = current_drift_time;
      break;
    }
    // if the flashtime found is not physical
    totalPEs.erase(totalPEs.begin() + max_flash_idx);
    flashtimes.erase(flashtimes.begin() + max_flash_idx);
  }
  // Double_t APA_distance = abs(truth_e_position.X());
  auto const detectorPropertiesData =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(
          event);
  Double_t electron_lifetime =
      detectorPropertiesData.ElectronLifetime(); // us, probably 10400
  // Double_t drift_velocity = detectorPropertiesData.DriftVelocity();
  // drift_time = APA_distance/drift_velocity; //using truth!!
  charge_corrected =
      charge_raw * ROOT::Math::exp(drift_time / electron_lifetime);

  // TODO: Calculate reco energy. This needs to be done after the
  //      charge-energy linear fit parameters are finalized.
}

void dune::PointResTree::calculate_electron_direction() {
  // std::cout << "Daughter Flipping" << std::endl;
  // std::cout << "NTrk: " << NTrks << std::endl;
  // std::cout << "Vector size: " << trk_start_x.size() << std::endl;
  // std::cout << "Primary Index:" << primary_trk_id << std::endl;
  if (primary_trk_id < 0)
    return; // if no trk is in reco, do nothing
  // double dist_thresh = 400;
  // using XYZVector=TVector3;
  // reset primary position after index update
  XYZVector primary_start(trk_start_x.at(primary_trk_id),
                          trk_start_y.at(primary_trk_id),
                          trk_start_z.at(primary_trk_id));

  XYZVector primary_end(trk_end_x.at(primary_trk_id),
                        trk_end_y.at(primary_trk_id),
                        trk_end_z.at(primary_trk_id));
  XYZVector primary_forward(trk_start_dir_x.at(primary_trk_id),
                            trk_start_dir_y.at(primary_trk_id),
                            trk_start_dir_z.at(primary_trk_id));
  XYZVector primary_backward(-trk_end_dir_x.at(primary_trk_id),
                             -trk_end_dir_y.at(primary_trk_id),
                             -trk_end_dir_z.at(primary_trk_id));

  Double_t sum_cos_forward = 0.0;
  Double_t sum_cos_backward = 0.0;
  XYZVector current_trk;
  XYZVector difference_start, difference_end;
  // loop through all tracks to find avg cos between primary track and all
  // daughter tracks.
  for (int i = 0; i < NTrks; i++) {
    if (i == primary_trk_id)
      continue; // don't include primary track
    current_trk.SetXYZ(trk_start_x.at(i), trk_start_y.at(i), trk_start_z.at(i));
    difference_start = current_trk - primary_start;
    difference_end = current_trk - primary_end;
    sum_cos_forward +=
        ROOT::Math::VectorUtil::CosTheta(difference_start, primary_forward);
    sum_cos_backward +=
        ROOT::Math::VectorUtil::CosTheta(difference_end, primary_backward);
  }
  // std::cout<<sum_cos_forward << std::endl;
  // std::cout<<sum_cos_backward<<std::endl;
  if (sum_cos_backward > sum_cos_forward) {
    reco_e_dir = primary_backward;
    reco_e_position = primary_end;
    // std::cout<<"Flipped"<<std::endl;
  } else {
    reco_e_dir = primary_forward;
    reco_e_position = primary_start;
  }
}

Bool_t dune::PointResTree::distance_cut(
    art::Event const &event, recob::Hit const &hit,
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints) {
  if (primary_trk_id < 0)
    return kFALSE;           // don't cut anything if there are no tracks
  double dist_th = 14.0 * 5; // 14 is radiation length
  Double_t distance = 0;
  // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  // auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService
  // const>()->DataFor(event); art::ServiceHandle<cheat::BackTrackerService>
  // bt;

  if (spacepoints.size()) {
    // average of all spacepoint recos
    auto spacepoint = spacepoints.at(0);
    auto pos = XYZVector(spacepoint->XYZ());
    distance += (pos - reco_e_position).Mag();
  } else { // no spacepoint found, do 2D reconstruction (X, Z)
    geo::GeometryCore const &geom = *(lar::providerFrom<geo::Geometry>());
    auto const detProperties =
        art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(
            event);
    auto wireID = hit.WireID();
    auto plane = geom.Plane(wireID);
    auto pos = plane.Wire(wireID).GetCenter<geo::Point_t>();
    plane.DriftPoint(pos,
                     -abs(detProperties.ConvertTicksToX(
                         hit.PeakTime(), wireID))); // drift away from plane
    distance = std::pow(pos.X() - reco_e_position.X(), 2) +
               std::pow(pos.Z() - reco_e_position.Z(), 2);
    distance = sqrt(distance);
  }

  return (distance > dist_th);
}

void dune::PointResTree::write_multi_distance(
    art::Event const &event, recob::Hit const &hit,
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints) {
  if (primary_trk_id < 0)
    return; // don't cut anything if there are no tracks
  Double_t distance = 0;
  // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  // auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService
  // const>()->DataFor(event); art::ServiceHandle<cheat::BackTrackerService>
  // bt;

  if (spacepoints.size()) {
    // average of all spacepoint recos
    auto spacepoint = spacepoints.at(0);
    auto pos = XYZVector(spacepoint->XYZ());
    distance += (pos - reco_e_position).Mag();
  } else { // no spacepoint found, do 2D reconstruction (X, Z)
    geo::GeometryCore const &geom = *(lar::providerFrom<geo::Geometry>());
    auto const detProperties =
        art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(
            event);
    auto wireID = hit.WireID();
    auto plane = geom.Plane(wireID);
    auto pos = plane.Wire(wireID).GetCenter<geo::Point_t>();
    plane.DriftPoint(pos,
                     -abs(detProperties.ConvertTicksToX(
                         hit.PeakTime(), wireID))); // drift away from plane
    distance = std::pow(pos.X() - truth_e_position.X(), 2) +
               std::pow(pos.Z() - truth_e_position.Z(), 2);
    distance = sqrt(distance);
  }

  // write charge to respective distance cut bins
  for (int i = 0; i < 20; i++) {
    double cut = 14.0 * (i + 1); // each bin signifies a cut applied at 14*bin
    if (distance < cut)
      charge_Z_dist.at(i) += hit.Integral();
  }
}

XYZVector dune::PointResTree::find_vertex_sp(art::Event const &event) {
  auto sp_handle = event.getValidHandle<std::vector<recob::SpacePoint>>(
      fSpacePointModuleLabel);
  auto sp_to_hit =
      art::FindManyP<recob::Hit>(sp_handle, event, fHitToSpacePointLabel);
  int N_spacepoint = sp_handle->size();
  // Using a priority queue to get the top k nmber of spacepoints. Runtime is
  // O(NlogK)
  using AVGAMP_SP = std::pair<double, recob::SpacePoint>;
  std::priority_queue<AVGAMP_SP, std::vector<AVGAMP_SP>,
                      std::greater<AVGAMP_SP>>
      sp_brightness_queue;
  size_t queue_size = 10; // number of spacepoints to keep
  for (int i = 0; i < N_spacepoint; i++) {
    auto hits = sp_to_hit.at(i);
    double avg_amplitude = 0;
    for (auto hit : hits) {
      // avg_amplitude += hit->PeakAmplitude();
      if (hit->View() == geo::kZ)
        avg_amplitude = hit->PeakAmplitude();
    }
    // avg_amplitude /= hits.size();
    if (sp_brightness_queue.size() < queue_size) {
      sp_brightness_queue.push(AVGAMP_SP(avg_amplitude, sp_handle->at(i)));
    } else if (sp_brightness_queue.top().first < avg_amplitude) {
      sp_brightness_queue.pop();
      sp_brightness_queue.push(AVGAMP_SP(avg_amplitude, sp_handle->at(i)));
    }
  }
  std::vector<recob::SpacePoint> brightest_sp;
  while (!sp_brightness_queue.empty()) {
    auto current = sp_brightness_queue.top().second;
    brightest_sp.push_back(current);
    sp_brightness_queue.pop();
  }
  // selects a spacepoint that is close to most of the other spacepoints.
  // In the case of a spacepoint caused by something other than the SN stub,
  // it should be far away from the majority of the other hits.
  double min_nearest = 1e9;
  recob::SpacePoint vertex;
  for (auto sp : brightest_sp) {
    double nearest = 1e9;
    // double nearest = 0;
    XYZVector p1(sp.XYZ());
    for (auto other_sp : brightest_sp) {
      XYZVector p2(other_sp.XYZ());
      if (p1 == p2)
        continue;
      auto dist = (p1 - p2).Mag();
      nearest = std::min(dist, nearest);
      nearest += dist;
    }
    if (nearest < min_nearest) {
      min_nearest = nearest;
      vertex = sp;
      has_sp_vertex = true;
    }
  }
  XYZVector res(vertex.XYZ());
  if (has_sp_vertex)
    vertex_reco_distance->Fill((res - truth_e_position).Mag());
  return res;
}
