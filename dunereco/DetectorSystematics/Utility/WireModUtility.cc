#include "WireModUtility.hh"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"


//For debugging purposes
//--- MakeTrackIDtoPDGMap ---
void sys::WireModUtility::MakeTrackIDtoPDGMap(
        const std::vector<simb::MCParticle>& mcParticles)
{
    trackID_to_PDG.clear();
    for (const auto& p : mcParticles) {
        trackID_to_PDG[p.TrackId()] = p.PdgCode();
    }
}

//--- CalcROIProperties ---
sys::WireModUtility::ROIProperties_t sys::WireModUtility::CalcROIProperties(recob::Wire const& wire, size_t const& roi_idx)
{
  // get the ROI. Need to verify ROI tick conventions and signal calibration units.
  recob::Wire::RegionsOfInterest_t::datarange_t const& roi = wire.SignalROI().get_ranges()[roi_idx];

  // initialize the return value
  ROIProperties_t roi_vals;
  roi_vals.channel = wire.Channel(); // DUNE_MOD: confirm channel mapping consistency
  roi_vals.view    = wire.View(); // DUNE_MOD: DUNE may have different plane/view conventions
  roi_vals.begin   = roi.begin_index();
  roi_vals.end     = roi.end_index();
  roi_vals.center  = 0;
  roi_vals.total_q = 0;
  roi_vals.sigma   = 0;

  // loop over the roi and find the charge-weighted center and total charge
  auto const& roi_data = roi.data();
  for (size_t i_t = 0; i_t < roi_data.size(); ++i_t)
  {
    roi_vals.center += roi_data[i_t]*(i_t+roi_vals.begin);
    roi_vals.total_q += roi_data[i_t];
  }

  //Do we need to protect against division by 0?
  roi_vals.center = roi_vals.center/roi_vals.total_q;

  // get the width (again charge-weighted)
  //   // if the ROI is only one tick set the cent to the middle of the tick and width to 0.5
  for (size_t i_t = 0; i_t<roi_data.size(); ++i_t)
  {
    roi_vals.sigma += roi_data[i_t]*(i_t+roi_vals.begin-roi_vals.center)*(i_t+roi_vals.begin-roi_vals.center);
  }
  //Do we need to protect against division by 0?
  roi_vals.sigma = std::sqrt(roi_vals.sigma/roi_vals.total_q);
  if (roi_vals.end-roi_vals.begin == 1)
  {
    roi_vals.center += 0.5;
    roi_vals.sigma   = 0.5;
  }
  
  // return the calc'd properies
  return roi_vals;
}


//--- GetTargetROIs ---
std::vector<std::pair<unsigned int, unsigned int>>
sys::WireModUtility::GetTargetROIs(sim::SimEnergyDeposit const& shifted_edep, double offset)
{
// DUNE_MOD: This function is highly geometry-dependent.
// // Validate PositionToTPCptr behavior with DUNE multi-TPC layouts.
  std::vector<std::pair<unsigned int, unsigned int>> target_roi_vec;

  geo::TPCGeo const* curTPCGeomPtr = geometry->PositionToTPCptr(shifted_edep.MidPoint());
  if (curTPCGeomPtr == nullptr)
  {
    return target_roi_vec;
  }
// DUNE_MOD: Iterate<PlaneGeo> assumes plane ordering compatible with SBN.
// // Verify plane indexing and orientation in DUNE.
  for (auto const& plane : wireReadout->Iterate<geo::PlaneGeo>(curTPCGeomPtr->ID())) {

    int wireNumber = int(0.5 + plane.WireCoordinate(shifted_edep.MidPoint()));
    if ((wireNumber < 0) || (wireNumber >= (int)plane.Nwires()))
    {
      continue;
    }

    geo::WireID edep_wireID = plane.NearestWireID(shifted_edep.MidPoint());

    // DUNE_MOD: planeXInWindow and planeXToTick rely on detector timing model.
    // // Confirm tickOffset and drift/timing calibration in DUNE.
    if (planeXInWindow(shifted_edep.X(), plane, *curTPCGeomPtr, offset + tickOffset))
    {
      target_roi_vec.emplace_back(
        wireReadout->PlaneWireToChannel(edep_wireID),
        std::round(planeXToTick(shifted_edep.X(), plane, *curTPCGeomPtr, offset + tickOffset))
    );
    }
  }

  return target_roi_vec;
}

//--- GetHitTargetROIs ---
std::vector<std::pair<unsigned int, unsigned int>>
sys::WireModUtility::GetHitTargetROIs(recob::Hit const& hit)
{
  // DUNE_MOD: Verify hit timing reference frame matches wire tick convention.
  std::vector<std::pair<unsigned int, unsigned int>> target_roi_vec;

  int hit_wire = hit.Channel();
  int hit_tick = int(round(hit.PeakTime()));

  // DUNE_MOD: readoutWindowTicks and tickOffset must match DUNE readout window.
  if (hit_tick < tickOffset || hit_tick >= readoutWindowTicks + tickOffset)
    return target_roi_vec;

  target_roi_vec.emplace_back((unsigned int) hit_wire, (unsigned int) hit_tick);
  return target_roi_vec;
}

//--- FillROIMatchedEdepMap ---
void sys::WireModUtility::FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const& edepVec, std::vector<recob::Wire> const& wireVec, double offset)
{
  // DUNE_MOD: Assumes channel numbering consistent between geometry and wireVec.
  ROIMatchedEdepMap.clear();

  

  std::unordered_map<unsigned int,unsigned int> wireChannelMap;
  for (size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    wireChannelMap[wireVec[i_w].Channel()] = i_w;

  event_counter++;
  float xmin=1e6, xmax=-1e6, zmin=1e6, zmax=1e-6;


  TCanvas* c1 = new TCanvas("c1","Edep XZ Display",800,600);
  TH2F* hMatched = new TH2F("hMatched","Matched Edeps X vs Z;X [cm];Z [cm]", 1500, -350, 350, 3000, 0, 1500);  
  TH2F* hUnMatched = new TH2F("hUnMatched","Non matched Edeps X vs Z;X [cm];Z [cm]", 1500, -350, 350, 3000, 0, 1500);
  TH1F* hMatchedE = new TH1F("hMatchedE","Matched Edeps energy; Energy; Counts", 100, 0, 2);
  TH1F* hUnMatchedE = new TH1F("hUnMatchedE","Non matched Edeps energy; Energy; Counts", 100, 0, 2);
  TH1F* hMatchedPDG = new TH1F("hMatchedPDG", "Matched edeps PDG code; PDG code; Counts",12, -0.5, 11.5);
  TH1F* hUnMatchedPDG = new TH1F("hUnMatchedPDG", "Non matched edeps PDG code; PDG code; Counts",12, -0.5, 11.5);

  hMatchedPDG->GetXaxis()->SetBinLabel(1, "PDG=0");
  hMatchedPDG->GetXaxis()->SetBinLabel(2, "e^{#pm}");
  hMatchedPDG->GetXaxis()->SetBinLabel(3, "#mu^{#pm}");
  hMatchedPDG->GetXaxis()->SetBinLabel(10, "#gamma");
  hMatchedPDG->GetXaxis()->SetBinLabel(4, "p");
  hMatchedPDG->GetXaxis()->SetBinLabel(5, "n");
  hMatchedPDG->GetXaxis()->SetBinLabel(6, "#pi^{#pm}");
  hMatchedPDG->GetXaxis()->SetBinLabel(7, "#pi^{0}");
  hMatchedPDG->GetXaxis()->SetBinLabel(8, "K^{#pm}");
  hMatchedPDG->GetXaxis()->SetBinLabel(9, "K^{0}");
  hMatchedPDG->GetXaxis()->SetBinLabel(11, "Nuclei");
  hMatchedPDG->GetXaxis()->SetBinLabel(12,"Other");

  hUnMatchedPDG->GetXaxis()->SetBinLabel(1, "PDG=0");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(2, "e^{#pm}");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(3, "#mu^{#pm}");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(10, "#gamma");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(4, "p");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(5, "n");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(6, "#pi^{#pm}");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(7, "#pi^{0}");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(8, "K^{#pm}");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(9, "K^{0}");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(11, "Nuclei");
  hUnMatchedPDG->GetXaxis()->SetBinLabel(12,"Other");

  hMatched->SetMarkerStyle(2);
  hMatched->SetMarkerColor(kBlue);
  hMatched->SetStats(0);
  hMatchedE->SetLineColor(kBlue);
  hMatchedE->SetStats(0);
  hMatchedPDG->SetLineColor(kBlue);
  hMatchedPDG->SetStats(0);

  hUnMatched->SetMarkerStyle(5);
  hUnMatched->SetMarkerColor(kRed);
  hUnMatched->SetStats(0);
  hUnMatchedE->SetLineColor(kRed);
  hUnMatchedE->SetStats(0);
  hUnMatchedPDG->SetLineColor(kRed);
  hUnMatchedPDG->SetStats(0);

  for (size_t i_e = 0; i_e < edepVec.size(); ++i_e)
  {
    bool isMatched=false;
    auto const& edep = edepVec[i_e];
    
    int trackID=std::abs(edep.TrackID());
    int pdg=0, pdg_rebin=11;
    auto it = trackID_to_PDG.find(trackID);
        if (it != trackID_to_PDG.end()) pdg = it->second; 
    if (pdg==-11 || pdg==11) pdg_rebin=1; //electron, positron
    if (pdg==-13 || pdg==13) pdg_rebin=2; //muons
    if (pdg==-2212 || pdg==2212) pdg_rebin=3; //protons
    if (pdg==-2112 || pdg==2112) pdg_rebin=4; //neutrons
    if (pdg==-211 || pdg==211) pdg_rebin=5; //charged pions
    if (pdg==-111 || pdg==111) pdg_rebin=6; //neutral pions
    if (pdg==-321 || pdg==321) pdg_rebin=7; //charged kaons
    if (pdg==310 || pdg==130) pdg_rebin=8; //neutral kaons
    if (pdg==22) pdg_rebin=9; //gamma
    if (pdg>1e9) pdg_rebin=10; //nuclei
    if (pdg==0) pdg_rebin=0;
    auto target_rois = GetTargetROIs(edepVec[i_e], offset);
    for (auto const& target_roi : target_rois)
    {
      
      if (wireChannelMap.find(target_roi.first) == wireChannelMap.end())
      {  
        continue;
      }
      auto const& target_wire = wireVec.at(wireChannelMap[target_roi.first]);

      if (not target_wire.SignalROI().is_valid())
      {
        //std::cout<<"invalid Signal ROI for this wire"<<std::endl;
        continue;
      }
      if (target_wire.SignalROI().empty())
      {
        //std::cout<<"empty signal ROI for this wire"<<std::endl;
        continue;
      }
      if (target_wire.SignalROI().n_ranges() == 0)
      {
        /*std::cout<<"0 ranges for this wire's signal ROI"<<std::endl;
        for (int delta = -2; delta <= 2; ++delta) {
          int neigh_channel = target_roi.first + delta;
          auto const& neigh_wire = wireVec.at(wireChannelMap[neigh_channel]);  // however you access wires
          if (neigh_wire.SignalROI().n_ranges() >0)
          {
            std::cout<<"Neighbour channel at "<<delta<<" has signal ROIs"<<std::endl;
          }
        }*/
       
        continue;
      }


      if (target_wire.SignalROI().size() <= target_roi.second)
      {
        //std::cout<<"signa ROI size smaller than target_roi tick"<<std::endl;
        continue;
      }
      if (target_wire.SignalROI().is_void(target_roi.second))
      {
        continue;
      }

      auto range_number = target_wire.SignalROI().find_range_iterator(target_roi.second) - target_wire.SignalROI().begin_range();
    
      ROIMatchedEdepMap[std::make_pair(target_wire.Channel(),range_number)].push_back(i_e);
      isMatched=true;
    }
    if (edep.X()>xmax) xmax=edep.X();
    else if(edep.X()<xmin) xmin=edep.X();
    if (edep.Z()>zmax) zmax=edep.Z();
    else if(edep.Z()<zmin) zmin=edep.Z();    
    if (isMatched){
      hMatched->Fill(edep.X(), edep.Z());
      hMatchedE->Fill(edep.Energy());
      hMatchedPDG->Fill(pdg_rebin);
      //std::cout<<"matched edep PDG: "<<pdg<<std::endl;
    }
    else{
      hUnMatched->Fill(edep.X(), edep.Z()); 
      hUnMatchedE->Fill(edep.Energy());
      hUnMatchedPDG->Fill(pdg_rebin);
      //std::cout<<"unmatched edep PDG: "<<pdg<<std::endl;
    }
  }

  std::cout<<"xmin = "<<xmin<<", xmax = "<<xmax<<", zmin = "<<zmin<<", zmax = "<<zmax<<std::endl;

  hMatched->GetXaxis()->SetRangeUser(xmin-5, xmax-5);
  hMatched->GetYaxis()->SetRangeUser(zmin-5, zmax-5);
  hUnMatched->GetXaxis()->SetRangeUser(xmin-5, xmax-5);
  hUnMatched->GetYaxis()->SetRangeUser(zmin-5, zmax-5);

  std::ostringstream filename;
  filename << "Edep_XZ_event_" << event_counter << ".pdf";
  std::string fname = filename.str(); 

  hMatched->Draw("P");
  hUnMatched->Draw("Psame");
  c1->SaveAs((fname+"(").c_str());
  
  hMatchedE->Draw();
  hUnMatchedE->Draw("same");
  c1->SaveAs(fname.c_str());

  c1->SetLogy();
  hMatchedPDG->Draw();
  hUnMatchedPDG->Draw("same");
  c1->SaveAs((fname+")").c_str());
  delete c1;
}

// DUNE_MOD: Remaining functions follow the same pattern.
// // Key areas requiring DUNE validation:
// // * Tick conversion using DetectorPropertiesData
// // * Assumption of exactly 3 planes (arrays sized [3])
// // * Plane indexing via plane.ID().Plane
// // * Geometry lookups using PositionToTPC
// // * Spline/graph arrays indexed by plane
// // * Timing offsets (tickOffset, readoutWindowTicks)
// // * Channel/view mapping consistency
//
// // In particular, search for:
// // for(size_t i_p = 0; i_p < 3; ++i_p)
// // These assume a 3-plane detector and may need to be generalized
// // or validated for DUNE geometry.


//--- FillROIMatchedHitMap ---
void sys::WireModUtility::FillROIMatchedHitMap(std::vector<recob::Hit> const& hitVec, std::vector<recob::Wire> const& wireVec)
{
  // clear the map in case it was already set
  ROIMatchedHitMap.clear();

  // get the channel from each wire and set up a map between them
  std::unordered_map<unsigned int,unsigned int> wireChannelMap;
  for (size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    wireChannelMap[wireVec[i_w].Channel()] = i_w;

  // loop over hits
  for (size_t i_h = 0; i_h < hitVec.size(); ++i_h)
  {
    // get the ROIs
    //     // <channel number, tick time>
    std::vector<std::pair<unsigned int, unsigned int>> target_rois = GetHitTargetROIs(hitVec[i_h]);

    // loop over ROI and match the energy deposits with wires
    for (auto const& target_roi : target_rois)
    {
      // if we can't find the wire, skip it
      if (wireChannelMap.find(target_roi.first) == wireChannelMap.end())
        continue;

      auto const& target_wire = wireVec.at(wireChannelMap[target_roi.first]);

      // if there are no ticks in in the wire skip it
      //       // likewise if there's nothing in the region of interst
      if (not target_wire.SignalROI().is_valid()              ||
          target_wire.SignalROI().empty()                     ||
          target_wire.SignalROI().n_ranges() == 0             ||
          target_wire.SignalROI().size() <= target_roi.second ||
          target_wire.SignalROI().is_void(target_roi.second)   )
        continue;

      // which range is it?
      auto range_number = target_wire.SignalROI().find_range_iterator(target_roi.second) - target_wire.SignalROI().begin_range();

      // pupluate the map
      ROIMatchedHitMap[std::make_pair(target_wire.Channel(),range_number)].push_back(i_h);
    }
  }
}

//--- CalcSubROIProperties ---
std::vector<sys::WireModUtility::SubROIProperties_t> sys::WireModUtility::CalcSubROIProperties(sys::WireModUtility::ROIProperties_t const& roi_properties, std::vector<const recob::Hit*> const& hitPtrVec)
{
  std::vector<sys::WireModUtility::SubROIProperties_t> subroi_properties_vec;
  sys::WireModUtility::SubROIProperties_t subroi_properties;
  subroi_properties.channel = roi_properties.channel;
  subroi_properties.view    = roi_properties.view;

  // if this ROI doesn't contain any hits, define subROI based on ROI properities
  //   // otherwise, define subROIs based on hits
  if (hitPtrVec.size() == 0)
  {
    subroi_properties.key     = std::make_pair(roi_properties.key, 0);
    subroi_properties.total_q = roi_properties.total_q;
    subroi_properties.center  = roi_properties.center;
    subroi_properties.sigma   = roi_properties.sigma;
    subroi_properties_vec.push_back(subroi_properties);
  } else
  {
    for (unsigned int i_h=0; i_h < hitPtrVec.size(); ++i_h)
    {
      auto hit_ptr = hitPtrVec[i_h];
      subroi_properties.key     = std::make_pair(roi_properties.key, i_h);
      subroi_properties.total_q = hit_ptr->Integral();
      subroi_properties.center  = hit_ptr->PeakTime();
      subroi_properties.sigma   = hit_ptr->RMS();
      subroi_properties_vec.push_back(subroi_properties);
    }
  }

  return subroi_properties_vec;
}

//--- MatchEdepsToSubROIs ---
std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> sys::WireModUtility::MatchEdepsToSubROIs(std::vector<sys::WireModUtility::SubROIProperties_t> const& subROIPropVec, 
                                                                                                                  std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec, double offset)
{
  // for each TrackID, which EDeps are associated with it? keys are TrackIDs
  std::map<int, std::vector<const sim::SimEnergyDeposit*>> TrackIDMatchedEDepMap;

  // total energy of EDeps matched to the ROI (not strictly necessary, but useful for understanding/development
  double total_energy = 0.0;
  
  // loop over edeps, fill TrackIDMatchedEDepMap and calculate total energy
  for (auto const& edep_ptr : edepPtrVec)
  {
    TrackIDMatchedEDepMap[edep_ptr->TrackID()].push_back(edep_ptr);
    total_energy += edep_ptr->E();
  }

  // calculate EDep properties by TrackID
  std::map<int, sys::WireModUtility::TruthProperties_t> TrackIDMatchedPropertyMap;
  for (auto const& track_edeps : TrackIDMatchedEDepMap)
    TrackIDMatchedPropertyMap[track_edeps.first] = CalcPropertiesFromEdeps(track_edeps.second, offset);

  // map it all out
  std::map<unsigned int, std::vector<unsigned int>> EDepMatchedSubROIMap;        // keys are indexes of edepPtrVec, values are vectors of indexes of subROIPropVec
  std::map<int, std::unordered_set<unsigned int>>   TrackIDMatchedSubROIMap;     // keys are TrackIDs, values are sets of indexes of subROIPropVec
  std::map<unsigned int, std::vector<unsigned int>> SubROIMatchedEDepMap;        // keys are indexes of subROIPropVec, values are vectors of indexes of edepPtrVec
  std::map<unsigned int, std::map<int, double>>     SubROIMatchedTrackEnergyMap; // keys are indexes of subROIPropVec, values are maps of TrackIDs to matched energy (in MeV)

  // loop over EDeps
  for (unsigned int i_e = 0; i_e < edepPtrVec.size(); ++i_e)
  {
    // get EDep properties
    auto edep_ptr  = edepPtrVec[i_e];
    const geo::TPCGeo& curTPCGeom = geometry->PositionToTPC(edep_ptr->MidPoint());
    const auto plane0 = wireReadout->FirstPlane(curTPCGeom.ID());
    double ticksPercm = detPropData.GetXTicksCoefficient(); // this should be by TPCID, but isn't building right now
    double zeroTick = detPropData.ConvertXToTicks(0, plane0.ID());
    auto edep_tick = ticksPercm * edep_ptr->X() + (zeroTick + offset) + tickOffset;
    edep_tick = detPropData.ConvertXToTicks(edep_ptr->X(), plane0.ID()) + offset + tickOffset;

    // loop over subROIs
    unsigned int closest_hit = std::numeric_limits<unsigned int>::max();
    float min_dist = std::numeric_limits<float>::max();
    for (unsigned int i_h = 0; i_h < subROIPropVec.size(); ++i_h)
    {
      auto subroi_prop = subROIPropVec[i_h];
      if (edep_tick > subroi_prop.center-subroi_prop.sigma && edep_tick < subroi_prop.center+subroi_prop.sigma)
      {
        EDepMatchedSubROIMap[i_e].push_back(i_h);
        TrackIDMatchedSubROIMap[edep_ptr->TrackID()].emplace(i_h);
      }
      float hit_dist = std::abs(edep_tick - subroi_prop.center) / subroi_prop.sigma;
      if (hit_dist < min_dist)
      {
        closest_hit = i_h;
        min_dist = hit_dist;
      }
    }

    // if EDep is less than 2.5 units away from closest subROI, assign it to that subROI
    if (min_dist < 5) // try 5 for testing purposes
    {
      auto i_h = closest_hit;
      SubROIMatchedEDepMap[i_h].push_back(i_e);
      SubROIMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] += edep_ptr->E();
    }
  }

  // convert to desired format (possibly a better way to do this...?)
  std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> ReturnMap;
  for (auto it_h = SubROIMatchedEDepMap.begin(); it_h != SubROIMatchedEDepMap.end(); ++it_h)
  {
    auto key = subROIPropVec[it_h->first].key;
    for (auto const& i_e : it_h->second)
    {
      ReturnMap[key].push_back(edepPtrVec[i_e]);
    }
  }

  return ReturnMap;
}


//--- CalcPropertiesFromEdeps ---
sys::WireModUtility::TruthProperties_t sys::WireModUtility::CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec, double offset)
{
  //split the edeps by TrackID
  std::map< int, std::vector<const sim::SimEnergyDeposit*> > edepptrs_by_trkid;
  std::map< int, double > energy_per_trkid;
  for(auto const& edep_ptr : edepPtrVec)
  {
    edepptrs_by_trkid[edep_ptr->TrackID()].push_back(edep_ptr);
    energy_per_trkid[edep_ptr->TrackID()]+=edep_ptr->E();
  }

  int trkid_max     = std::numeric_limits<int>::min();
  double energy_max = std::numeric_limits<double>::min();
  for(auto const& e_p_id : energy_per_trkid)
  {
    if(e_p_id.second > energy_max)
    {
      trkid_max = e_p_id.first;
      energy_max = e_p_id.second;
    }
  }

  auto edepPtrVecMaxE = edepptrs_by_trkid[trkid_max];
  
  //first, let's loop over all edeps and get an average weight scale...
  sys::WireModUtility::TruthProperties_t edep_props;
  double total_energy_all = 0.0;
  sys::WireModUtility::ScaleValues_t scales_e_weighted[3];
  for(size_t i_p = 0; i_p < 3; ++i_p)
  {
    scales_e_weighted[i_p].r_Q     = 0.0;
    scales_e_weighted[i_p].r_sigma = 0.0;
  }

  for(auto const edep_ptr : edepPtrVec)
  {
    if (edep_ptr->StepLength() == 0)
      continue;

    edep_props.x = edep_ptr->X();
    edep_props.y = edep_ptr->Y();
    edep_props.z = edep_ptr->Z();

    edep_props.dxdr = (edep_ptr->EndX() - edep_ptr->StartX()) / edep_ptr->StepLength();
    edep_props.dydr = (edep_ptr->EndY() - edep_ptr->StartY()) / edep_ptr->StepLength();
    edep_props.dzdr = (edep_ptr->EndZ() - edep_ptr->StartZ()) / edep_ptr->StepLength();

    edep_props.dedr = edep_ptr->E() / edep_ptr->StepLength();

    total_energy_all += edep_ptr->E();

    const geo::TPCGeo& curTPCGeom = geometry->PositionToTPC(edep_ptr->MidPoint());
    for (auto const& plane : wireReadout->Iterate<geo::PlaneGeo>(curTPCGeom.ID())) {
      int i_p = plane.ID().Plane;
      auto scales = GetViewScaleValues(edep_props, plane.View());
      scales_e_weighted[i_p].r_Q     += edep_ptr->E()*scales.r_Q;
      scales_e_weighted[i_p].r_sigma += edep_ptr->E()*scales.r_sigma; 
    }
  }

  for(size_t i_p = 0; i_p < 3; ++i_p)
  {
    if (total_energy_all > 0)
    {
      scales_e_weighted[i_p].r_Q     = scales_e_weighted[i_p].r_Q / total_energy_all;
      scales_e_weighted[i_p].r_sigma = scales_e_weighted[i_p].r_sigma / total_energy_all;
    }
    if (scales_e_weighted[i_p].r_Q == 0)
      scales_e_weighted[i_p].r_Q = 1;
    if (scales_e_weighted[i_p].r_sigma == 0)
      scales_e_weighted[i_p].r_sigma = 1;
  }

  TruthProperties_t edep_col_properties;

  //copy in the scales that were calculated before
  for(size_t i_p = 0; i_p < 3; ++i_p)
  {
    edep_col_properties.scales_avg[i_p].r_Q = scales_e_weighted[i_p].r_Q;
    edep_col_properties.scales_avg[i_p].r_sigma = scales_e_weighted[i_p].r_sigma;
  }

  // calculations happen here
  edep_col_properties.x              = 0.;
  edep_col_properties.x_rms          = 0.;
  edep_col_properties.x_rms_noWeight = 0.;
  edep_col_properties.x_min          = std::numeric_limits<float>::max();
  edep_col_properties.x_max          = std::numeric_limits<float>::min();

  edep_col_properties.y    = 0.;
  edep_col_properties.z    = 0.;
  edep_col_properties.dxdr = 0.;
  edep_col_properties.dydr = 0.;
  edep_col_properties.dzdr = 0.;

  edep_col_properties.dedr = 0.;

  double total_energy = 0.0;
  for (auto const& edep_ptr : edepPtrVecMaxE)
  {
    edep_col_properties.x += edep_ptr->X()*edep_ptr->E();
    edep_col_properties.x_min = (edep_ptr->X() < edep_col_properties.x_min) ? edep_ptr->X() : edep_col_properties.x_min;
    edep_col_properties.x_max = (edep_ptr->X() > edep_col_properties.x_max) ? edep_ptr->X() : edep_col_properties.x_max;
    total_energy += edep_ptr->E();
    edep_col_properties.y += edep_ptr->Y()*edep_ptr->E();
    edep_col_properties.z += edep_ptr->Z()*edep_ptr->E();

    if (edep_ptr->StepLength() == 0)
      continue;

    edep_col_properties.dxdr += edep_ptr->E()*(edep_ptr->EndX() - edep_ptr->StartX()) / edep_ptr->StepLength();
    edep_col_properties.dydr += edep_ptr->E()*(edep_ptr->EndY() - edep_ptr->StartY()) / edep_ptr->StepLength();
    edep_col_properties.dzdr += edep_ptr->E()*(edep_ptr->EndZ() - edep_ptr->StartZ()) / edep_ptr->StepLength();

    edep_col_properties.dedr += edep_ptr->E()*edep_ptr->E() / edep_ptr->StepLength();
  }
  if (total_energy > 0)
  {
    edep_col_properties.x    = edep_col_properties.x / total_energy;
    edep_col_properties.y    = edep_col_properties.y / total_energy;
    edep_col_properties.z    = edep_col_properties.z / total_energy;
    edep_col_properties.dxdr = edep_col_properties.dxdr / total_energy;
    edep_col_properties.dydr = edep_col_properties.dydr / total_energy;
    edep_col_properties.dzdr = edep_col_properties.dzdr / total_energy;

    edep_col_properties.dedr = edep_col_properties.dedr / total_energy;
  }

  for (auto const& edep_ptr : edepPtrVecMaxE)
  {
    edep_col_properties.x_rms          += (edep_ptr->X()-edep_col_properties.x)*(edep_ptr->X()-edep_col_properties.x)*edep_ptr->E();
    edep_col_properties.x_rms_noWeight += (edep_ptr->X()-edep_col_properties.x)*(edep_ptr->X()-edep_col_properties.x);
  }
  edep_col_properties.x_rms_noWeight = std::sqrt(edep_col_properties.x_rms_noWeight);

  if (total_energy > 0)
    edep_col_properties.x_rms = std::sqrt(edep_col_properties.x_rms/total_energy);

  const geo::TPCGeo& tpcGeom = geometry->PositionToTPC({edep_col_properties.x, edep_col_properties.y, edep_col_properties.z});
  const auto plane0 = wireReadout->FirstPlane(tpcGeom.ID());
  double ticksPercm = detPropData.GetXTicksCoefficient(); // this should be by TPCID, but isn't building right now
  edep_col_properties.tick              = detPropData.ConvertXToTicks(edep_col_properties.x    , plane0.ID()) + offset + tickOffset;
  edep_col_properties.tick_rms          = ticksPercm*edep_col_properties.x_rms;
  edep_col_properties.tick_rms_noWeight = ticksPercm*edep_col_properties.x_rms_noWeight;
  edep_col_properties.tick_min          = detPropData.ConvertXToTicks(edep_col_properties.x_min, plane0.ID()) + offset + tickOffset;
  edep_col_properties.tick_max          = detPropData.ConvertXToTicks(edep_col_properties.x_max, plane0.ID()) + offset + tickOffset;
  edep_col_properties.total_energy      = total_energy;
  
  return edep_col_properties; 
}

//--- GetScaleValues ---
sys::WireModUtility::ScaleValues_t sys::WireModUtility::GetScaleValues(sys::WireModUtility::TruthProperties_t const& truth_props, sys::WireModUtility::ROIProperties_t const& roi_vals)
{
  sys::WireModUtility::ScaleValues_t scales;
  sys::WireModUtility::ScaleValues_t channelScales = GetChannelScaleValues(truth_props, roi_vals.channel);
  sys::WireModUtility::ScaleValues_t viewScales    = GetViewScaleValues(truth_props, roi_vals.view);
  scales.r_Q     = channelScales.r_Q     * viewScales.r_Q;
  scales.r_sigma = channelScales.r_sigma * viewScales.r_sigma;
  return scales;
}

//--- GetChannelScaleValues ---
// Rescaling depending only on the channel (noise, gain systematics etc.) To modify for DUNE
sys::WireModUtility::ScaleValues_t sys::WireModUtility::GetChannelScaleValues(sys::WireModUtility::TruthProperties_t const& truth_props, raw::ChannelID_t const& channel)
{
  // initialize return
  sys::WireModUtility::ScaleValues_t scales;
  scales.r_Q     = 1.0;
  scales.r_sigma = 1.0;
  
  // try to get geo
  //   // if not in a TPC return default values
  double const truth_coords[3] = {truth_props.x, truth_props.y, truth_props.z};
  geo::TPCGeo const* curTPCGeomPtr = geometry->PositionToTPCptr(geo::vect::makePointFromCoords(truth_coords));
  if (curTPCGeomPtr == nullptr)
    return scales;

  if (applyGainScale)
  {
    scales.r_Q *= gainScale; //changing only amplitude for electronics gain. Should I also change sigma?
  }
  return scales;
}

//--- GetViewScaleValues ---
// Rescaling depending on truth info of the event (recombination, attenuation etc.). To modify for DUNE using analytical formula
sys::WireModUtility::ScaleValues_t sys::WireModUtility::GetViewScaleValues(sys::WireModUtility::TruthProperties_t const& truth_props, geo::View_t const& view)
{
  // initialize return
  sys::WireModUtility::ScaleValues_t scales;
  scales.r_Q     = 1.0;
  scales.r_sigma = 1.0;
  
  // try to get geo
  //   // if not in a TPC return default values
  double const truth_coords[3] = {truth_props.x, truth_props.y, truth_props.z};
  geo::TPCGeo const* curTPCGeomPtr = geometry->PositionToTPCptr(geo::vect::makePointFromCoords(truth_coords));
  if (curTPCGeomPtr == nullptr)
    return scales;

  //empty for now. Will then apply scaling for recombination, attenuation etc. here
  return scales;
}

//--- ModifyROI ---
void sys::WireModUtility::ModifyROI(std::vector<float> & roi_data,
                                    sys::WireModUtility::ROIProperties_t const& roi_prop,
                                    std::vector<sys::WireModUtility::SubROIProperties_t> const& subROIPropVec,
                                    std::map<sys::WireModUtility::SubROI_Key_t, sys::WireModUtility::ScaleValues_t> const& subROIScaleMap)
{
  // do you want a bunch of messages?
  bool verbose = false;
  
  // initialize some values
  double q_orig = 0.0;
  double q_mod  = 0.0;
  double scale_ratio = 1.0;

  // loop over the ticks
  for(size_t i_t = 0; i_t < roi_data.size(); ++i_t)
  {
    // reset your values
    q_orig = 0.0;
    q_mod  = 0.0;
    scale_ratio = 1.0;

    // loop over the subs
    for (auto const& subroi_prop : subROIPropVec)
    {
      // get your scale vals
      auto scale_vals = subROIScaleMap.find(subroi_prop.key)->second;

      q_orig += gausFunc(i_t + roi_prop.begin, subroi_prop.center,                      subroi_prop.sigma,                  subroi_prop.total_q);
      q_mod  += gausFunc(i_t + roi_prop.begin, subroi_prop.center, scale_vals.r_sigma * subroi_prop.sigma, scale_vals.r_Q * subroi_prop.total_q);

      if (verbose)
        std::cout << "    Incrementing q_orig by gausFunc(" << i_t + roi_prop.begin << ", " << subroi_prop.center << ", " <<                      subroi_prop.sigma << ", " <<                  subroi_prop.total_q << ")" << '\n'
                  << "    Incrementing q_mod  by gausFunc(" << i_t + roi_prop.begin << ", " << subroi_prop.center << ", " << scale_vals.r_sigma * subroi_prop.sigma << ", " << scale_vals.r_Q * subroi_prop.total_q << ")" << std::endl;
    }

    // do some sanity checks
    if        (isnan(q_orig))
    {
      if (verbose)
        std::cout << "WARNING: obtained q_orig = NaN... setting scale to 1" << std::endl;
      scale_ratio = 1.0;
    } else if (isnan(q_mod)) {
      if (verbose)
        std::cout << "WARNING: obtained q_mod = NaN... setting scale to 0" << std::endl;
      scale_ratio = 0.0;
    } else {
      scale_ratio = q_mod / q_orig;
    }
    if(isnan(scale_ratio) || isinf(scale_ratio))
    {
      if (verbose)
        std::cout << "WARNING: obtained scale_ratio = " << q_mod << " / " << q_orig << " = NAN/Inf... setting to 1" << std::endl;
      scale_ratio = 1.0;
    }

    roi_data[i_t] = scale_ratio * roi_data[i_t];
    if (verbose)
      std::cout << "\t tick " << i_t << ":"
                <<  " data="   << roi_data[i_t]
                << ", q_orig=" << q_orig
                << ", q_mod="  << q_mod
                << ", ratio="  << scale_ratio << std::endl;
  }

  // we're done now
  return;
}
