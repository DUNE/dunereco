#include "FDSelectionUtils.h"

double FDSelectionUtils::CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& selected_hits, const std::vector<art::Ptr<recob::Hit> >& all_hits, int track_id){
  int num_matches_in_sel_hits = 0;
  int num_matches_in_all_hits = 0;
  for (unsigned int i_hit = 0; i_hit < selected_hits.size(); i_hit++){
    int matched_id = TruthMatchUtils::TrueParticleID(clockData, selected_hits[i_hit], 1); 
    if (matched_id==track_id) num_matches_in_sel_hits++;
  }
  for (unsigned int i_hit = 0; i_hit < all_hits.size(); i_hit++){
    int matched_id = TruthMatchUtils::TrueParticleID(clockData, all_hits[i_hit], 1); 
    if (matched_id==track_id) num_matches_in_all_hits++;
  }

  double completeness = 0;
  if (num_matches_in_all_hits > 0) completeness = 1.*num_matches_in_sel_hits/num_matches_in_all_hits;
  return completeness;
}

double FDSelectionUtils::HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& selected_hits, int track_id){
  int num_matches_in_sel_hits = 0;
  for (unsigned int i_hit = 0; i_hit < selected_hits.size(); i_hit++){
    int matched_id = TruthMatchUtils::TrueParticleID(clockData, selected_hits[i_hit], 1); 
    if (matched_id==track_id) num_matches_in_sel_hits++;
  }
  double hit_purity = 0;
  if (selected_hits.size() > 0) hit_purity = 1.*num_matches_in_sel_hits/selected_hits.size();
  return hit_purity;
}

bool FDSelectionUtils::IsInsideTPC(TVector3 position, double distance_buffer){
  geo::Point_t vtx(position.X(), position.Y(), position.Z());
  bool inside = false;
  art::ServiceHandle<geo::Geometry> geom;
  geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

  if (geom->HasTPC(idtpc))
  {
    const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

    //need to loop over cryostats?

    for (size_t cIndex = 0; cIndex < geom->Ncryostats(); cIndex++)
    {
      geo::CryostatID c(cIndex);

      const geo::CryostatGeo& cryostat = geom->Cryostat(c);
      for (size_t t = 0; t < cryostat.NTPC(); t++)
      {
        const geo::TPCGeo& tpcg = cryostat.TPC(t);
        if (tpcg.MinX() < minx) minx = tpcg.MinX();
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
        if (tpcg.MinY() < miny) miny = tpcg.MinY();
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
      }
    }

    //x
    double dista = fabs(minx - position.X());
    double distb = fabs(position.X() - maxx);
    if ((position.X() > minx) && (position.X() < maxx) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    //y
    dista = fabs(maxy - position.Y());
    distb = fabs(position.Y() - miny);
    if (inside && (position.Y() > miny) && (position.Y() < maxy) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    else inside = false;
    //z
    dista = fabs(maxz - position.Z());
    distb = fabs(position.Z() - minz);
    if (inside && (position.Z() > minz) && (position.Z() < maxz) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    else inside = false;
  }

  return inside;
}

double FDSelectionUtils::CalculateTrackLength(const art::Ptr<recob::Track> track){
  double length = 0;
  if (track->NumberTrajectoryPoints()==1) return length; //Nothing to calculate if there is only one point

  for (size_t i_tp = 0; i_tp < track->NumberTrajectoryPoints()-1; i_tp++){ //Loop from the first to 2nd to last point
    TVector3 this_point(track->TrajectoryPoint(i_tp).position.X(),track->TrajectoryPoint(i_tp).position.Y(),track->TrajectoryPoint(i_tp).position.Z());
    if (!FDSelectionUtils::IsInsideTPC(this_point,0)){
      std::cout<<"FDSelectionUtils::CalculateTrackLength - Current trajectory point not in the TPC volume.  Skip over this point in the track length calculation"<<std::endl;
      continue;
    }
    TVector3 next_point(track->TrajectoryPoint(i_tp+1).position.X(),track->TrajectoryPoint(i_tp+1).position.Y(),track->TrajectoryPoint(i_tp+1).position.Z());
    if (!FDSelectionUtils::IsInsideTPC(next_point,0)){
      std::cout<<"FDSelectionUtils::CalculateTrackLength - Next trajectory point not in the TPC volume.  Skip over this point in the track length calculation"<<std::endl;
      continue;
    }

    length+=(next_point-this_point).Mag();
  }
  return length;
}



