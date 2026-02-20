#include <vector>

#include "TFile.h"
#include "TSpline.h"
#include "TGraph2D.h"
#include "TNtuple.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <limits>
#include <memory>

namespace sys {
  class WireModUtility{
    public:

      const geo::GeometryCore* geometry;                  // save the TPC geometry
      const geo::WireReadoutGeom* wireReadout;            // new for LarSoft v10
      const detinfo::DetectorPropertiesData& detPropData; // save the detector property data
      bool applyGainScale; //for electronics gain systematics
      double lifetime_var; //new lifetime value for given variation, in microseconds
      double readoutWindowTicks;                          // how many ticks are in the readout window?
      double tickOffset;                                  // do we want an offset in the ticks?

      //removed members related to splines

      //graphs could still be needed, to be tuned
      std::vector<TGraph2D*> graph2Ds_Charge_XXZAngle;    // the graphs for the charge correction in X vs XZ angle
      std::vector<TGraph2D*> graph2Ds_Sigma_XXZAngle;     // the graphs for the width correction in X vs XZ angle
      std::vector<TGraph2D*> graph2Ds_Charge_XdQdX;       // the graphs for charge correction in X vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Sigma_XdQdX;        // the graphs for width correction in X vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Charge_XZAngledQdX; // the graphs for charge correction in XZ angle vs dQ/dX
      std::vector<TGraph2D*> graph2Ds_Sigma_XZAngledQdX;  // the graphs for width correction in XZ angle vs dQ/dX

      //for debugging purposes: adding a map between the trackID and the PDG code
      std::map<int,int> trackID_to_PDG;
      unsigned int event_counter; 
   
    //constructor
      WireModUtility(const geo::GeometryCore* geom, 
                     const geo::WireReadoutGeom* wireRead,
                     const detinfo::DetectorPropertiesData& detProp,
                     const bool& arg_ApplyGainScale = false,
                     const double& arg_gainScale = 1, //relative gain scale. Default is 1, no bias.
                     const double& arg_TickOffset = 0)
      : geometry(geom),
        wireReadout(wireRead),
        detPropData(detProp),
        applyGainScale(arg_ApplyGainScale),
        gainScale(arg_gainScale),
        readoutWindowTicks(detProp.ReadOutWindowSize()),                                               // the default A2795 (ICARUS TPC readout board) readout window is 4096 samples
        tickOffset(arg_TickOffset),                                                                     // tick offset is for MC truth, default to zero and set only as necessary
        event_counter(0)
      {
      }
      
      typedef std::pair<unsigned int,unsigned int>  ROI_Key_t;
      typedef std::pair<ROI_Key_t, unsigned int> SubROI_Key_t;
  
      typedef struct ROIProperties
      {
        ROI_Key_t key;
        raw::ChannelID_t channel;
        geo::View_t view;
        float begin;
        float end;
        float total_q;
        float center;   //charge weighted center of ROI
        float sigma;    //charge weighted RMS of ROI
      } ROIProperties_t;

      typedef struct SubROIProperties
      {
        SubROI_Key_t key;
        raw::ChannelID_t channel;
        geo::View_t view;
        float total_q;
        float center;
        float sigma;
      } SubROIProperties_t;

      typedef struct ScaleValues
      {
        double r_Q;
        double r_sigma;
      } ScaleValues_t;

      typedef struct TruthProperties
      {
        float x;
        float x_rms;
        float x_rms_noWeight;
        float tick;
        float tick_rms;
        float tick_rms_noWeight;
        float total_energy;
        float x_min;
        float x_max;
        float tick_min;
        float tick_max;
        float y;
        float z;
        //Direction cosines
        double dxdr;
        double dydr;
        double dzdr;
        double dqdr;
        double dedr;
        ScaleValues_t scales_avg[3];
      } TruthProperties_t;
    
      std::map< ROI_Key_t,std::vector<size_t> > ROIMatchedEdepMap;
      std::map< ROI_Key_t,std::vector<size_t> > ROIMatchedHitMap;


      //useful functions
      //geometry functions
      double planeXToTick(double xPos, const geo::PlaneGeo& plane, const geo::TPCGeo& tpcGeom, double offset = 0) {
          return detPropData.ConvertXToTicks(xPos, plane.ID()) + offset;
      }

      bool planeXInWindow(double xPos, const geo::PlaneGeo& plane, const geo::TPCGeo& tpcGeom, double offset = 0)
      {
        double tick = planeXToTick(xPos, plane, tpcGeom, offset);
        return (tick > 0 && tick <= detPropData.ReadOutWindowSize());
      } 

      double gausFunc(double t, double mean, double sigma, double a = 1.0)
      {
        return (a / (sigma * std::sqrt(2 * util::pi()))) * std::exp(-0.5 * std::pow((t - mean)/sigma, 2));
      }

      //Folds angle in [0,90deg]
      double FoldAngle(double theta)
      {
        return (std::abs(theta) > 0.5 * util::pi()) ? util::pi() - std::abs(theta) : std::abs(theta);
      }


      //Functions to calculte fold angle related to plane. Needed if planes have a certain angle. planeAngle in degrees!
      double ThetaXZ_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double planeAngleRad = planeAngle * (util::pi() / 180.0);
        double sinPlaneAngle = std::sin(planeAngleRad);
        double cosPlaneAngle = std::cos(planeAngleRad);

        double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;
        
        double theta = std::atan2(dxdr, dzdrPlaneRel);
        return FoldAngle(theta);
       }

      double ThetaYZ_PlaneRel(double dxdr, double dydr, double dzdr, double planeAngle)
      {
        double planeAngleRad = planeAngle * (util::pi() / 180.0);
        double sinPlaneAngle = std::sin(planeAngleRad);
        double cosPlaneAngle = std::cos(planeAngleRad);

        double dydrPlaneRel = dydr * cosPlaneAngle - dzdr * sinPlaneAngle;
        double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;

        double theta = std::atan2(dydrPlaneRel, dzdrPlaneRel);
        return FoldAngle(theta);
      }

      //Defined in the .cc
      void MakeTrackIDtoPDGMap(const std::vector<simb::MCParticle>& mcParticles);      

      ROIProperties_t CalcROIProperties(recob::Wire const&, size_t const&);

      std::vector<std::pair<unsigned int, unsigned int>> GetTargetROIs(sim::SimEnergyDeposit const&, double offset);
      std::vector<std::pair<unsigned int, unsigned int>> GetHitTargetROIs(recob::Hit const&);

      void FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const&, std::vector<recob::Wire> const&, double offset);
      void FillROIMatchedHitMap(std::vector<recob::Hit> const&, std::vector<recob::Wire> const&);

      std::vector<SubROIProperties_t> CalcSubROIProperties(ROIProperties_t const&, std::vector<const recob::Hit*> const&);

      std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> MatchEdepsToSubROIs(std::vector<SubROIProperties_t> const&, std::vector<const sim::SimEnergyDeposit*> const&, double offset);

      TruthProperties_t CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const&, double offset);

      //product of view scale and channel scale
      ScaleValues_t GetScaleValues(TruthProperties_t const&, ROIProperties_t const&);
      //scaling factors not dependent on thruth properties (e.g., electronics gain)
      ScaleValues_t GetChannelScaleValues(TruthProperties_t const&, raw::ChannelID_t const&);
      //scaling factors dependent on truth properties (e.g., recombination, attenuation)
      ScaleValues_t GetViewScaleValues(TruthProperties_t const&, geo::View_t const&);

      void ModifyROI(std::vector<float> &,
                     ROIProperties_t const &,
                     std::vector<SubROIProperties_t> const&,
                     std::map<SubROI_Key_t, ScaleValues_t> const&);
  }; // end class
} // end namespace
