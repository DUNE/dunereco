// std inlcudes
#include <string>
#include <vector>

//ROOT includes
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "Utility/WireModUtility.hh"

//namespace
namespace wiremod
{

  class WireModifier : public art::EDProducer
  {
    public:
      explicit WireModifier(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

   private:
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
      const geo::WireReadoutGeom* fWireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());
      std::string fRatioFileName; // there is where we try to grab the splines/graphs (if they exist)
      //TODO: Not using splines here
      TSpline3*              fSpline_charge_Channel;
      TSpline3*              fSpline_sigma_Channel;
      std::vector<TSpline3*> fSpline_charge_X;
      std::vector<TSpline3*> fSpline_sigma_X;
      std::vector<TSpline3*> fSpline_charge_XZAngle;
      std::vector<TSpline3*> fSpline_sigma_XZAngle;
      std::vector<TSpline3*> fSpline_charge_YZAngle;
      std::vector<TSpline3*> fSpline_sigma_YZAngle;
      std::vector<TSpline3*> fSpline_charge_dEdX;
      std::vector<TSpline3*> fSpline_sigma_dEdX;
      std::vector<TGraph2D*> fGraph_charge_YZ; 
      std::vector<TGraph2D*> fGraph_sigma_YZ;
      art::InputTag fWireLabel;     // which wires are we pulling in?
      art::InputTag fHitLabel;      // which hits are we pulling in?
      art::InputTag fEDepOrigLabel; // which are the unshifted EDeps?
      art::InputTag fEDepShftLabel; // which are the shifted EDeps?
      bool fApplyGainScale;
      double fGainScale;
      bool fApplyLifetimeVar;
      double fLifetimeVar;
      bool fApplyModBoxVar;
      double fModBoxAlphaVar;
      double fModBoxBetaVar;
      bool fApplyLongitudinalDiffusion;
      double fDLnew;
      bool fSaveChargeRatioPlots; //plot showing the ratio of ROI total charge after/before modification
      TH2F *hRat = nullptr;
      TGraph *grChargeRatExp = nullptr;
      //The three following are probably not needed
      bool fSaveHistsByChannel;     // save modified signals by channel?
      bool fSaveHistsByWire;        // save modified signals by wire?
      bool fIsData;                 // DEBUG: get data wires

  }; // end WireModifier class

  WireModifier::WireModifier(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  void WireModifier::reconfigure(fhicl::ParameterSet const& pset)
  {
    fWireLabel     = pset.get<art::InputTag>("WireLabel", "tpcrawdecoder:gauss");
    fHitLabel      = pset.get<art::InputTag>("HitLabel", "gaushit");
    fEDepOrigLabel = pset.get<art::InputTag>("EDepOrigLabel", "IonAndScint");
    fEDepShftLabel = pset.get<art::InputTag>("EDepShftLabel", "largeant:LArG4DetectorServicevolTPCActiveInner");

    // what, if anything, are we putting in the histogram files
    fSaveHistsByChannel = pset.get<bool>("SaveByChannel", false);
    //fSaveHistsByChannel = pset.get<bool>("SaveByChannel", true);
    fSaveHistsByWire    = pset.get<bool>("SaveByWire"   , false);
    fIsData             = pset.get<bool>("IsData"       , false);

    // try to read in the graphs/splines from a file
    //     // if that file does not exist then fake them
    fRatioFileName = pset.get<std::string>("RatioFileName", "NOFILE");
    fApplyGainScale    = pset.get<bool>("ApplyGainScale"   , false);
    if (fApplyGainScale)
    {
      fGainScale = pset.get<double>("ElectronicsGainScale");
    }
    fApplyLifetimeVar    = pset.get<bool>("ApplyLifetimeVar"   , false);
    if (fApplyLifetimeVar)
    {
      fLifetimeVar = pset.get<double>("LifetimeVar");
    }
    fApplyModBoxVar    = pset.get<bool>("ApplyModBoxVar"   , false);
    if (fApplyModBoxVar)
    {
      fModBoxAlphaVar = pset.get<double>("ModBoxAlphaVar");
      fModBoxBetaVar = pset.get<double>("ModBoxBetaVar");
    }
    fApplyLongitudinalDiffusion = pset.get<bool>("ApplyLongitudinalDiffusionVar"   , false);
    if (fApplyLongitudinalDiffusion){
      fDLnew = pset.get<double>("DLvar");
    }

    fSaveChargeRatioPlots = pset.get<bool>("SaveChargeRatioPlot"   , false);
    if (!fApplyModBoxVar and !fApplyLifetimeVar) fSaveChargeRatioPlots = false; //for now only consider recombination and attenuation for charge ration plots
    if (fSaveChargeRatioPlots){
      const int N=200;
      double xs[N], ys[N];
      if (fApplyModBoxVar){ 
        hRat = new TH2F("hRatRecomb", "Total ROI charge after/before scaling vs dE/dx; dE/dx [MeV/cm]; Charge ratio", 200, 1, 3, 200, 0, 1);
        double beta_nom=0.212;
        double alpha_nom=0.93;
        double E=0.5;
        double rho=1.39295;
        double Bnom=beta_nom/rho/E;
        double Bvar=fModBoxBetaVar/rho/E;
        double xmin=1., step=2./N;
        for (int i = 0; i < N; i++) {
          xs[i] = xmin + i * step;
          ys[i] = std::log(fModBoxAlphaVar+Bvar*xs[i])/std::log(alpha_nom+Bnom*xs[i])*beta_nom/fModBoxBetaVar;
        }
      }
      else{ 
        hRat = new TH2F("hRatRecomb", "Total ROI charge after/before scaling vs drift distance; Drift distance [cm]; Charge ratio", 200, 0, 350, 200, 0, 1);   
        double factor=1./160.563*(1./3.-1/10.4);
        double xmin=0., step=350./N;
        for (int i = 0; i < N; i++) {
          xs[i] = xmin + i * step;
          ys[i] = std::exp(-xs[i] * factor);
        }   
      }
      grChargeRatExp = new TGraph(N, xs, ys);   
    }
    // we make these things
    produces<std::vector<recob::Wire      >>();
  }

 

  void WireModifier::produce(art::Event& evt)
  {

    // here's where the "magic" happens
    art::ServiceHandle<art::TFileService> tfs;

    // get a clock and det props for the event
    const detinfo::DetectorPropertiesData detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);


    // get the things to do the things on
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel(fWireLabel, wireHandle);
    auto const& wireVec(*wireHandle);

    art::FindManyP<raw::RawDigit> digit_assn(wireHandle, evt, fWireLabel);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepOrigHandle;
    evt.getByLabel(fEDepOrigLabel, edepOrigHandle);
    auto const& edepOrigVec(*edepOrigHandle);
      
    art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
    evt.getByLabel(fEDepShftLabel, edepShiftedHandle);
    auto const& edepShiftedVec(*edepShiftedHandle);

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitLabel, hitHandle);
    auto const& hitVec(*hitHandle);

    //for debugging purposes, to know if event is CC or NC
    auto const& mctruths = *evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    bool isCC = false;
    if (!mctruths.empty()) {
      isCC = (mctruths[0].GetNeutrino().CCNC() == 0);
    }
    if (isCC) std::cout<<"This is a CC events"<<std::endl;

    // put the new stuff somewhere
    std::unique_ptr<std::vector<recob::Wire      >> new_wires(new std::vector<recob::Wire      >());


    //TODO: modify for DUNE implementation! (correct flags, correct if loops etc.)
    sys::WireModUtility wmUtil(fGeometry, fWireReadout, detProp); // detector geometry & properties

    //Get the simulated particles to get the PDG ID for debugging purposes
    /*art::Handle<std::vector<simb::MCParticle>> mcPartHandle;
    evt.getByLabel("largeant", mcPartHandle);  // adjust label if needed
    auto const& mcParticles(*mcPartHandle);
    wmUtil.MakeTrackIDtoPDGMap(mcParticles); //map is created and will be used to get Edeps PDG plot
    */

    wmUtil.applyGainScale    = fApplyGainScale;
    if (wmUtil.applyGainScale)
    {
      std::cout<<"Will apply gain scale"<<std::endl;
      wmUtil.gainScale = fGainScale;
    }
    wmUtil.applyLifetimeVar    = fApplyLifetimeVar;
    if (wmUtil.applyLifetimeVar)
    {
      std::cout<<"Will apply lifetime variation scaling"<<std::endl;
      wmUtil.lifetime_var = fLifetimeVar;
    }
    wmUtil.applyModBoxVar = fApplyModBoxVar;
    if (wmUtil.applyModBoxVar)
    {
      std::cout<<"Will apply modified box variation scaling"<<std::endl;
      wmUtil.ModBoxAlphaVar = fModBoxAlphaVar;
      wmUtil.ModBoxBetaVar = fModBoxBetaVar;
    }
    wmUtil.applyLongitudinalDiffusion = fApplyLongitudinalDiffusion;
    if (wmUtil.applyLongitudinalDiffusion){
      std::cout<<"Will apply longitudinal diffusion scaling"<<std::endl;
      wmUtil.DLnew = fDLnew;
    }
    // add some debugging here
    MF_LOG_VERBATIM("WireModifier")
      << "DUMP CONFIG:" << '\n'
      << "---------------------------------------------------" << '\n'
      << "  applyGainScale:     " << wmUtil.applyGainScale     << '\n'
      << "  applyLifetimeVar:     " << wmUtil.applyLifetimeVar     << '\n'
      << "  readoutWindowTicks: " << wmUtil.readoutWindowTicks << '\n'
      << "  tickOffset:         " << wmUtil.tickOffset         << '\n'
      << "---------------------------------------------------";

    // do the things
    double offset_ADC = 0; // don't use an offset atm
    MF_LOG_VERBATIM("WireModifier")
      << "Get Edep Map";
    //std::cout<<"Total number of shifted Edeps: "<<edepShiftedVec.size()<<std::endl;
    wmUtil.FillROIMatchedEdepMap(edepShiftedVec, wireVec, offset_ADC);
    MF_LOG_VERBATIM("WireModifier")
      << "Got Edep Map." << '\n'
      << "Get Hit Map";
    wmUtil.FillROIMatchedHitMap(hitVec, wireVec);
    MF_LOG_VERBATIM("WireModifier")
      << "Got Hit Map.";


    int nROIs=0;
    int nROIs_mod=0;
    int nROIs_lowQ_mod=0;
    int nROIs_lowQ=0;
    int nROIs_hit=0;
    int nROIs_hit_mod=0;
    int nROIs_hit_highQ=0;
    int nROIs_hit_highQ_mod=0;



    // loop-de-loop
    for(size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    {
      MF_LOG_DEBUG("WireModifier")
        << "Checking wire " << i_w;

      auto const& wire = wireVec.at(i_w);


      recob::Wire::RegionsOfInterest_t new_rois;
      new_rois     .resize(wire.SignalROI().size());

      unsigned int my_plane = geo::kUnknown;
      if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 0)).View())
      {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on plane 0, view " << wire.View();
        my_plane = 0;
      } else if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 1)).View()) {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on plane 1, view " << wire.View();
        my_plane = 1;
      } else if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 2)).View()) {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on plane 2, view " << wire.View();
        my_plane = 2;
      }

      if (my_plane == geo::kUnknown)
      {
        //nNoPlane++;
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on unsupported plane. Skip.";
      }

      // keep track of if this wire is modified
      bool isModified = false;
      
      for(size_t i_r = 0; i_r < wire.SignalROI().get_ranges().size(); ++i_r)
      {
        nROIs++;
        MF_LOG_DEBUG("WireModifier")
          << "  Checking ROI " << i_r;
        auto const& range = wire.SignalROI().get_ranges()[i_r];
        sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(), i_r);

        std::vector<float> modified_data(range.data());

        //The following line can be moved back later (no need for ROI properties if it is not matched): need roi_properties of unmatched ROIs for debigging purposes
        auto roi_properties = wmUtil.CalcROIProperties(wire, i_r);
        bool hasHighQ = true;
        if (roi_properties.total_q<80){ 
          nROIs_lowQ++;
          hasHighQ = false;
        }
     
        bool hasHit = false;
        //The following line can be moved back later: putting it here to see if the ROI has hit before checking if it has matched edep
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key); 
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ){
          nROIs_hit++;
          if (hasHighQ) nROIs_hit_highQ++;
          hasHit = true;
        } 
    
        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          new_rois     .add_range(range.begin_index(), modified_data);
          MF_LOG_DEBUG("WireModifier")
            << "    Could not find matching Edep. Skip";
          //nNoMatch++;
          //hUnMatched->Fill(roi_properties.end-roi_properties.begin, roi_properties.total_q);
          //hUnMatched->Fill(roi_properties.sigma, roi_properties.total_q);
          continue;
        }
        //hMatched->Fill(roi_properties.sigma, roi_properties.total_q);
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          new_rois     .add_range(range.begin_index(), modified_data);
          MF_LOG_DEBUG("WireModifier")
            << "    No indices for Edep. Skip";
          //nNoIndex++;
          continue;
        }
        std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
        std::vector<const sim::SimEnergyDeposit*> matchedShiftedEdepPtrVec;
        for(auto i_e : matchedEdepIdxVec)
        {
          matchedEdepPtrVec.push_back(&edepOrigVec[i_e]);
          matchedShiftedEdepPtrVec.push_back(&edepShiftedVec[i_e]);
        }
        MF_LOG_DEBUG("WireModifier")
          << "  Found " << matchedShiftedEdepPtrVec.size() << " shifted Edeps";

        std::vector<const recob::Hit*> matchedHitPtrVec;
        //auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }

        MF_LOG_DEBUG("WireModifier")
          << "    Found " << matchedHitPtrVec.size() << " matching hits";

        //Moving the following line before to have roi_properties also for unmatched ROIs
        //auto roi_properties = wmUtil.CalcROIProperties(wire, i_r);
        MF_LOG_DEBUG("WireModifier")
          << "    ROI Properties:" << '\n'
          << "                    key:     (" << roi_properties.key.first << ", " << roi_properties.key.second << ")" << '\n'
          << "                    view:    " << roi_properties.view << '\n'
          << "                    begin:   " << roi_properties.begin << '\n'
          << "                    end:     " << roi_properties.end << '\n'
          << "                    total_q: " << roi_properties.total_q << '\n'
          << "                    center:  " << roi_properties.center << '\n'
          << "                    sigma:   " << roi_properties.sigma;

        auto subROIPropVec = wmUtil.CalcSubROIProperties(roi_properties, matchedHitPtrVec);

        MF_LOG_DEBUG("WireModifier")
          << "    have " << subROIPropVec.size() << " SubROT";

        auto SubROIMatchedShiftedEdepMap = wmUtil.MatchEdepsToSubROIs(subROIPropVec, matchedShiftedEdepPtrVec, offset_ADC);
        MF_LOG_DEBUG("WireModifier")
          << "    size of SubROIMatchedShiftedEdepMap: " << SubROIMatchedShiftedEdepMap.size();
        std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> SubROIMatchedEdepMap;
        for ( auto const& key_edepPtrVec_pair : SubROIMatchedShiftedEdepMap ) {
          auto key = key_edepPtrVec_pair.first;
          for ( auto const& shifted_edep_ptr : key_edepPtrVec_pair.second ) {
            for ( unsigned int i_e=0; i_e < matchedShiftedEdepPtrVec.size(); i_e++ ) {
              if ( shifted_edep_ptr == matchedShiftedEdepPtrVec[i_e] ) {
                MF_LOG_DEBUG("WireModifier")
                  << "    found matching shifted Edep!";
                SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
                break;
              }
            }
          }
        }

        MF_LOG_DEBUG("WireModifier")
          << "    size of SubROIMatchedEdepMap: " << SubROIMatchedEdepMap.size();

        std::map<sys::WireModUtility::SubROI_Key_t, sys::WireModUtility::ScaleValues_t> SubROIMatchedScalesMap;
        double drift_distance=0;
        double dedx_avg=0;
        double total_E=0;
        for ( auto const& subroi_prop : subROIPropVec ) {
          //nsubROIs++;
          sys::WireModUtility::ScaleValues_t scale_vals;
          auto key = subroi_prop.key;
          auto key_it =  SubROIMatchedEdepMap.find(key);

          if ( key_it != SubROIMatchedEdepMap.end() && key_it->second.size() > 0 ) {
            auto truth_vals = wmUtil.CalcPropertiesFromEdeps(key_it->second, offset_ADC);
            total_E+=truth_vals.total_energy;
            drift_distance+=truth_vals.total_energy*truth_vals.x;
            dedx_avg+=truth_vals.total_energy*truth_vals.dedx;
            /*if ( truth_vals.total_energy < 0.3 && subroi_prop.total_q > 80 ) {
              scale_vals.r_Q     = 1.;
              scale_vals.r_sigma = 1.;
            } 
            else {*/
              scale_vals = wmUtil.GetScaleValues(truth_vals, roi_properties);
              mf::LogDebug("WireModifier")
                << "Scaling! Q scale: " << scale_vals.r_Q
                << "     sigma scale: " << scale_vals.r_sigma;
              isModified = true;
            //}
          }
          else {
            scale_vals.r_Q     = 1.;
            scale_vals.r_sigma = 1.;
          }
          
          SubROIMatchedScalesMap[key] = scale_vals;
        }
        drift_distance=abs(drift_distance/total_E);
        dedx_avg=dedx_avg/total_E;
        if (isModified){ 
          nROIs_mod++;
          if (roi_properties.total_q<80) nROIs_lowQ_mod++;
          if (hasHit){ 
            nROIs_hit_mod++;
            if (hasHighQ) nROIs_hit_highQ_mod++;
          }
        }
        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap);
        double charge = 0;
        for (auto const& adc : modified_data)
          charge += adc;
        if (fSaveChargeRatioPlots){ 
          if (fApplyLifetimeVar) hRat->Fill(drift_distance, charge/roi_properties.total_q);
          else hRat->Fill(dedx_avg, charge/roi_properties.total_q);
        }
        new_rois     .add_range(roi_properties.begin, modified_data);
      }

        

      new_wires->emplace_back(new_rois,      wire.Channel(), wire.View());

      if (fSaveHistsByChannel && isModified)
      {
        readout::ROPID ropID = fWireReadout->ChannelToROP(wire.Channel());  
        std::string titleStr =  "Cryo-"         + std::to_string(ropID.Cryostat)
                             + "_TPCset-"       + std::to_string(ropID.TPCset)
                             + "_ReadOutPlane-" + std::to_string(ropID.ROP)
                             + "_Channel-"      + std::to_string(wire.Channel());
        TH1F* oldChannelHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        TH1F* newChannelHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
        {
          float oldSample = (tick <     wire         .Signal().size() ) ?     wire         .Signal().at(tick) : 0;
          float newSample = (tick < new_wires->back().Signal().size() ) ? new_wires->back().Signal().at(tick) : 0;
          oldChannelHist->SetBinContent(tick + 1, oldSample);
          newChannelHist->SetBinContent(tick + 1, newSample);
        }
          
        TH1F* oldChannelHist_toSave = tfs->make<TH1F>(*oldChannelHist);
        TH1F* newChannelHist_toSave = tfs->make<TH1F>(*newChannelHist);
        mf::LogDebug("WireModifier")
          << "Saved histograms " << oldChannelHist_toSave->GetName() << '\n'
          << "             and " << newChannelHist_toSave->GetName();
      }

      if (fSaveHistsByWire && isModified)
      {
        std::vector<geo::WireID> wireIDs = fWireReadout->ChannelToWire(wire.Channel());
        mf::LogDebug("WireModifier")
          << "Channel " << wire.Channel() << " has " << wireIDs.size() << " wire(s)";
        for (auto const& wireID : wireIDs)
        {
          std::string titleStr =  "Cryo-"  + std::to_string(wireID.Cryostat)
                               + "_TPC-"   + std::to_string(wireID.TPC)
                               + "_Plane-" + std::to_string(wireID.Plane)
                               + "_Wire-"  + std::to_string(wireID.Wire);
          TH1F* oldWireHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          TH1F* newWireHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
          {
            float oldSample = (tick <     wire         .Signal().size() ) ?     wire         .Signal().at(tick) : 0;
            float newSample = (tick < new_wires->back().Signal().size() ) ? new_wires->back().Signal().at(tick) : 0;
            oldWireHist->SetBinContent(tick + 1, oldSample);
            newWireHist->SetBinContent(tick + 1, newSample);
          }

          TH1F* oldWireHist_toSave = tfs->make<TH1F>(*oldWireHist);
          TH1F* newWireHist_toSave = tfs->make<TH1F>(*newWireHist);
          mf::LogDebug("WireModifier")
            << "Saved histograms " << oldWireHist_toSave->GetName() << '\n'
            << "             and " << newWireHist_toSave->GetName();
        }
      }
    } // end loop over wires
    std::cout<<"Modified "<<nROIs_mod<<" ROIs over a total of "<<nROIs<<std::endl;
    //std::cout<<"Modified "<<nsubROIs_mod<<" subROIs over a total of "<<nsubROIs<<std::endl;
    //std::cout<<nNoMatch<<" ROIs could not be associated to an edep"<<std::endl;
    std::cout<<"Modified "<<nROIs_lowQ_mod<<" low charge ROIs over a total of "<<nROIs_lowQ<<std::endl;
    std::cout<<"Modified "<<nROIs_hit_mod<<" ROIs over a total of "<<nROIs_hit<<" leading to an hit"<<std::endl;
    std::cout<<"Modified "<<nROIs_hit_highQ_mod<<" ROIs over a total of "<<nROIs_hit_highQ<<" ROIS with high charge leading to an hit"<<std::endl;
    //std::cout<<nNoIndex<<" no index found for the edep"<<std::endl;
    /*hMatched->SetMarkerStyle(2);
    hMatched->SetMarkerColor(kBlue);
    hMatched->SetStats(0);
    hUnMatched->SetMarkerStyle(5);
    hUnMatched->SetMarkerColor(kRed);
    hUnMatched->SetStats(0);
    hMatched->Draw("P");
    hUnMatched->Draw("Psame");
    c1->SaveAs("ROI_properties.pdf");*/
    if (fSaveChargeRatioPlots){
      TCanvas *c1 = new TCanvas("c1", "ROI properties", 800, 600);
      grChargeRatExp->SetLineWidth(2);
      grChargeRatExp->SetLineColor(kRed);
      hRat->SetStats(0);
      hRat->SetLineColor(kBlack);
      grChargeRatExp->SetTitle(hRat->GetTitle());
      grChargeRatExp->GetXaxis()->SetTitle(hRat->GetXaxis()->GetTitle());
      grChargeRatExp->GetYaxis()->SetTitle(hRat->GetYaxis()->GetTitle());
      grChargeRatExp->Draw("AL");
      hRat->Draw("Psame");
      TLegend *leg = new TLegend(0.15, 0.12, 0.85, 0.3);
      if (fApplyModBoxVar) leg->SetHeader(Form("#splitline{Nominal Mod. Box A: %.2f, Nominal Mod. Box B: %.3f (g.kV)/(MeV.cm^{2})}{Varied Mod. Box A: %.2f, Varied Mod. Box B: %.3f (g.kV)/(MeV.cm^{2})}", 0.93, 0.212, fModBoxAlphaVar, fModBoxBetaVar), "C");
      else leg->SetHeader(Form("Nominal lifetime: %.2f ms, Varied lifetime: %.2f ms", 10.4, fLifetimeVar/1000.), "C");
      leg->AddEntry(hRat, "ROI total charge ratio");
      leg->AddEntry(grChargeRatExp, "Analytical ratio");
      leg->Draw("same");  
      TLatex latex;
      latex.SetNDC();                 // use normalized coordinates (0 → 1)
      latex.SetTextSize(0.04);        // adjust size
      latex.SetTextFont(42);          // nice standard font
      latex.DrawLatex(0.15, 0.85, "DUNE Work in Progress");
      c1->SaveAs("ROI_charge_modification.pdf");
      delete c1;
    }
    evt.put(std::move(new_wires));
  }
  DEFINE_ART_MODULE(WireModifier)
} // end namespace

