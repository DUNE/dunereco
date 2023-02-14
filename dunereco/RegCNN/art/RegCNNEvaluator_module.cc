////////////////////////////////////////////////////////////////////////
// \file    RegCNNEvaluator_module.cc
// \brief   Producer module creating RegCNN results modified from CVNEvaluator_module.cc
// \author  Ilsoo Seong - iseong@uci.edu
//
// Modifications to interface for numu energy estimation
// Modifications to interface for direction reconstruction (electron and muon)
//  - Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include "TFile.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVectorD.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "dunereco/RegCNN/func/RegCNNResult.h"
#include "dunereco/RegCNN/func/RegPixelMap.h"
#include "dunereco/RegCNN/art/TFRegNetHandler.h"
#include "dunereco/RegCNN/art/RegCNNVtxHandler.h"
#include "dunereco/RegCNN/art/RegCNNNumuHandler.h"

namespace cnn {

  class RegCNNEvaluator : public art::EDProducer {

    public:

      explicit RegCNNEvaluator(fhicl::ParameterSet const& pset);
      ~RegCNNEvaluator();

      void produce(art::Event& evt);
      void beginJob();
      void endJob();

    private:

      void PrepareEvent(const art::Event& event);
      bool insideContVol(const double posX, const double posY, const double posZ);

      art::ServiceHandle<geo::Geometry> fGeom;

      // Module label for input pixel maps
      std::string fPixelMapInput;
      std::string fResultLabel;
      std::string fCNNType;
      std::string fTarget;

      cnn::TFRegNetHandler fTFHandler;
      cnn::RegCNNVtxHandler fRegCNNVtxHandler;
      cnn::RegCNNNumuHandler fRegCNNNumuHandler;

      std::string fHitsModuleLabel;
      std::string fTrackModuleLabel;

      double fContVolCut;

      // for checking track status
      bool fLongestTrackContained;

      void getCM(const RegPixelMap& pm, std::vector<float> &cm_list);
  }; // class RegCNNEvaluator

  //.......................................................................
  RegCNNEvaluator::RegCNNEvaluator(fhicl::ParameterSet const& pset):
    EDProducer(pset),
    fPixelMapInput     (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel       (pset.get<std::string>         ("ResultLabel")),
    fCNNType           (pset.get<std::string>         ("CNNType")),
    fTarget            (pset.get<std::string>         ("Target")),
    fTFHandler         (pset.get<fhicl::ParameterSet> ("TFNetHandler")),
    fRegCNNVtxHandler  (pset.get<fhicl::ParameterSet> ("RegCNNVtxHandler")),
    fRegCNNNumuHandler (pset.get<fhicl::ParameterSet> ("RegCNNNumuHandler")),
    fHitsModuleLabel   (pset.get<std::string>         ("HitsModuleLabel")),
    fTrackModuleLabel  (pset.get<std::string>         ("TrackModuleLabel")),
    fContVolCut        (pset.get<double>              ("ContVolCut"))
  {
    produces< std::vector<cnn::RegCNNResult> >(fResultLabel);
  }

  //......................................................................
  RegCNNEvaluator::~RegCNNEvaluator()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void RegCNNEvaluator::beginJob()
  {  
  }

  //......................................................................
  void RegCNNEvaluator::endJob()
  {
  }

  //......................................................................
  void RegCNNEvaluator::getCM(const RegPixelMap& pm, std::vector<float> &cm_list)
  {
    //std::cout << pm.fBound.fFirstWire[0]+pm.fNWire/2 << std::endl;
    //std::cout << pm.fBound.fFirstTDC[0]+pm.fNTdc*pm.fNTRes/2 << std::endl;
    for (int ii = 0; ii < 3; ii++){
        float mean_wire = pm.fBound.fFirstWire[ii]+pm.fNWire*pm.fNWRes/2;
        float mean_tdc  = pm.fBound.fFirstTDC[ii]+pm.fNTdc*pm.fNTRes/2;
        cm_list[2*ii] = mean_tdc;
        cm_list[2*ii+1] = mean_wire;
    }
  }

  //......................................................................
  void RegCNNEvaluator::produce(art::Event& evt)
  {

    this->PrepareEvent(evt);

    /// Define containers for the things we're going to produce
    std::unique_ptr< std::vector<RegCNNResult> >
                                  resultCol(new std::vector<RegCNNResult>);

    /// Load in the pixel maps
    std::vector< art::Ptr< cnn::RegPixelMap > > pixelmaplist;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto pixelmapListHandle = evt.getHandle< std::vector< cnn::RegPixelMap > >(itag1);
    if (pixelmapListHandle){
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);
    }

    /// Load 3D pixel map for direction reco.
    std::vector< art::Ptr< cnn::RegPixelMap3D > > pixelmap3Dlist;
    art::InputTag itag2(fPixelMapInput, fPixelMapInput);
    auto pixelmap3DListHandle = evt.getHandle< std::vector< cnn::RegPixelMap3D > >(itag2);
    if (pixelmap3DListHandle) {
        art::fill_ptr_vector(pixelmap3Dlist, pixelmap3DListHandle);
    }

    /// Make sure we have a valid name for the CNN type
    if(fCNNType == "TF" || fCNNType == "Tensorflow" || fCNNType == "TensorFlow"){
        // If we have a pixel map then use the TF interface to give us a prediction
        if(pixelmaplist.size() > 0){
            std::vector<float> networkOutput;
            if (fTarget == "nueenergy"){
                networkOutput = fTFHandler.Predict(*pixelmaplist[0]);
                //std::cout << "-->" << networkOutput[0] << std::endl;
            }
            else if (fTarget == "nuevertex"){
                std::vector<float> center_of_mass(6,0);
                getCM(*pixelmaplist[0], center_of_mass);
                std::cout << "cm: " << center_of_mass[0] << " " << center_of_mass[1] << " " << center_of_mass[2] << std::endl;
                networkOutput = fTFHandler.Predict(*pixelmaplist[0], center_of_mass);
                std::cout << "cnn nuevertex : "<<networkOutput[0] << " " << networkOutput[1] << " " << networkOutput[2] << std::endl;
            }
            else if (fTarget == "nuevertex_on_img"){
                auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
                auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
                networkOutput = fRegCNNVtxHandler.GetVertex(clockData, detProp, evt, *pixelmaplist[0]);
            } 
            else if (fTarget == "numuenergy") {
                networkOutput = fRegCNNNumuHandler.Predict(*pixelmaplist[0], fLongestTrackContained);
            }
            else {
                std::cout << "Wrong Target with 2D 3-view pixel maps" << std::endl;
                abort();
            }

            // cnn::Result can now take a vector of floats and works out the number of outputs
            resultCol->emplace_back(networkOutput);
        }
    } else {
        mf::LogError("RegCNNEvaluator::produce") << "CNN Type not in the allowed list: Tensorflow, Torch" << std::endl;
        mf::LogError("RegCNNEvaluator::produce") << "Exiting without processing events" << std::endl;
        return;
    } // end fCNNType

    evt.put(std::move(resultCol), fResultLabel);
  }

  void RegCNNEvaluator::PrepareEvent(const art::Event& evt) {
      // Hits
      auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);

      // Tracks
      std::vector<art::Ptr<recob::Track> > tracklist;
      auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
      if (trackListHandle)
          art::fill_ptr_vector(tracklist, trackListHandle);

      // Associations
      art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
      art::FindManyP<recob::SpacePoint> fmhs(hitListHandle, evt, fTrackModuleLabel);

      int ntracks = tracklist.size();

      double fMaxTrackLength = -1.0;
      int iLongestTrack = -1;
      // loop over tracks to find the longest track
      for (int i = 0; i < ntracks; ++i){
          if(tracklist[i]->Length() > fMaxTrackLength){
              fMaxTrackLength = tracklist[i]->Length();
              iLongestTrack = i;
          }
      }

      fLongestTrackContained = true;
      if (iLongestTrack >= 0 && iLongestTrack <= ntracks-1) {
          if (fmth.isValid()) {
              std::vector< art::Ptr<recob::Hit> > vhit = fmth.at(iLongestTrack);
              for (size_t h = 0; h < vhit.size(); ++h) {
                  if (vhit[h]->WireID().Plane == 2) {
                      std::vector< art::Ptr<recob::SpacePoint> > spts = fmhs.at(vhit[h].key());
                      if (spts.size()) {
                          if (!insideContVol(spts[0]->XYZ()[0], spts[0]->XYZ()[1], spts[0]->XYZ()[2]))
                              fLongestTrackContained = false;
                      }
                  }
              }
          }
      } // End of search longestTrack
  }

  bool RegCNNEvaluator::insideContVol(const double posX, const double posY, const double posZ) {
      geo::Point_t const vtx{posX, posY, posZ};
      bool inside = false;

      geo::TPCID idtpc = fGeom->FindTPCAtPosition(vtx);

      if (fGeom->HasTPC(idtpc)) {
          const geo::TPCGeo& tpcgeo = fGeom->GetElement(idtpc);
          double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
          double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
          double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

          for (auto const& tpcg : fGeom->Iterate<geo::TPCGeo>()) {
                  if (tpcg.MinX() < minx) minx = tpcg.MinX();
                  if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
                  if (tpcg.MinY() < miny) miny = tpcg.MinY();
                  if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
                  if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
                  if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
          }

          //x
          double dista = fabs(minx - posX);
          double distb = fabs(posX - maxx);
          if ((posX > minx) && (posX < maxx) &&
                  (dista > fContVolCut) && (distb > fContVolCut)) inside = true;
          //y
          dista = fabs(maxy - posY);
          distb = fabs(posY - miny);
          if (inside && (posY > miny) && (posY < maxy) &&
                  (dista > fContVolCut) && (distb > fContVolCut)) inside = true;
          else inside = false;
          //z
          dista = fabs(maxz - posZ);
          distb = fabs(posZ - minz);
          if (inside && (posZ > minz) && (posZ < maxz) &&
                  (dista > fContVolCut) && (distb > fContVolCut)) inside = true;
          else inside = false;
      }

      return inside;
  }

  DEFINE_ART_MODULE(cnn::RegCNNEvaluator)
} // end namespace cnn
