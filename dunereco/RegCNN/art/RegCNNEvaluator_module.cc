////////////////////////////////////////////////////////////////////////
// \file    RegCNNEvaluator_module.cc
// \brief   Producer module creating RegCNN neural net results modified from CVNEvaluator_module.cc
// \author  Ilsoo Seong - iseong@uci.edu
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

// NOvASoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"


#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/RegCNN/func/RegCNNResult.h"
#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/art/TFRegNetHandler.h"
#include "dune/RegCNN/art/RegCNNVtxHandler.h"

namespace cnn {

  class RegCNNEvaluator : public art::EDProducer {
  public:
    explicit RegCNNEvaluator(fhicl::ParameterSet const& pset);
    ~RegCNNEvaluator();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:

    /// Module label for input pixel maps
    std::string fPixelMapInput;
    std::string fResultLabel;

    std::string fCNNType;
    std::string fTarget;

    cnn::TFRegNetHandler fTFHandler;
    cnn::RegCNNVtxHandler fRegCNNVtxHandler;

    /// Number of outputs fron neural net
    //unsigned int fNOutput;

    void getCM(const RegPixelMap& pm, std::vector<float> &cm_list);
  };

  //.......................................................................
  RegCNNEvaluator::RegCNNEvaluator(fhicl::ParameterSet const& pset):
    EDProducer(pset),
    fPixelMapInput    (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel      (pset.get<std::string>         ("ResultLabel")),
    fCNNType          (pset.get<std::string>         ("CNNType")),
    fTarget           (pset.get<std::string>         ("Target")),
    fTFHandler        (pset.get<fhicl::ParameterSet> ("TFNetHandler")),
    fRegCNNVtxHandler (pset.get<fhicl::ParameterSet> ("RegCNNVtxHandler"))
  {
    produces< std::vector<cnn::RegCNNResult>   >(fResultLabel);

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
  {  }

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
        float mean_wire = pm.fBound.fFirstWire[ii]+pm.fNWire/2;
        float mean_tdc  = pm.fBound.fFirstTDC[ii]+pm.fNTdc*pm.fNTRes/2;
        cm_list[2*ii] = mean_tdc;
        cm_list[2*ii+1] = mean_wire;
    }
  }

  //......................................................................
  void RegCNNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr< std::vector<RegCNNResult> >
                                  resultCol(new std::vector<RegCNNResult>);

    /// Load in the pixel maps
    art::Handle< std::vector< cnn::RegPixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cnn::RegPixelMap > > pixelmaplist;
    //if (evt.getByLabel(fPixelMapInput, "regcnnmap", pixelmapListHandle))
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmapListHandle)){
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);
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
            //std::cout << "cm: " << center_of_mass[0] << " " << center_of_mass[1] << " " << center_of_mass[2] << std::endl;
            networkOutput = fTFHandler.Predict(*pixelmaplist[0], center_of_mass);
            //std::cout << networkOutput[0] << " " << networkOutput[1] << " " << networkOutput[2] << std::endl;
	} else if (fTarget == "nuevertex_on_img"){
	    networkOutput = fRegCNNVtxHandler.GetVertex(evt, *pixelmaplist[0]);
        } else {
            std::cout << "Wrong Target" << std::endl;
            abort();
        }

        // cnn::Result can now take a vector of floats and works out the number of outputs
        resultCol->emplace_back(networkOutput);
      }
    }
    else{
      mf::LogError("RegCNNEvaluator::produce") << "CNN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("RegCNNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }


    evt.put(std::move(resultCol), fResultLabel);

  }

  DEFINE_ART_MODULE(cnn::RegCNNEvaluator)
} // end namespace cnn
////////////////////////////////////////////////////////////////////////







