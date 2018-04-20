////////////////////////////////////////////////////////////////////////
// \file    CVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
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
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
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

#include "larsim/MCCheater/BackTracker.h"

#include "dune/CVN/func/Result.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/art/CaffeNetHandler.h"
#include "dune/CVN/art/TFNetHandler.h"


namespace cvn {

  class CVNEvaluator : public art::EDProducer {
  public:
    explicit CVNEvaluator(fhicl::ParameterSet const& pset);
    ~CVNEvaluator();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:

    /// Module label for input pixel maps
    std::string fPixelMapInput;
    std::string fResultLabel;

    /// Can use Caffe or Tensorflow
    std::string fCVNType;

    cvn::CaffeNetHandler fCaffeHandler;
    cvn::TFNetHandler fTFHandler;

    /// Number of outputs fron neural net
    unsigned int fNOutput;

  };

  //.......................................................................
  CVNEvaluator::CVNEvaluator(fhicl::ParameterSet const& pset):
    fPixelMapInput (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel (pset.get<std::string>         ("ResultLabel")),
    fCVNType     (pset.get<std::string>         ("CVNType")),
    fCaffeHandler       (pset.get<fhicl::ParameterSet> ("CaffeNetHandler")),
    fTFHandler       (pset.get<fhicl::ParameterSet> ("TFNetHandler")),
    fNOutput       (fCaffeHandler.NOutput())
  {
    produces< std::vector<cvn::Result>   >(fResultLabel);
  }
  //......................................................................
  CVNEvaluator::~CVNEvaluator()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNEvaluator::beginJob()
  {  }

  //......................................................................
  void CVNEvaluator::endJob()
  {
  }

  //......................................................................
  void CVNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr< std::vector<Result> >
                                  resultCol(new std::vector<Result>);

    /// Load in the pixel maps
    art::Handle< std::vector< cvn::PixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cvn::PixelMap > > pixelmaplist;
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmapListHandle))
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    /// Make sure we have a valid name for the CVN type
    if(fCVNType == "Caffe"){
      if(pixelmaplist.size()>0){
        std::pair<const float*, const float*> pairedoutput= fCaffeHandler.Predict(*pixelmaplist[0]);
        const float* output=pairedoutput.first;
        //const float* features=pairedoutput.second;

        resultCol->emplace_back(output, fNOutput);

      }
    }
    else if(fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow"){
      // If we have a pixel map then use the TF interface to give us a prediction
      if(pixelmaplist.size() > 0){
        
        std::vector<float> networkOutput = fTFHandler.Predict(*pixelmaplist[0]);

        // cvn::Result can now take a vector of floats and works out the number of outputs
        resultCol->emplace_back(networkOutput);
      }
    }
    else{
      mf::LogError("CVNEvaluator::produce") << "CVN Type not in the allowed list: Tensorflow, Caffe" << std::endl;
      mf::LogError("CVNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }

    mf::LogInfo("CVNEvaluator::produce") << " Predicted: " << (*resultCol)[0].PredictedInteractionType() << std::endl; 
    std::cout << " Predicted: " << (*resultCol)[0].PredictedInteractionType() << std::endl; 

    evt.put(std::move(resultCol), fResultLabel);

  }

  DEFINE_ART_MODULE(cvn::CVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







