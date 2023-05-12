////////////////////////////////////////////////////////////////////////
// \file    CVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
//          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

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

#include "dunereco/CVN/func/Result.h"
#include "dunereco/CVN/func/PixelMap.h"
//#include "dunereco/CVN/art/CaffeNetHandler.h"
#include "dunereco/CVN/art/TFNetHandler.h"
#include "dunereco/CVN/func/AssignLabels.h"
#include "dunereco/CVN/func/InteractionType.h"

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

    //cvn::CaffeNetHandler fCaffeHandler;
    cvn::TFNetHandler fTFHandler;

    /// Number of outputs fron neural net
    //unsigned int fNOutput;

    /// If there are multiple pixel maps per event can we use them?
    bool fMultiplePMs;

    unsigned int fTotal;
    unsigned int fCorrect;
    unsigned int fFullyCorrect;

    unsigned int fTotNue;
    unsigned int fSelNue;
    unsigned int fTotNumu;
    unsigned int fSelNumu;

    // Three binned versions of above
    std::vector<unsigned int> fTotNueBins;
    std::vector<unsigned int> fSelNueBins;
    std::vector<unsigned int> fTotNumuBins;
    std::vector<unsigned int> fSelNumuBins;
  };

  //.......................................................................
  CVNEvaluator::CVNEvaluator(fhicl::ParameterSet const& pset): EDProducer{pset},
    fPixelMapInput (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel (pset.get<std::string>         ("ResultLabel")),
    fCVNType     (pset.get<std::string>         ("CVNType")),
    //fCaffeHandler       (pset.get<fhicl::ParameterSet> ("CaffeNetHandler")),
    fTFHandler       (pset.get<fhicl::ParameterSet> ("TFNetHandler")),
    //fNOutput       (fCaffeHandler.NOutput()),
    fMultiplePMs (pset.get<bool> ("MultiplePMs"))
  {
    produces< std::vector<cvn::Result>   >(fResultLabel);
    fTotal = 0;
    fCorrect = 0;
    fFullyCorrect = 0;
    fTotNue = 0;
    fSelNue = 0;
    fTotNumu = 0;
    fSelNumu = 0;

    fTotNueBins = {0,0,0};
    fSelNueBins = {0,0,0};

    fTotNumuBins = {0,0,0};
    fSelNumuBins = {0,0,0};
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
    /*
    float tot = static_cast<float>(fTotal);
    std::cout << "Total of " << fTotal << " events passed the fiducial volume cut and were classified" << std::endl;
    if(fTotal > 0)   std::cout << "Correct flavour (interaction) identification = " << 100.*(fCorrect / tot) << "% (" << 100.*(fFullyCorrect/tot) << "%)" << std::endl;
    if(fTotNue > 0)  std::cout << "Nue efficiency  = " << 100.*(fSelNue/static_cast<float>(fTotNue)) << "%" << std::endl;
    for(unsigned int i = 0; i < fTotNueBins.size(); ++i){
      if(fTotNueBins[i] > 0) std::cout << "Nue efficiency " << i << "= " << 100.*(fSelNueBins[i]/static_cast<float>(fTotNueBins[i])) << "%" << std::endl;
    }
    if(fTotNumu > 0) std::cout << "Numu efficiency = " << 100.*(fSelNumu/static_cast<float>(fTotNumu)) << "%" << std::endl;
    for(unsigned int i = 0; i < fTotNumuBins.size(); ++i){
      if(fTotNumuBins[i] > 0) std::cout << "Numu efficiency " << i << "= " << 100.*(fSelNumuBins[i]/static_cast<float>(fTotNumuBins[i])) << "%" << std::endl;
    }
    */
  }

  //......................................................................
  void CVNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr< std::vector<Result> >
                                  resultCol(new std::vector<Result>);

    /// Load in the pixel maps
    std::vector<cvn::PixelMap > pixelmaplist;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto pixelmapListHandle = evt.getHandle< std::vector< cvn::PixelMap > >(itag1);
    if (pixelmapListHandle)
      pixelmaplist = *pixelmapListHandle;

    /// Make sure we have a valid name for the CVN type
    /*
    if(fCVNType == "Caffe"){
      if(pixelmaplist.size()>0){
        std::pair<const float*, const float*> pairedoutput= fCaffeHandler.Predict(*pixelmaplist[0]);
        const float* output=pairedoutput.first;
        //const float* features=pairedoutput.second;

        resultCol->emplace_back(output, fNOutput);

      }
    }*/
    if(fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow"){
      // If we have a pixel map then use the TF interface to give us a prediction
      if(pixelmaplist.size() > 0){
        
        std::vector< std::vector<float> > networkOutput = fTFHandler.Predict(*pixelmaplist[0]);
        // cvn::Result can now take a vector of floats and works out the number of outputs
        resultCol->emplace_back(networkOutput);

        /*
        for(auto const& resaux: (*resultCol))
        {
            std::cout << "Is antineutrino: " << resaux.GetIsAntineutrinoProbability() << std::endl;
            std::cout << "CC Numu: " << resaux.GetNumuProbability() << std::endl;
            std::cout << "CC Nue: " << resaux.GetNueProbability() << std::endl;
            std::cout << "CC Nutau: " << resaux.GetNutauProbability() << std::endl;
            std::cout << "NC: " << resaux.GetNCProbability() << std::endl;
            std::cout << "CC QE: " << resaux.GetQEProbability() << std::endl;
            std::cout << "CC Res: " << resaux.GetResProbability() << std::endl;
            std::cout << "CC DIS: " << resaux.GetDISProbability() << std::endl;
            std::cout << "CC Other: " << resaux.GetOtherProbability() << std::endl;
            std::cout << "0 Protons: " << resaux.Get0protonsProbability() << std::endl;
            std::cout << "1 Protons: " << resaux.Get1protonsProbability() << std::endl;
            std::cout << "2 Protons: " << resaux.Get2protonsProbability() << std::endl;
            std::cout << ">2 Protons: " << resaux.GetNprotonsProbability() << std::endl;
            std::cout << "0 Pions: " << resaux.Get0pionsProbability() << std::endl;
            std::cout << "1 Pions: " << resaux.Get1pionsProbability() << std::endl;
            std::cout << "2 Pions: " << resaux.Get2pionsProbability() << std::endl;
            std::cout << ">3 Pions: " << resaux.GetNpionsProbability() << std::endl;
            std::cout << "0 Pizeros: " << resaux.Get0pizerosProbability() << std::endl;
            std::cout << "1 Pizeros: " << resaux.Get1pizerosProbability() << std::endl;
            std::cout << "2 Pizeros: " << resaux.Get2pizerosProbability() << std::endl;
            std::cout << ">2 Pizeros: " << resaux.GetNpizerosProbability() << std::endl;
            std::cout << "0 Neutrons: " << resaux.Get0neutronsProbability() << std::endl;
            std::cout << "1 Neutrons: " << resaux.Get1neutronsProbability() << std::endl;
            std::cout << "2 Neutrons: " << resaux.Get2neutronsProbability() << std::endl;
            std::cout << ">2 Neutrons: " << resaux.GetNneutronsProbability() << std::endl << std::endl;

            std::cout << "Is antineutrino: " << resaux.PredictedIsAntineutrino() << std::endl;  
            std::cout << "Predicted flavour: " << resaux.PredictedFlavour() << std::endl; 
            std::cout << "Predicted interaction: " << resaux.PredictedInteraction() << std::endl; 
            std::cout << "Predicted protons: " << resaux.PredictedProtons() << std::endl; 
            std::cout << "Predicted pions: " << resaux.PredictedPions() << std::endl; 
            std::cout << "Predicted pizeros: " << resaux.PredictedPizeros() << std::endl; 
            std::cout << "Predicted neutrons: " << resaux.PredictedNeutrons() << std::endl; 
        }
        */

        // Classify other pixel maps if they exist
        if(fMultiplePMs){
          for(unsigned int p = 1; p < pixelmaplist.size(); ++p){
            std::vector< std::vector<float> > output = fTFHandler.Predict(*pixelmaplist[p]);
            resultCol->emplace_back(output);
          }
        }

      }
    }
    else{
      mf::LogError("CVNEvaluator::produce") << "CVN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("CVNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }
    evt.put(std::move(resultCol), fResultLabel);

  }

  DEFINE_ART_MODULE(cvn::CVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







