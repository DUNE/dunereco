////////////////////////////////////////////////////////////////////////
// \file    CVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
//          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
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

#include "dune/CVN/func/Result.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/art/CaffeNetHandler.h"
#include "dune/CVN/art/TFNetHandler.h"
#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/InteractionType.h"

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
  CVNEvaluator::CVNEvaluator(fhicl::ParameterSet const& pset):
    fPixelMapInput (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel (pset.get<std::string>         ("ResultLabel")),
    fCVNType     (pset.get<std::string>         ("CVNType")),
    fCaffeHandler       (pset.get<fhicl::ParameterSet> ("CaffeNetHandler")),
    fTFHandler       (pset.get<fhicl::ParameterSet> ("TFNetHandler")),
    fNOutput       (fCaffeHandler.NOutput()),
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
      mf::LogError("CVNEvaluator::produce") << "CVN Type not in the allowed list: Tensorflow, Caffe" << std::endl;
      mf::LogError("CVNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }
/*
// Truth level debug code
//    mf::LogInfo("CVNEvaluator::produce") << " Predicted: " << (*resultCol)[0].PredictedInteractionType() << std::endl; 

    // Leigh: temporary testing code for performance

    InteractionType interaction = kOther;

    // * monte carlo
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if(evt.getByLabel("generator",mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    //unsigned short nmclist=  mclist.size();

    //std::cout<<"mctruth: "<<nmclist<<std::endl;

    art::Ptr<simb::MCTruth> truth = mclist[0];
    simb::MCNeutrino truthN=truth->GetNeutrino();
    //truth = mclist[0];

    interaction = GetInteractionType(truthN);

    // Get the interaction vertex from the end point of the neutrino. This is 
    // because the start point of the lepton doesn't make sense for taus as they
    // are decayed by the generator and not GEANT
    TVector3 vtx = truthN.Nu().EndPosition().Vect();
    bool isFid = (fabs(vtx.X())<310 && fabs(vtx.Y())<550 && vtx.Z()>50 && vtx.Z()<1244);


//    if(isFid && pixelmaplist.size() > 0){
    
      unsigned int correctedInt = static_cast<unsigned int>(interaction);
      if(correctedInt == 13) correctedInt = 12;

    if(isFid){
      std::cout << "This fiducial event is a true " << truthN.Nu().PdgCode() << " and is it CC? " << (truthN.CCNC()==0) << " (" << correctedInt << ")" << std::endl;
      float nueProb = (*resultCol)[0].GetNueProbability();
      float numuProb = (*resultCol)[0].GetNumuProbability();
      float nutauProb = (*resultCol)[0].GetNutauProbability();
      float ncProb = (*resultCol)[0].GetNCProbability();
      std::cout << "Summed probabilities " << numuProb << ", " << nueProb << ", " << nutauProb << ", " << ncProb << std::endl;
    }

      unsigned int predInt = static_cast<unsigned int>((*resultCol)[0].PredictedInteractionType());

      ++fTotal;
//      std::cout << " Truth :: Predicted = " << correctedInt << " :: " << predInt << " (" << numuProb << ", " << nueProb << ", " << nutauProb << ", " << ncProb  << ")" << std::endl; 
      if((correctedInt >= 0 && correctedInt <=3) && (predInt >= 0 && predInt <= 3)){
        ++fCorrect;
      }
      if((correctedInt >= 4 && correctedInt <= 7) && (predInt >= 4 && predInt <= 7)){
        ++fCorrect;
      }
      if((correctedInt >= 8 && correctedInt <=11) && (predInt >= 8 && predInt <= 11)){
        ++fCorrect;
      }
      if(correctedInt == 12 && predInt == 12){
        ++fCorrect;
      }

      if(correctedInt == predInt){
        ++fFullyCorrect;
      }

      if(abs(truthN.Nu().PdgCode()) == 12 && truthN.CCNC()==0){
        ++fTotNue;
        if(truthN.Nu().E() < 2.0) ++fTotNueBins[0];
        if(truthN.Nu().E() >= 2.0 && truthN.Nu().E() <= 6.0) ++fTotNueBins[1];
        if(truthN.Nu().E() > 6.0) ++fTotNueBins[2];

        if(nueProb > 0.7){
          ++fSelNue;
          if(truthN.Nu().E() < 2.0) ++fSelNueBins[0];
          if(truthN.Nu().E() >= 2.0 && truthN.Nu().E() <= 6.0) ++fSelNueBins[1];
          if(truthN.Nu().E() > 6.0) ++fSelNueBins[2];
        }
      }
      if(abs(truthN.Nu().PdgCode()) == 14 && truthN.CCNC()==0){
        ++fTotNumu;
        if(truthN.Nu().E() < 2.0) ++fTotNumuBins[0];
        if(truthN.Nu().E() >= 2.0 && truthN.Nu().E() <= 6.0) ++fTotNumuBins[1];
        if(truthN.Nu().E() > 6.0) ++fTotNumuBins[2];

        if(numuProb > 0.5){
          ++fSelNumu;
        if(truthN.Nu().E() < 2.0) ++fSelNumuBins[0];
        if(truthN.Nu().E() >= 2.0 && truthN.Nu().E() <= 6.0) ++fSelNumuBins[1];
        if(truthN.Nu().E() > 6.0) ++fSelNumuBins[2];
        }
      } 
  
    }
*/ // End of truth level debug code

    evt.put(std::move(resultCol), fResultLabel);

  }

  DEFINE_ART_MODULE(cvn::CVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







