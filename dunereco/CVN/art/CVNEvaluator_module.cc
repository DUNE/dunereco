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

    cvn::CaffeNetHandler fHandler;

    /// Number of outputs fron neural net
    unsigned int fNOutput;

  };



  //.......................................................................
  CVNEvaluator::CVNEvaluator(fhicl::ParameterSet const& pset):
    fPixelMapInput (pset.get<std::string>         ("PixelMapInput")),
    fResultLabel (pset.get<std::string>         ("ResultLabel")),
    fHandler       (pset.get<fhicl::ParameterSet> ("CaffeNetHandler")),
    fNOutput       (fHandler.NOutput())
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

    art::Handle< std::vector< cvn::PixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cvn::PixelMap > > pixelmaplist;
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmapListHandle))
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    if(pixelmaplist.size()>0){
      std::pair<const float*, const float*> pairedoutput= fHandler.Predict(*pixelmaplist[0]);
      const float* output=pairedoutput.first;
      //const float* features=pairedoutput.second;

      resultCol->emplace_back(output, fNOutput);

    }

    evt.put(std::move(resultCol), fResultLabel);

  }

  DEFINE_ART_MODULE(cvn::CVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







