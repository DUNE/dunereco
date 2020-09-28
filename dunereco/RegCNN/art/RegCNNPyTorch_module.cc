// C/C++ includes
#include <iostream>
#include <sstream>

//#include <torch/script.h>
//#include <torch/torch.h>

#include "cetlib/getenv.h"

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

#include "dune/RegCNN/func/RegCNNResult.h"

namespace cnn {

    class RegCNNPyTorch : public art::EDProducer {

        public:
            explicit RegCNNPyTorch(fhicl::ParameterSet const& pset);
            ~RegCNNPyTorch();

            void produce(art::Event& evt);
            void beginJob();
            void endJob();

        private:

            std::string fLibPath;
            std::string fNetwork;
            std::string fResultLabel;
    }; // class RegCNNPyTorch

    RegCNNPyTorch::RegCNNPyTorch(fhicl::ParameterSet const& pset):
        EDProducer(pset),
        fLibPath    (cet::getenv(pset.get<std::string>("LibPath", ""))),
        fNetwork    (fLibPath + "/" + pset.get<std::string> ("Network")),
        fResultLabel(pset.get<std::string> ("ResultLabel"))
    {
        produces<std::vector<cnn::RegCNNResult> >(fResultLabel);
    }

    RegCNNPyTorch::~RegCNNPyTorch() {

    }

    void RegCNNPyTorch::beginJob() {
        std::cout<<"job begins ...... "<<std::endl;
        std::cout<<"loading model "<<fNetwork<<" ......"<<std::endl;
    }

    void RegCNNPyTorch::endJob() {
    }

    void RegCNNPyTorch::produce(art::Event& evt) {
        std::cout<<"hello"<<std::endl;
        /// Define containers for the things we're going to produce
        std::unique_ptr< std::vector<RegCNNResult> >
                                      resultCol(new std::vector<RegCNNResult>);

        //torch::jit::script::Module model{ torch::jit::load(fNetwork) };

        std::vector<float> networkOutput(3);
        resultCol->emplace_back(networkOutput);
        evt.put(std::move(resultCol), fResultLabel);
    }

    DEFINE_ART_MODULE(cnn::RegCNNPyTorch)
} // end namespace cnn
