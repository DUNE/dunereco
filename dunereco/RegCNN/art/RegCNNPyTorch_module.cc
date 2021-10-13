// C/C++ includes
#include <iostream>
#include <sstream>

#include <torch/script.h>
#include <torch/torch.h>

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
#include "dune/RegCNN/func/RegPixelMap3D.h"

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
            std::string fPixelMapInput;
            std::string fResultLabel;
        
            torch::jit::script::Module module;
    }; // class RegCNNPyTorch

    RegCNNPyTorch::RegCNNPyTorch(fhicl::ParameterSet const& pset):
        EDProducer(pset),
        fLibPath       (cet::getenv(pset.get<std::string>      ("LibPath", ""))),
        fNetwork       (fLibPath + "/" + pset.get<std::string> ("Network")),
        fPixelMapInput (pset.get<std::string>                  ("PixelMapInput")),
        fResultLabel   (pset.get<std::string>                  ("ResultLabel"))
    {
        produces<std::vector<cnn::RegCNNResult> >(fResultLabel);
    }

    RegCNNPyTorch::~RegCNNPyTorch() {

    }

    void RegCNNPyTorch::beginJob() {
        std::cout<<"regcnn_torch job begins ...... "<<std::endl;
        try {
            // Deserialize the ScriptModule from a file using torch::jit::load().
            module = torch::jit::load(fNetwork);
        }
        catch (const c10::Error& e) {
            std::cerr<<"error loading the model\n";
            return;
        }
        mf::LogDebug("RegCNNPyTorch::beginJob")<<"loaded model "<<fNetwork<<" ... ok\n";
    }

    void RegCNNPyTorch::endJob() {
    }

    void RegCNNPyTorch::produce(art::Event& evt) {
        /// Define containers for the things we're going to produce
        std::unique_ptr< std::vector<RegCNNResult> >
                                      resultCol(new std::vector<RegCNNResult>);

        /// Load 3D pixel map for direction reco.
        std::vector< art::Ptr< cnn::RegPixelMap3D > > pixelmap3Dlist;
	art::InputTag itag1(fPixelMapInput, fPixelMapInput);
        auto pixelmap3DListHandle = evt.getHandle< std::vector< cnn::RegPixelMap3D > >(itag1);
        if (pixelmap3DListHandle) {
            art::fill_ptr_vector(pixelmap3Dlist, pixelmap3DListHandle);
        }

        if (pixelmap3Dlist.size() > 0) {
            mf::LogDebug("RegCNNPyTorch::produce")<<"3D pixel map was made for this event, loading it as the input";

            RegPixelMap3D pm = *pixelmap3Dlist[0];

            // Convert RegPixelMap3D to at::Tensor, which is the actual input of the network
            // Currently we have two configurations for the 3D pixel map
            // 100*100*100: a pixel map centered at the vertex, need longer evaluation time
            // 32*32*32: cropped pixel map around the vertex of the interaction from 100*100*100 pm, faster
            at::Tensor t_pm;
            if (pm.IsCroppedPM()) {
                t_pm = torch::from_blob(pm.fPECropped.data(), {1,1,32,32,32});
            } else {
                t_pm = torch::from_blob(pm.fPE.data(), {1,1,100,100,100});
            }

            std::vector<torch::jit::IValue> inputs_pm;
            inputs_pm.push_back(t_pm);
            at::Tensor torchOutput = module.forward(inputs_pm).toTensor();

            // Output of the network
            // Currently only direction reconstruction utilizes 3D CNN, which has 3 output represent 3 components
            // Absolute value of 3 output are meaningless, their combination is the direction of the prong
            std::vector<float> networkOutput;
            for (unsigned int i= 0; i< 3; ++i) {
                networkOutput.push_back(torchOutput[0][i].item<float>());
                std::cout<<torchOutput[0][i].item<float>()<<std::endl;
            }

            resultCol->emplace_back(networkOutput);

            //std::cout<<output.slice(/*dim=*/1, /*start=*/0, /*end=*/10) << '\n';
            //std::cout<<output[0]<<std::endl;
        }

        evt.put(std::move(resultCol), fResultLabel);
    }

    DEFINE_ART_MODULE(cnn::RegCNNPyTorch)
} // end namespace cnn
