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
#include "canvas/Persistency/Common/Ptr.h"

#include "dunereco/TransformerCVN/func/TransformerCVNResult.h"
#include "dunereco/TransformerCVN/func/TransformerPixelMap.h"

namespace cnn {

    class TransformerCVNEvaluator : public art::EDProducer {

        public:
            explicit TransformerCVNEvaluator(fhicl::ParameterSet const& pset);
            ~TransformerCVNEvaluator();

            void produce(art::Event& evt);
            void beginJob();
            void endJob();

        private:

            std::string fLibPath;
            std::string fNetwork;
            std::string fPixelMapInput;
            std::string fEventPixelMapInput;
            std::string fEventResultLabel;
            std::string fProngPixelMapInput;
            std::string fProngResultLabel;
            std::vector<int> fOutputPtypes;
            int fMaxProngs;

            torch::jit::script::Module module;

            std::vector<float> PixelMapToFlatVector(const TransformerPixelMap pm); 
    }; // class TransformerCVNEvaluator

    TransformerCVNEvaluator::TransformerCVNEvaluator(fhicl::ParameterSet const& pset):
        EDProducer(pset),
        fLibPath       (cet::getenv(pset.get<std::string>        ("LibPath", ""))),
        fNetwork       (fLibPath + "/" + pset.get<std::string>   ("Network")),
        fPixelMapInput (pset.get<std::string>                    ("PixelMapInput")),
        fEventPixelMapInput (pset.get<std::string>               ("EventPixelMapInput")),
        fEventResultLabel   (pset.get<std::string>               ("EventResultLabel")),
        fProngPixelMapInput (pset.get<std::string>               ("ProngPixelMapInput")),
        fProngResultLabel      (pset.get<std::string>            ("ProngResultLabel")),
        fMaxProngs          (pset.get<int>                       ("MaxProngs"))
    {
        produces<std::vector<cnn::TransformerCVNResult> >(fEventResultLabel);
        produces<std::vector<cnn::TransformerCVNResult   > >(fProngResultLabel);

        for(auto i = 0u; i < 8; i++)
	   fOutputPtypes.push_back(i);
    }

    TransformerCVNEvaluator::~TransformerCVNEvaluator() {

    }

    void TransformerCVNEvaluator::beginJob() {
        std::cout << "Loading model " << fNetwork << std::endl;
        try {
            // Deserialize the ScriptModule from a file using torch::jit::load().
            module = torch::jit::load(fNetwork);
        }
        catch (const c10::Error& e) {
            std::cerr<<"error loading the model\n";
            abort();
            return;
        }
        mf::LogDebug("TransformerCVNEvaluator::beginJob")<<"loaded model "<<fNetwork<<" ... ok\n";
    }

    void TransformerCVNEvaluator::endJob() {
    }

    void TransformerCVNEvaluator::produce(art::Event& evt) {
        /// Define containers for the things we're going to produce
        std::unique_ptr< std::vector<TransformerCVNResult> >
                                      eventResultCol(new std::vector<TransformerCVNResult>);
        std::unique_ptr< std::vector<TransformerCVNResult> >
                                      prongResultCol(new std::vector<TransformerCVNResult>);

        /// Load event pixel map.
        std::vector< art::Ptr< cnn::TransformerPixelMap > > event_pixelmaplist;
	    art::InputTag itag1(fPixelMapInput, fEventPixelMapInput);
        auto event_pixelmapListHandle = evt.getHandle< std::vector< cnn::TransformerPixelMap > >(itag1);
        if (event_pixelmapListHandle) {
            art::fill_ptr_vector(event_pixelmaplist, event_pixelmapListHandle);
        }
        
        /// Load all track and shower pixel maps.
        std::vector< art::Ptr< cnn::TransformerPixelMap > > prong_pixelmaplist;
	    art::InputTag itag2(fPixelMapInput, fProngPixelMapInput);
        auto prong_pixelmapListHandle = evt.getHandle< std::vector< cnn::TransformerPixelMap > >(itag2);
        if (prong_pixelmapListHandle) {
            art::fill_ptr_vector(prong_pixelmaplist, prong_pixelmapListHandle);
        }

        if (event_pixelmaplist.size() > 0) {
            mf::LogDebug("TransformerCVNEvaluator::produce")<<"pixel maps were made for this event, loading it as the input";

            // Prepare network inputs
            std::vector<float> inputs;

            // Prepare event input
            TransformerPixelMap pm = *event_pixelmaplist[0];
            std::vector<float> flatvector;
            flatvector = PixelMapToFlatVector(pm);
            inputs.insert(inputs.end(), flatvector.begin(), flatvector.end());

            int valid_prongs = prong_pixelmaplist.size();
            // Prepare prong inputs
            for( int iProng = 0; iProng < valid_prongs; ++iProng ) {
                pm = *prong_pixelmaplist[iProng];
                flatvector = PixelMapToFlatVector(pm);
            }
            std::cout << valid_prongs << " prongs in event\n";

            // Convert vector inputs to at::Tensor, which is the actual input of the network
            at::Tensor t_pm = torch::from_blob(inputs.data(), {3*400*280*(valid_prongs+1)}, torch::TensorOptions().dtype(torch::kFloat32));

            std::vector<torch::jit::IValue> inputs_pm;
            inputs_pm.push_back(t_pm);
            torch::jit::IValue torchOutput = module.forward(inputs_pm);

            at::Tensor event_preds = torchOutput.toTuple()->elements()[0].toTensor();
            at::Tensor prong_preds = torchOutput.toTuple()->elements()[1].toTensor();

            std::vector<float> event_preds_vec;
            for (int i = 0; i<4; i++)
            {
                event_preds_vec.push_back(event_preds.select(0, i).item().to<float>());
            }
            eventResultCol->emplace_back(event_preds_vec);

            std::vector<float> prong_preds_vec;
            for( int iProng = 0; iProng<valid_prongs; ++iProng )
            {
                if (iProng >= fMaxProngs) break;

                for(int i = 0; i < 8; i++)
                {
                    prong_preds_vec.push_back(prong_preds.select(0, iProng).select(0, i).item().to<float>());
                }
            }
            prongResultCol->emplace_back(prong_preds_vec);
        }

        evt.put(std::move(eventResultCol), fEventResultLabel);
        evt.put(std::move(prongResultCol), fProngResultLabel);
    }

    std::vector<float> TransformerCVNEvaluator::PixelMapToFlatVector(const TransformerPixelMap pm) {
        std::vector<float> flatvector;

        flatvector.insert(flatvector.end(), pm.fPEX.begin(), pm.fPEX.end());
        flatvector.insert(flatvector.end(), pm.fPEY.begin(), pm.fPEY.end());
        flatvector.insert(flatvector.end(), pm.fPEZ.begin(), pm.fPEZ.end());

        return flatvector;
    }

    DEFINE_ART_MODULE(cnn::TransformerCVNEvaluator)
} // end namespace cnn
