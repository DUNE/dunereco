////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       CTPGraph
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//              S.Alonso-Monsalve,              from DUNE, CERN, Aug. 2018
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "CTPGraph.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

#include "tensorflow/core/public/session_options.h"

// -------------------------------------------------------------------
tf::CTPGraph::CTPGraph(const char* graph_file_name, const std::vector<std::string> & outputs, bool & success, int ninputs, int noutputs)
{

//    std::cout << "Starting to build the graph" << std::endl;
    success = false; // until all is done correctly

    n_inputs = ninputs;
    n_outputs = noutputs;

    // Force tf to only use a single core so it doesn't eat batch farms
    tensorflow::SessionOptions options;
    tensorflow::ConfigProto &config = options.config;
    config.set_inter_op_parallelism_threads(1);
    config.set_intra_op_parallelism_threads(1);
    config.set_use_per_session_threads(false);

//    std::cout << "Starting tf session" << std::endl;
    auto status = tensorflow::NewSession(options, &fSession);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

//    std::cout << "Session started... reading network architecture" << std::endl;
    tensorflow::GraphDef graph_def;
    status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), graph_file_name, &graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

//    std::cout << "Extracting input names" << std::endl;
    size_t ng = graph_def.node().size();

    // fill input names (TODO: generic)
    for(int i=0; i<n_inputs; ++i)
    {
        fInputNames.push_back(graph_def.node()[i].name());
    }

//    std::cout << "Extracting output names" << std::endl;
    // last node as output if no specific name provided
    if (outputs.empty()) 
    {
        for(int i=n_outputs; i>0; --i)
        {
            fOutputNames.push_back(graph_def.node()[ng - i].name()); 
        }

    }
    else // or last nodes with names containing provided strings
    {
        std::string last, current, basename, name;
        for (size_t n = 0; n < ng; ++n)
        {
            name = graph_def.node()[n].name();
            auto pos = name.find("/");
            if (pos != std::string::npos) { basename = name.substr(0, pos); }
            else { continue; }

            bool found = false;
            for (const auto & s : outputs)
            {
                if (name.find(s) != std::string::npos) { found = true; break; }
            }
            if (found)
            {
                if (!last.empty() && (basename != current))
                {
                    fOutputNames.push_back(last);
                }
                current = basename;
                last = name;
            }
        }
        if (!last.empty()) { fOutputNames.push_back(last); }
    }
    if (fOutputNames.empty())
    {
        std::cout << "Output nodes not found in the graph." << std::endl;
        return;
    }

//    std::cout << "About to create graph" << std::endl;

    status = fSession->Create(graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    success = true; // ok, graph loaded from the file

//    std::cout << "Graph success? " << success << std::endl;
}

tf::CTPGraph::~CTPGraph()
{
    auto status = fSession->Close();
    if (!status.ok()) {
      std::cout << "tf::CTPGraph::dtor: " << "Close failed." << std::endl;
    }
    delete fSession;
}

// -------------------------------------------------------------------

std::vector< std::vector< std::vector<float> > > tf::CTPGraph::run(
	const std::vector<  std::vector< std::vector<float> > >  & input)
{
  // Number of objects to classify
  const unsigned int nSamples = input.size();

  // There are two inputs to our network...
  // 1) The 100 element dE/dx array
  const unsigned int nEls = input.front().at(0).size();
  tensorflow::Tensor dEdxTensor(tensorflow::DT_FLOAT,tensorflow::TensorShape({nSamples,nEls,1}));

  // 2) The 7 additional classification variables
  const unsigned int nVars = input.front().at(1).size();
  // NB: this input doesn't need a defined depth as it doesn't use convolutions
  tensorflow::Tensor varsTensor(tensorflow::DT_FLOAT,tensorflow::TensorShape({nSamples,nVars}));

//  std::cout << "Input shapes: " << input.size() << ", " << input.front().size() << ", " << input.front().at(0).size() << ", " << input.front().at(1).size() << std::endl;

  // Fill the tensors
  auto dedxInputMap = dEdxTensor.tensor<float,3>();
  auto varsInputMap = varsTensor.tensor<float,2>();
  for(unsigned int s = 0; s < nSamples; ++s){
    // dEdx first
    for(unsigned int e = 0; e < nEls; ++e){
//      std::cout << "Adding element " << e << " with value " << input.at(s).at(0).at(e) << std::endl;
      dedxInputMap(s,e,0) = input.at(s).at(0).at(e);
    }
    // And the other variables
    for(unsigned int v = 0; v < nVars; ++v){
//      std::cout << "Adding variable " << v << " with value " << input.at(s).at(1).at(v) << std::endl;
      varsInputMap(s,v) = input.at(s).at(1).at(v);
    }
  }

  std::vector<tensorflow::Tensor> inputTensors;
  inputTensors.push_back(dEdxTensor);
  inputTensors.push_back(varsTensor);

//  std::cout << "Input tensors arranged inside the interface" << std::endl;

  return run(inputTensors);
}

// -------------------------------------------------------------------

std::vector< std::vector< std::vector< float > > > tf::CTPGraph::run(const std::vector< tensorflow::Tensor > & x)
{
    // Pair up the inputs with their names in the network
    std::vector< std::pair<std::string, tensorflow::Tensor> > inputs;
    for(int i=0; i<n_inputs; ++i){
//        std::cout << "Pairing up input with name " << fInputNames[i] << " with input number " << i << std::endl;
        inputs.push_back({fInputNames[i], x[i]});
    }

    // The output from TF has dimensions nOutputs, nSamples, nNodes
    std::vector<tensorflow::Tensor> outputs;
    auto status = fSession->Run(inputs, fOutputNames, {}, &outputs);

//    std::cout << "Sorting out the outputs inside the interface" << std::endl;

    // Dimensions we want to return are  nSamples, nOutputs, nNodes
    std::vector< std::vector< std::vector<float> > > result;
    if(status.ok()){
      // The first dimension of the output vector is the number of outputs
      const unsigned int nOut = outputs.size();
      unsigned int nSamples = 0;
      // Get the number of samples and check it is the same for all outputs
      for(unsigned int o = 0; o < nOut; ++o){
        if (o == 0){
          nSamples = outputs[o].dim_size(0);
        }
        else if (nSamples != outputs[o].dim_size(0))
        {
          throw std::string("TF outputs size inconsistent.");
        }
      }

      result.resize(nSamples,std::vector< std::vector<float> >(nOut));
      for(unsigned int s = 0; s < nSamples; ++s){
        for(unsigned int o = 0; o < nOut; ++o){
          // Get the 2D tensor (nSamples,nNodes) for this output
          auto output_map = outputs[o].tensor<float,2>();
          const unsigned int nNodes = outputs[o].dim_size(1);
          result[s][o].resize(nNodes);
          for(unsigned int n = 0; n < nNodes; ++n){
            result[s][o][n] = output_map(s,n);
          }
        }
      }  
    }
    else{
      std::cout << "Processing error in Tensorflow. Returning empty output." << std::endl;
      std::cout << "Error = " << status.ToString() << std::endl;
    }

    return result;

}
// -------------------------------------------------------------------

