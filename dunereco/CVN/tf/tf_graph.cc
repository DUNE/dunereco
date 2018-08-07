////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       Graph
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//              S.Alonso-Monsalve,              from DUNE, CERN, Aug. 2018
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tf_graph.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

// -------------------------------------------------------------------
tf::Graph::Graph(const char* graph_file_name, const std::vector<std::string> & outputs, bool & success, int ninputs, int noutputs)
{
    success = false; // until all is done correctly

    n_inputs = ninputs;
    n_outputs = noutputs;

    auto status = tensorflow::NewSession(tensorflow::SessionOptions(), &fSession);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    tensorflow::GraphDef graph_def;
    status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), graph_file_name, &graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    size_t ng = graph_def.node().size();

    // fill input names (TODO: generic)
    for(int i=0; i<n_inputs; ++i)
    {
        fInputNames.push_back(graph_def.node()[i].name());
    }

    /*
     std::cout << "Input names: " << std::endl;
     for(int i=0; i<n_inputs; ++i)
         std::cout << fInputNames[i] << std::endl;
    */

    // last node as output if no specific name provided
    if (outputs.empty()) 
    {
        for(int i=n_outputs; i>0; --i)
        {
            fOutputNames.push_back(graph_def.node()[ng - i].name()); 
        }

        /*
        std::cout << "Output names: " << std::endl;
        for(int i=0; i<n_outputs; ++i)
            std::cout << fOutputNames[i] << std::endl;
        */
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

    status = fSession->Create(graph_def);
    if (!status.ok())
    {
        std::cout << status.ToString() << std::endl;
        return;
    }

    success = true; // ok, graph loaded from the file
}

tf::Graph::~Graph()
{
    fSession->Close();
    delete fSession;
}

// -------------------------------------------------------------------

std::vector< std::vector< std::vector<float> > > tf::Graph::run(
	const std::vector<  std::vector<  std::vector< std::vector<float> > > > & x,
	long long int samples)
{
    if ((samples == 0) || x.empty() || x.front().empty() || x.front().front().empty() || x.front().front().front().empty())
        return std::vector< std::vector< std::vector<float> > >();

    if ((samples == -1) || (samples > (long long int)x.size())) { samples = x.size(); }

    long long int
              rows = x.front().size(),
              cols = x.front().front().size(),
              depth = x.front().front().front().size();

    std::vector< tensorflow::Tensor > _x;

    // Single-output network
    if (n_inputs == 1)
    {
        _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, depth })));
        auto input_map = _x[0].tensor<float, 4>();
        for (long long int s = 0; s < samples; ++s) {
            const auto & sample = x[s];
            for (long long int r = 0; r < rows; ++r) {
                const auto & row = sample[r];
                for (long long int c = 0; c < cols; ++c) {
                    const auto & col = row[c];
                    for (long long int d = 0; d < depth; ++d) {
                        input_map(s, r, c, d) = col[d];
                    }
                }
            }
        }
    }
    // Multi-output network
    else
    {
        for(int i=0; i<depth; ++i){
            _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, 1 })));
        }

        //tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, depth }));

        for(int view=0; view<depth; ++view){
            auto input_map = _x[view].tensor<float, 4>();
            for (long long int s = 0; s < samples; ++s) {
                const auto & sample = x[s];
                for (long long int r = 0; r < rows; ++r) {
                    const auto & row = sample[r];
                    for (long long int c = 0; c < cols; ++c) {
                        const auto & col = row[c];
                        long long int d = view;
                        input_map(s, r, c, d) = col[d];
                    }
                }
            }
        }
    }

    return run(_x);
}

// -------------------------------------------------------------------

std::vector< std::vector< std::vector< float > > > tf::Graph::run(const std::vector< tensorflow::Tensor > & x)
{
    std::vector< std::pair<std::string, tensorflow::Tensor> > inputs;
    for(int i=0; i<n_inputs; ++i){
        inputs.push_back({fInputNames[i], x[i]});
    }

    /*
    // print input/outputs
    for(int i = 0; i<n_inputs; ++i)
        std::cout << inputs[i].first << std::endl;
    for(int i = 0; i<n_outputs; ++i)
        std::cout << fOutputNames[i] << std::endl;
    */
    //std::cout << "run session" << std::endl;

    std::vector<tensorflow::Tensor> outputs;
    auto status = fSession->Run(inputs, fOutputNames, {}, &outputs);

    //std::cout << "out size " << outputs.size() << std::endl;

    if (status.ok())
    {
        size_t samples = 0;

        for (size_t o = 0; o < outputs.size(); ++o)
        {
            if (o == 0) { samples = outputs[o].dim_size(0); }
            else if ((int)samples != outputs[o].dim_size(0))
            {
                throw std::string("TF outputs size inconsistent.");
            }
        }

        std::vector< std::vector< std::vector< float > > > result;
        result.resize(samples, std::vector< std::vector< float > >(outputs.size()));

        for (size_t s = 0; s < samples; ++s) 
        {
            for (size_t o = 0; o < outputs.size(); ++o)
            {
                size_t n = outputs[o].dim_size(1);                
                auto output_map = outputs[o].tensor<float, 2>();

                result[s][o].resize(outputs[o].dim_size(1));

                std::vector< float > & vs = result[s][o];
                for (size_t i = 0; i < n; ++i) 
                {
                    vs[i] = output_map(s, i);
                }            
            }
        }

        return result;
    }
    else
    {
        std::cout << status.ToString() << std::endl;
        return std::vector< std::vector< std::vector< float > > >();
    }
}
// -------------------------------------------------------------------

