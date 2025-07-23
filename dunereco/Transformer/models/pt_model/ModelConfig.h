#pragma once

#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>

/*
 * # Model Configuration File
 *
 * Model configuration file is a json dictionary, containing
 * a list of variable names whose values will be used to construct
 * an input tensor for that node.
 *
 * Here is an example of the config file:
 * ```
 * {
 *      "scalar"  : [ "scalar1", "scalar2", "scalar3" ],
 *      "vector"  : [ "vector1", "vector2", "vector3" ]
 *      "target" : ["target1", "target2"]
 *  }
 * ```
 *
 * This config file describes a PyTorch graph that has three scalar variables:
 * [ "scalar1", "scalar2", "scalar3" ], three vector variables:
 * [ "vector1", "vector2", "vector3" ], and two output nodes, ["target1", "target2"].
 *
 *
 * # `ModelConfig`
 *
 * Model configuration file can be parsed with a help of `ModelConfig` class.
 * Constructor of `ModelConfig` receives names of config file keys that specify
 * input/output configuration.
 *
 * The code below can be used to parse an example configuration file above:
 */

namespace pt = boost::property_tree;

namespace transformer{

const std::string CONFIG_NAME = "input.json";
const std::string MODEL_NAME  = "cpu_model_norm.pt";

class ModelConfig
{
public:
    ModelConfig(
        const std::string &savedir,
        const std::string &scalarInputKey,
        const std::string &vectorInputKey,
        const std::string &outputKey
    );

    // Expand environment variables in the path
    static std::string expandPath(const std::string &path);

    void load();
    bool isLoaded() const;

    std::string getConfigPath() const;
    std::string getModelPath()  const;
    std::string getSavedir()    const;

    // Get list of configurations of scalar input nodes
    const std::vector<std::string>& getScalarInputs() const;

    // Get list of configurations of vector input nodes
    const std::vector<std::string>& getVectorInputs() const;

    // Get list of output node names
    const std::vector<std::string>& getOutputNodes()  const;

private:
    std::string savedir;
    bool loaded;

    void parseConfig(const std::string &fname);

    static std::vector<std::string> parseInput(
        pt::ptree &tree, const std::string &key
    );

    std::string scalarInputKey;
    std::string vectorInputKey;
    std::string outputKey;

    std::vector<std::string> scalarNames;
    std::vector<std::string> vectorNames;
    std::vector<std::string> outputNames;
};

}