#pragma once

#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>

/*
 * # Model Configuration File
 *
 * Model configuration file is a json dictionary, containing names of
 * input/output nodes of a tensorflow graph. For each input node, it also
 * contains a list of variable names whose values will be used to construct
 * an input tensor for that node.
 *
 * Here is an example of the config file:
 * ```
 * {
 *      "input_node"  : "input_1",
 *      "input_vars"  : [ "var1", "var2", "var3" ],
 *      "output_node" : "target/BiasAdd"
 *  }
 * ```
 *
 * This config file describes a tensorflow graph has one input node "input_1"
 * and one output node "target/BiasAdd". The names "input_1" and
 * "target/BiasAdd" are node names in the tensorflow graph. The input tensor
 * for the "input_1" node should be constructed from values of three variables:
 * [ "var1", "var2", "var3" ].
 *
 *
 * # `ModelConfig`
 *
 * Model configuration file can be parsed with a help of `ModelConfig` class.
 * Constructor of `ModelConfig` receives names of config file keys that specify
 * input/output configuration.
 *
 * The code below can be used to parse an example configuration file above:
 * ```
 * // Keys in config file that specify "input_1" node configuration:
 * const InputConfigKeys configKeys1("input_node", "input_vars");
 *
 * const std::string savedir;   // Directory where config file is stored
 * const std::vector<InputConfigKeys> scalarInputKeys({ configKeys1 });
 * const std::vector<InputConfigKeys> vectorInputKeys({ });
 * const std::vector<std::string>     outputKeys     ({ "output_node" });
 *
 * ModelConfig config(savedir, scalarInputKeys, vectorInputKeys, outputKeys);
 * config.load(); // Parse config file
 * ```
 */

namespace pt = boost::property_tree;

const std::string CONFIG_NAME = "config.json";
const std::string MODEL_NAME  = "model.pb";

/*
 * `InputConfig` -- Input configuration for a single tensorflow graph node.
 *
 * `nodeName` -- name of the tensorflow graph node.
 * `varNames` -- list of variable names.
 *      Input tensor for `nodeName` will be constructed from variables
 *      listed in `varNames`.
 */
struct InputConfig
{
    std::string              nodeName;
    std::vector<std::string> varNames;
};

/*
 * `InputConfigKeys` -- names of the config file keys specifying input
 *      configuration for a single tensorflow graph node:
 *
 *  `nodeName` -- name of the key in a config file which value specifies
 *      name of the tensorflow graph node.
 *  `varNames` -- name of the key in a config file which value is a list
 *      of variable names.
 *
 */
struct InputConfigKeys
{
    std::string nodeName;
    std::string varNames;
};

class ModelConfig
{
public:
    ModelConfig(
        const std::string &savedir,
        const std::vector<InputConfigKeys> &scalarInputKeys,
        const std::vector<InputConfigKeys> &vectorInputKeys,
        const std::vector<std::string>     &outputKeys
    );

    // Expand environment variables in the path
    static std::string expandPath(const std::string &path);

    void load();
    bool isLoaded() const;

    std::string getConfigPath() const;
    std::string getModelPath()  const;
    std::string getSavedir()    const;

    // Get list of configurations of scalar input nodes
    const std::vector<InputConfig>& getScalarInputs() const;

    // Get list of configurations of vector input nodes
    const std::vector<InputConfig>& getVectorInputs() const;

    // Get list of output node names
    const std::vector<std::string>& getOutputNodes()  const;

private:
    std::string savedir;
    bool loaded;

    void parseConfig(const std::string &fname);

    static InputConfig parseInput(
        pt::ptree &tree, const InputConfigKeys &keys
    );

    std::vector<InputConfigKeys> scalarInputKeys;
    std::vector<InputConfigKeys> vectorInputKeys;
    std::vector<std::string>     outputKeys;

    std::vector<InputConfig> scalarInputs;
    std::vector<InputConfig> vectorInputs;
    std::vector<std::string> outputNodes;
};

