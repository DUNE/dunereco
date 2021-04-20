#include "ModelConfig.h"

#include <wordexp.h>
#include <boost/property_tree/json_parser.hpp>

std::string ModelConfig::expandPath(const std::string &path)
{
    wordexp_t p;

    int s = wordexp(path.c_str(), &p, WRDE_NOCMD | WRDE_SHOWERR);

    if ((s != 0) || (p.we_wordc != 1)) {
        wordfree(&p);
        return path;
    }

    std::string result(p.we_wordv[0]);
    wordfree(&p);

    return result;
}

ModelConfig::ModelConfig(
    const std::string &savedir,
    const std::vector<InputConfigKeys> &scalarInputKeys,
    const std::vector<InputConfigKeys> &vectorInputKeys,
    const std::vector<std::string>     &outputKeys
) : savedir(ModelConfig::expandPath(savedir)),
    loaded(false),
    scalarInputKeys(scalarInputKeys),
    vectorInputKeys(vectorInputKeys),
    outputKeys(outputKeys)
{
    scalarInputs.reserve(scalarInputKeys.size());
    vectorInputs.reserve(vectorInputKeys.size());
    outputNodes.reserve(outputKeys.size());
}

bool ModelConfig::isLoaded() const
{
    return loaded;
}

void ModelConfig::load()
{
    if (loaded) {
        return;
    }

    const std::string configPath = getConfigPath();
    try {
        parseConfig(configPath);
    }
    catch (const std::exception& e) {
        throw std::runtime_error(
            "Failed to parse config file: " + configPath + ". : " + e.what()
        );
    }

    loaded = true;
}

InputConfig ModelConfig::parseInput(
    pt::ptree &tree, const InputConfigKeys &keys
)
{
    std::string nodeName = tree.get<std::string>(keys.nodeName);
    std::vector<std::string> varNames;

    for (auto &stree : tree.get_child(keys.varNames)) {
        varNames.push_back(stree.second.get_value<std::string>());
    }

    return InputConfig{ std::move(nodeName), std::move(varNames) };
}

void ModelConfig::parseConfig(const std::string &fname)
{
    pt::ptree tree;
    pt::read_json(fname, tree);

    for (const auto &inputKeys : scalarInputKeys) {
        scalarInputs.emplace_back(parseInput(tree, inputKeys));
    }

    for (const auto &inputKeys : vectorInputKeys) {
        vectorInputs.emplace_back(parseInput(tree, inputKeys));
    }

    for (const auto &outputKey : outputKeys) {
        outputNodes.emplace_back(tree.get<std::string>(outputKey));
    }
}

std::string ModelConfig::getSavedir() const
{
    return savedir;
}

std::string ModelConfig::getConfigPath() const
{
    return savedir + "/" + CONFIG_NAME;
}

std::string ModelConfig::getModelPath() const
{
    return savedir + "/" + MODEL_NAME;
}

const std::vector<InputConfig>& ModelConfig::getScalarInputs() const
{
    return scalarInputs;
}

const std::vector<InputConfig>& ModelConfig::getVectorInputs() const
{
    return vectorInputs;
}

const std::vector<std::string>& ModelConfig::getOutputNodes() const
{
    return outputNodes;
}

