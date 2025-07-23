#include "ModelConfig.h"

#include <wordexp.h>
#include <boost/property_tree/json_parser.hpp>

namespace transformer{

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
    const std::string &scalarInputKey,
    const std::string &vectorInputKey,
    const std::string &outputKey
) : savedir(ModelConfig::expandPath(savedir)),
    loaded(false),
    scalarInputKey(scalarInputKey),
    vectorInputKey(vectorInputKey),
    outputKey(outputKey){}

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
        std::throw_with_nested(
            std::runtime_error(
                "Failed to parse config file: " + configPath + ". : " + e.what()
            ));
    }

    loaded = true;
}

std::vector<std::string> ModelConfig::parseInput(
    pt::ptree &tree, const std::string &key
)
{
    std::vector<std::string> varNames;

    for (auto &stree : tree.get_child(key)) {
        varNames.push_back(stree.second.get_value<std::string>());
    }

    return varNames;
}

void ModelConfig::parseConfig(const std::string &fname)
{
    pt::ptree tree;
    pt::read_json(fname, tree);

    scalarNames = parseInput(tree, scalarInputKey);
    vectorNames = parseInput(tree, vectorInputKey);
    outputNames = parseInput(tree, outputKey);
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

const std::vector<std::string>& ModelConfig::getScalarInputs() const
{
    return scalarNames;
}

const std::vector<std::string>& ModelConfig::getVectorInputs() const
{
    return vectorNames;
}

const std::vector<std::string>& ModelConfig::getOutputNodes() const
{
    return outputNames;
}

}