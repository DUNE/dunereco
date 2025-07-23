#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "dunereco/Transformer/art/var_extractors/DefaultInputVarExtractor.h"
#include "dunereco/Transformer/art/data_generators/utils.h"
#include "dunereco/Transformer/data/exporters/CSVExporter.h"
#include "dunereco/Transformer/models/zoo/TransformerEnergyModel.h"

namespace transformer {

class TransformerEnergyAnalyzer : public art::EDAnalyzer
{
public:
    explicit TransformerEnergyAnalyzer(const fhicl::ParameterSet &pset);

    void analyze(const art::Event &evt) override;
    void respondToOpenInputFile(const art::FileBlock &fb) override;

private:
    Format format;
    int    precision;

    DefaultInputVarExtractor inputVarExtractor;
    TransformerEnergyModel model;

    VarDict vars;
    std::unique_ptr<CSVExporter> exporter;
};

TransformerEnergyAnalyzer::TransformerEnergyAnalyzer(const fhicl::ParameterSet &pset)
  : EDAnalyzer(pset),
    precision(pset.get<int>("OutputPrecision")),
    inputVarExtractor("", pset.get<fhicl::ParameterSet>("ConfigInputVars")),
    model(pset.get<std::string>("ModelPath"))
{
    format = parseFormat(pset.get<std::string>("OutputFormat"));
}

void TransformerEnergyAnalyzer::respondToOpenInputFile(const art::FileBlock& fb)
{
    const std::string filename = convertFilename(fb.fileName(), "./", format);

    switch (format) {
    case Format::CSV:
        exporter = std::make_unique<CSVExporter>(filename);
        exporter->setPrecision(precision);
        break;
    }
}

void TransformerEnergyAnalyzer::analyze(const art::Event &evt)
{
    inputVarExtractor.extract(evt, vars);

    const TransformerEnergy energy = model.predict(vars);

    vars.scalar["transformer.energy.totalE"]     = energy.totalE;
    vars.scalar["transformer.energy.primaryE"]   = energy.primaryE;
    vars.scalar["transformer.energy.secondaryE"] = energy.totalE - energy.primaryE;

    exporter->exportVars(vars);
}

DEFINE_ART_MODULE(TransformerEnergyAnalyzer)

}

