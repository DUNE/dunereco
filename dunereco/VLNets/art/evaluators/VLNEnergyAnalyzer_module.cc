#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "dune/VLNets/art/var_extractors/DefaultInputVarExtractor.h"
#include "dune/VLNets/art/data_generators/utils.h"
#include "dune/VLNets/data/exporters/CSVExporter.h"
#include "dune/VLNets/models/zoo/VLNEnergyModel.h"

namespace VLN {

class VLNEnergyAnalyzer : public art::EDAnalyzer
{
public:
    explicit VLNEnergyAnalyzer(const fhicl::ParameterSet &pset);

    void analyze(const art::Event &evt) override;
    void respondToOpenInputFile(const art::FileBlock &fb) override;

private:
    Format format;
    int    precision;

    DefaultInputVarExtractor inputVarExtractor;
    VLNEnergyModel model;

    VarDict vars;
    std::unique_ptr<CSVExporter> exporter;
};

VLNEnergyAnalyzer::VLNEnergyAnalyzer(const fhicl::ParameterSet &pset)
  : EDAnalyzer(pset),
    precision(pset.get<int>("OutputPrecision")),
    inputVarExtractor("", pset.get<fhicl::ParameterSet>("ConfigInputVars")),
    model(pset.get<std::string>("ModelPath"))
{
    format = parseFormat(pset.get<std::string>("OutputFormat"));
}

void VLNEnergyAnalyzer::respondToOpenInputFile(const art::FileBlock& fb)
{
    const std::string filename = convertFilename(fb.fileName(), "./", format);

    switch (format) {
    case Format::CSV:
        exporter = std::make_unique<CSVExporter>(filename);
        exporter->setPrecision(precision);
        break;
    }
}

void VLNEnergyAnalyzer::analyze(const art::Event &evt)
{
    inputVarExtractor.extract(evt, vars);

    const VLNEnergy energy = model.predict(vars);

    vars.scalar["vln.energy.totalE"]     = energy.totalE;
    vars.scalar["vln.energy.primaryE"]   = energy.primaryE;
    vars.scalar["vln.energy.secondaryE"] = energy.totalE - energy.primaryE;

    exporter->exportVars(vars);
}

DEFINE_ART_MODULE(VLNEnergyAnalyzer)

}

