#include <stdexcept>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "dune/VLNets/art/var_extractors/DefaultInputVarExtractor.h"
#include "dune/VLNets/art/var_extractors/EventRecoEVarExtractor.h"
#include "dune/VLNets/art/var_extractors/EventMCVarExtractor.h"
#include "dune/VLNets/art/var_extractors/FiducialCutVarExtractor.h"
#include "dune/VLNets/data/exporters/CSVExporter.h"

#include "utils.h"

namespace VLN {

class VLNEnergyDataGen : public art::EDAnalyzer
{
public:
    explicit VLNEnergyDataGen(const fhicl::ParameterSet &pset);

    void analyze(const art::Event &evt) override;
    void respondToOpenInputFile(const art::FileBlock &fb) override;

protected:
    bool passesCut(const art::Event &evt, const VarDict &vars);

private:
    Flavor flavor;
    Format format;

    bool   applyFiducialCut;
    int    isCC;
    double maxEnergy;
    int    precision;

    DefaultInputVarExtractor inputVarExtractor;
    EventRecoEVarExtractor   recoEVarExtractor;
    EventMCVarExtractor      truthVarExtractor;
    FiducialCutVarExtractor  fiducialCutVarExtractor;

    VarDict vars;
    std::unique_ptr<CSVExporter> exporter;
};

VLNEnergyDataGen::VLNEnergyDataGen(const fhicl::ParameterSet &pset)
  : EDAnalyzer(pset),
    applyFiducialCut(pset.get<bool>("ApplyFiducialCut")),
    isCC(pset.get<int>("IsCC")),
    maxEnergy(pset.get<double>("MaxEnergy")),
    precision(pset.get<int>("OutputPrecision")),
    inputVarExtractor("", pset.get<fhicl::ParameterSet>("ConfigInputVars")),
    recoEVarExtractor(
        pset.get<std::string>("Flavor") + "e.",
        pset.get<std::string>("LabelRecoE")
    ),
    truthVarExtractor("mc.", pset.get<std::string>("LabelGenerator")),
    fiducialCutVarExtractor(
        "mc.",
        pset.get<fhicl::ParameterSet>("ConfigFiducialCut"),
        pset.get<std::string>("LabelGenerator")
    )
{
    flavor = parseFlavor(pset.get<std::string>("Flavor"));
    format = parseFormat(pset.get<std::string>("OutputFormat"));
}

void VLNEnergyDataGen::respondToOpenInputFile(const art::FileBlock& fb)
{
    const std::string filename = convertFilename(fb.fileName(), "./", format);

    switch (format) {
    case Format::CSV:
        exporter = std::make_unique<CSVExporter>(filename);
        exporter->setPrecision(precision);
        break;
    }
}

bool VLNEnergyDataGen::passesCut(const art::Event &evt, const VarDict &vars)
{
    if (
           (flavor != Flavor::Any)
        && (std::abs(vars.scalar.at("mc.pdg")) != static_cast<int>(flavor))
    ) {
        return false;
    }

    if ((isCC >= 0) && (vars.scalar.at("mc.isCC") != isCC)) {
        return false;
    }

    if (applyFiducialCut && (vars.scalar.at("mc.vtxContain") != 1)) {
        return false;
    }

    if ((maxEnergy > 0) && (vars.scalar.at("mc.nuE") > maxEnergy)) {
        return false;
    }

    /* TODO: find proper way to check containment for non numu events */
    if (
           (flavor == Flavor::NuMu)
        && (vars.scalar.at("numue.longestTrackContained") != 1)
    ) {
        return false;
    }

    return true;
}

void VLNEnergyDataGen::analyze(const art::Event &evt)
{
    truthVarExtractor.extract(evt, vars);
    recoEVarExtractor.extract(evt, vars);
    fiducialCutVarExtractor.extract(evt, vars);

    if (! passesCut(evt, vars)) {
        return;
    }

    inputVarExtractor.extract(evt, vars);

    exporter->exportVars(vars);
}

DEFINE_ART_MODULE(VLNEnergyDataGen)

}

