#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "dune/VLNets/art/var_extractors/DefaultInputVarExtractor.h"
#include "dune/VLNets/models/zoo/VLNEnergyModel.h"

namespace VLN {

class VLNEnergyProducer : public art::EDProducer
{
public:
    explicit VLNEnergyProducer(const fhicl::ParameterSet &pset);
    void produce(art::Event &evt) override;

private:
    DefaultInputVarExtractor inputVarExtractor;
    VLNEnergyModel model;
    VarDict vars;
};

VLNEnergyProducer::VLNEnergyProducer(const fhicl::ParameterSet &pset)
  : EDProducer(pset),
    inputVarExtractor("", pset.get<fhicl::ParameterSet>("ConfigInputVars")),
    model(pset.get<std::string>("ModelPath"))
{
    produces<VLNEnergy>();
}

void VLNEnergyProducer::produce(art::Event &evt)
{
    inputVarExtractor.extract(evt, vars);
    VLNEnergy energy = model.predict(vars);

    evt.put(std::make_unique<VLNEnergy>(std::move(energy)));
}

DEFINE_ART_MODULE(VLNEnergyProducer)

}

