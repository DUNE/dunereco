#include "dunereco/Transformer/models/zoo/TransformerEnergyModel.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "dunereco/Transformer/art/var_extractors/DefaultInputVarExtractor.h"

namespace transformer {

class TransformerEnergyProducer : public art::EDProducer
{
public:
    explicit TransformerEnergyProducer(const fhicl::ParameterSet &pset);
    void produce(art::Event &evt) override;

private:
    DefaultInputVarExtractor inputVarExtractor;
    TransformerEnergyModel model;
    VarDict vars;
};

TransformerEnergyProducer::TransformerEnergyProducer(const fhicl::ParameterSet &pset)
  : EDProducer(pset),
    inputVarExtractor("", pset.get<fhicl::ParameterSet>("ConfigInputVars")),
    model(pset.get<std::string>("ModelPath"))
{
    produces<TransformerEnergy>();
}

void TransformerEnergyProducer::produce(art::Event &evt)
{
    inputVarExtractor.extract(evt, vars);
    TransformerEnergy energy = model.predict(vars);

    evt.put(std::make_unique<TransformerEnergy>(std::move(energy)));
}

DEFINE_ART_MODULE(TransformerEnergyProducer)

}

