#pragma once

#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "VarExtractorBase.h"

namespace VLN {

class PFParticleVarExtractor : public VarExtractorBase
{
public:
    PFParticleVarExtractor(
        const std::string    &prefix,
        calo::CalorimetryAlg &algCalorimetry,
        const std::string    &labelPFPModule,
        const std::string    &labelPFPTrack,
        const std::string    &labelPFPShower,
        unsigned int         plane = 2
    );

    ~PFParticleVarExtractor() = default;

protected:
    void extractVars(const art::Event &evt, VarDict &vars) override;

private:
    void extractShowerVars(
        const art::Event                  &evt,
        const art::Ptr<recob::PFParticle> &particle,
        VarDict                           &vars
    );

    void extractTrackVars(
        const art::Event                  &evt,
        const art::Ptr<recob::PFParticle> &particle,
        VarDict                           &vars
    );

    void extractBasicVars(
        const art::Event                        &evt,
        const std::vector<art::Ptr<recob::Hit>> &hits,
        VarDict                                 &vars
    );

private:
    calo::CalorimetryAlg &algCalorimetry;
    std::string          labelPFPModule;
    std::string          labelPFPTrack;
    std::string          labelPFPShower;
    unsigned int         plane;
};

}

