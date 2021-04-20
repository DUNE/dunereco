#include "PFParticleVarExtractor.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"

#include "utils.h"

namespace VLN {

static const std::vector<std::string> SCALAR_VARS({ });
static const std::vector<std::string> VECTOR_VARS({
    "length",
    "is_shower",
    "start.x",
    "start.y",
    "start.z",
    "dir.x",
    "dir.y",
    "dir.z",
    "energy",
    "nHit",
    "charge",
    "calE",
});

PFParticleVarExtractor::PFParticleVarExtractor(
    const std::string    &prefix,
    calo::CalorimetryAlg &algCalorimetry,
    const std::string    &labelPFPModule,
    const std::string    &labelPFPTrack,
    const std::string    &labelPFPShower,
    unsigned int         plane
) : VarExtractorBase(prefix, SCALAR_VARS, VECTOR_VARS),
    algCalorimetry(algCalorimetry),
    labelPFPModule(labelPFPModule),
    labelPFPTrack(labelPFPTrack),
    labelPFPShower(labelPFPShower),
    plane(plane)
{ }

void PFParticleVarExtractor::extractTrackVars(
    const art::Event                  &evt,
    const art::Ptr<recob::PFParticle> &particle,
    VarDict                           &vars
)
{
    auto track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(
        particle, evt, labelPFPModule, labelPFPTrack
    );
    auto start = track->Start();
    auto dir   = track->StartDirection();

    appendToVectorVar(vars, "is_shower", 0);
    appendToVectorVar(vars, "length",    track->Length());
    appendToVectorVar(vars, "start.x",   start.x());
    appendToVectorVar(vars, "start.y",   start.y());
    appendToVectorVar(vars, "start.z",   start.z());
    appendToVectorVar(vars, "dir.x",     dir.x());
    appendToVectorVar(vars, "dir.y",     dir.y());
    appendToVectorVar(vars, "dir.z",     dir.z());
    appendToVectorVar(vars, "energy",    track->StartMomentum());

    const auto hits = dune_ana::DUNEAnaTrackUtils::GetHits(
        track, evt, labelPFPTrack
    );

    extractBasicVars(evt, hits, vars);
}

void PFParticleVarExtractor::extractShowerVars(
    const art::Event                  &evt,
    const art::Ptr<recob::PFParticle> &particle,
    VarDict                           &vars
)
{
    auto shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(
        particle, evt, labelPFPModule, labelPFPShower
    );
    auto start  = shower->ShowerStart();
    auto dir    = shower->Direction();
    auto energy = shower->Energy();

    appendToVectorVar(vars, "is_shower", 1);
    appendToVectorVar(vars, "length",    shower->Length());
    appendToVectorVar(vars, "start.x",   start.x());
    appendToVectorVar(vars, "start.y",   start.y());
    appendToVectorVar(vars, "start.z",   start.z());
    appendToVectorVar(vars, "dir.x",     dir.x());
    appendToVectorVar(vars, "dir.y",     dir.y());
    appendToVectorVar(vars, "dir.z",     dir.z());
    appendToVectorVar(
        vars, "energy",  (energy.size() < plane + 1) ? 0 : energy[plane]
    );

    const auto hits = dune_ana::DUNEAnaShowerUtils::GetHits(
        shower, evt, labelPFPShower
    );

    extractBasicVars(evt, hits, vars);
}

void PFParticleVarExtractor::extractBasicVars(
    const art::Event                        &evt,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    VarDict                                 &vars
)
{
    const auto chargeCalE = calcHitsChargeCalE(
        hits, evt, algCalorimetry, plane
    );

    appendToVectorVar(vars, "nHit",   hits.size());
    appendToVectorVar(vars, "charge", chargeCalE.first);
    appendToVectorVar(vars, "calE",   chargeCalE.second);
}

void PFParticleVarExtractor::extractVars(const art::Event &evt, VarDict &vars)
{
    const std::vector<art::Ptr<recob::PFParticle>> particles
        = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, labelPFPModule);

    for (const auto &particle : particles)
    {
        const bool isTrack = dune_ana::DUNEAnaPFParticleUtils::IsTrack(
            particle, evt, labelPFPModule, labelPFPTrack
        );

        const bool isShower = dune_ana::DUNEAnaPFParticleUtils::IsShower(
            particle, evt, labelPFPModule, labelPFPShower
        );

        if (isTrack == isShower) {
            /* TODO: Maybe handle this condition somehow? */
            continue;
        }

        if (isShower) {
            extractShowerVars(evt, particle, vars);
        }
        else if (isTrack) {
            extractTrackVars(evt, particle, vars);
        }
    }
}

}

