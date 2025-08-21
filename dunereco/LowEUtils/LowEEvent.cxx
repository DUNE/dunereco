#include "LowEEvent.h"

using namespace solar;

namespace lowe {
    LowEEvent::LowEEvent()
        : pset(), fClusters(), fFlashes(), fMatchedFlash(), fIsMatchedFlash(false)
    {
    }

    LowEEvent::LowEEvent(const fhicl::ParameterSet& p, const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes)
        : pset(p), fClusters(clusters), fFlashes(flashes), fMatchedFlash(), fIsMatchedFlash(false)
    {
    }

    LowEEvent::LowEEvent(const LowEEvent& toCopy)
        : pset(toCopy.pset), fClusters(toCopy.fClusters), fFlashes(toCopy.fFlashes), fMatchedFlash(toCopy.fMatchedFlash), fIsMatchedFlash(toCopy.fIsMatchedFlash)
    {
    }

    LowEEvent& LowEEvent::operator=(LowEEvent const& toCopy)
    {
        if (this != &toCopy) {
            pset = toCopy.pset;
            fClusters = toCopy.fClusters;
            fFlashes = toCopy.fFlashes;
            fMatchedFlash = toCopy.fMatchedFlash;
            fIsMatchedFlash = toCopy.fIsMatchedFlash;
        }
        return *this;
    }

    void LowEEvent::initialize(const fhicl::ParameterSet& p, const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes)
    {
        pset = fhicl::ParameterSet(p); // Copy the parameter set
        fClusters = clusters;
        fFlashes = flashes;
        fMatchedFlash = recob::OpFlash();
        fIsMatchedFlash = false; // Reset matched flash status
    }
} // namespace lowe
