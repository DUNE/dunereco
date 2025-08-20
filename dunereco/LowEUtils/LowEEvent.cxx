#include "LowEEvent.h"

using namespace solar;

namespace lowe {
    LowEEvent::LowEEvent()
        : fClusters(), fFlashes(), fMatchedFlash()
    {
    }

    LowEEvent::LowEEvent(const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes)
        : fClusters(clusters), fFlashes(flashes), fMatchedFlash()
    {
    }

    LowEEvent::LowEEvent(const LowEEvent& toCopy)
        : fClusters(toCopy.fClusters), fFlashes(toCopy.fFlashes), fMatchedFlash(toCopy.fMatchedFlash)
    {
    }

    LowEEvent& LowEEvent::operator=(LowEEvent const& toCopy)
    {
        if (this != &toCopy) {
            fClusters = toCopy.fClusters;
            fFlashes = toCopy.fFlashes;
            fMatchedFlash = toCopy.fMatchedFlash;
        }
        return *this;
    }

    void LowEEvent::initialize(const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes)
    {
        fClusters = clusters;
        fFlashes = flashes;
        fMatchedFlash = recob::OpFlash();
    }
} // namespace lowe
