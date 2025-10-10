#include "LowEEvent.h"

using namespace solar;

namespace lowe {
    LowEEvent::LowEEvent()
        : fClusters(), fFlashes(), fMatchedFlash(), fIsFlashMatched(false)
    {
    }

    LowEEvent::LowEEvent(const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes)
        : fClusters(clusters), fFlashes(flashes), fMatchedFlash(), fIsFlashMatched(false)
    {
    }

    LowEEvent::LowEEvent(const LowEEvent& toCopy)
        : fClusters(toCopy.fClusters), fFlashes(toCopy.fFlashes), fMatchedFlash(toCopy.fMatchedFlash), fIsFlashMatched(toCopy.fIsFlashMatched)
    {
    }

    LowEEvent& LowEEvent::operator=(LowEEvent const& toCopy)
    {
        if (this != &toCopy) {
            fClusters = toCopy.fClusters;
            fFlashes = toCopy.fFlashes;
            fMatchedFlash = toCopy.fMatchedFlash;
            fIsFlashMatched = toCopy.fIsFlashMatched;
        }
        return *this;
    }

    void LowEEvent::initialize(const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes)
    {
        fClusters = clusters;
        fFlashes = flashes;
        fMatchedFlash = recob::OpFlash();
        fIsFlashMatched = false; // Reset matched flash status
    }
} // namespace lowe
