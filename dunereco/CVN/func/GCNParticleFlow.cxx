////////////////////////////////////////////////////////////////////////
/// \file    GCNParticleFlow.cxx
/// \brief   MC truth particle flow for graph connection ground truth
/// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "dune/CVN/func/GCNParticleFlow.h"

namespace cvn
{
  /// Add true particle to map
  void GCNParticleFlow::AddParticle(unsigned int particle,
    unsigned int parent) {
    
    // Check this particle isn't already filled
    if (fTruthMap.count(particle)) {
      // If it is filled, just check the results are consistent
      if (fTruthMap[particle] != parent) {
        std::ostringstream err;
        err << "Asked to add parent " << parent << " for particle "
          << particle << ", but its parent was already assigned as "
          << fTruthMap[particle] << "!";
        throw std::runtime_error(err.str());
      }
      else return; // If the map is already filled, don't do anything
    }

    // Add particle parent
    fTruthMap[particle] = parent;
    return;

  } // function GCNParticleFlow::AddParticle

  unsigned int GCNParticleFlow::GetParent(unsigned int particle) {

    // Check corresponding map entry is filled
    if (!fTruthMap.count(particle)) {
      std::ostringstream err;
      err << "No parent defined for true particle " << particle << "!";
      throw std::runtime_error(err.str());
    }

    // Otherwise return
    return fTruthMap[particle];

  } // function GCNParticleFlow::GetParent

} // namespace cvn
