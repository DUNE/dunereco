////////////////////////////////////////////////////////////////////////
/// \file    GCNParticleFlow.h
/// \brief   MC truth particle flow for graph connection ground truth
/// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef CVN_GCNPARTICLEFLOW_H
#define CVN_GCNPARTICLEFLOW_H

#include <map>

namespace cvn
{

  /// GCNParticleFlow, a map of true particles to their parents
  class GCNParticleFlow
  {
  public:

    /// Default constructor
    GCNParticleFlow(){};
    /// Default destructor
    ~GCNParticleFlow(){};

    /// Add a new true particle
    void AddParticle(unsigned int particle, unsigned int parent);

    /// Retrieve truth information
    std::map<unsigned int, unsigned int> GetMap() { return fTruthMap; };
    unsigned int GetParent(unsigned int particle);

  private:

    /// Map of true particle ID to parent particle ID
    std::map<unsigned int, unsigned int> fTruthMap;

  };

}

#endif  // CVN_GCNPARTICLEFLOW_H
