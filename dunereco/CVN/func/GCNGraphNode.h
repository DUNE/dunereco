////////////////////////////////////////////////////////////////////////
/// \file    GCNGraphNode.h
/// \brief   Node for GCN
/// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CVN_GCNGRAPHNODE_H
#define CVN_GCNGRAPHNODE_H

#include <ostream>
#include <vector>

namespace cvn
{

  class GCNGraphNode
  {
  public:

    /// Default constructor
    GCNGraphNode();
    /// Constructor with position and feature vectors
    GCNGraphNode(std::vector<float> position, std::vector<float> features);
    /// Destructor
    ~GCNGraphNode(){};
    
    /// Get the node position
    const std::vector<float> GetPosition() const;

    /// Get the node features
    const std::vector<float> GetFeatures() const;

    /// Add a node position coordinate
    void AddPositionCoordinate(float pos);

    /// Add a node feature
    void AddFeature(float feature);

    /// Get the number of features
    const unsigned int GetNumberOfFeatures() const;

    /// Get the number of position coordinates
    const unsigned int GetNumberOfCoordinates() const;

    /// Get feature - zero indexed - and returns -999. if feature doesn't exist
    const float GetFeature(const unsigned int feature) const;

  private:
    std::vector<float> fPosition;
    std::vector<float> fFeatures;
  };

}

#endif  // CVN_GCNGRAPHNODE_H
