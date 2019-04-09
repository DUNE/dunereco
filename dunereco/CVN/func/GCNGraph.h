////////////////////////////////////////////////////////////////////////
/// \file    GCNGraph.h
/// \brief   GCNGraph for GCN
/// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CVN_GCNGRAPH_H
#define CVN_GCNGRAPH_H

#include <ostream>
#include <vector>
#include "dune/CVN/func/GCNGraphNode.h"

namespace cvn
{

  /// GCNGraph, basic input for the GCN
  class GCNGraph
  {
  public:

    /// Default constructor
    GCNGraph();
    /// Constructor with position and feature vectors
    GCNGraph(std::vector<std::vector<float>> positions, std::vector<std::vector<float>> features);
    GCNGraph(std::vector<GCNGraphNode> nodes);
    /// Destructor
    ~GCNGraph(){};

    /// Add a new node
    void AddNode(std::vector<float> position, std::vector<float> features);
    void AddNode(GCNGraphNode node);
    
    /// Get the number of nodes
    const unsigned int GetNumberOfNodes() const;

    /// Access nodes
    const GCNGraphNode GetNode(const unsigned int index) const;

    /// Return minimum and maximum coordinate values ((xmin,xmax),(ymin,ymax),(zmin,zmax))
    const std::vector<std::pair<float,float>> GetMinMaxPositions() const;
    const std::pair<float,float> GetMinMaxX() const;
    const std::pair<float,float> GetMinMaxY() const;
    const std::pair<float,float> GetMinMaxZ() const;

    /// Get the extent in each dimension (x,y,z)
    const std::vector<float> GetSpacialExtent() const;
    const float GetSpacialExtentX() const;
    const float GetSpacialExtentY() const;
    const float GetSpacialExtentZ() const;

    /// Function to linearise the graph to a vector for zlib file creation
    const std::vector<float> ConvertGraphToVector() const;

  private:

    /// Store the nodes
    std::vector<GCNGraphNode> fNodes;

  };

  std::ostream& operator<<(std::ostream& os, const GCNGraph& m);

}

#endif  // CVN_GCNGRAPH_H
