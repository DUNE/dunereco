////////////////////////////////////////////////////////////////////////
/// \file    GCNGraph.cxx
/// \brief   GCNGraph for GCN
/// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include "dune/CVN/func/GCNGraph.h"

namespace cvn
{

  GCNGraph::GCNGraph()
  {}

  GCNGraph::GCNGraph(std::vector<std::vector<float>> positions,std::vector<std::vector<float>> features):
  fNodePositions(positions),
  fNodeFeatures(features)
  {
    if(fNodePositions.size() != fNodeFeatures.size()){
      std::cerr << "The number of nodes must be the same for the position and feature vectors" << std::endl;
      assert(0);
    }
  }

  // Add a new node
  void GCNGraph::AddNode(std::vector<float> position, std::vector<float> features){
    fNodePositions.push_back(position);
    fNodeFeatures.push_back(features);
  }

  // Get the number of nodes
  const unsigned int GCNGraph::GetNumberOfNodes() const{
    return fNodePositions.size();
  }

  // Get number of features
  const unsigned int GCNGraph::GetNumberOfNodeFeatures() const{
    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNumberOfNodeFeatures(): Graph object has no nodes, returning 0" << std::endl;
      return 0;
    }
    else return fNodeFeatures.size();
  }

  // Access nodes
  const std::pair<const std::vector<float>,const std::vector<float>> GCNGraph::GetNode(const unsigned int index) const{

    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNode(): Can't access node with index " << index << std::endl;
      assert(0);
    }

    std::pair<const std::vector<float>,const std::vector<float>> node = std::make_pair(fNodePositions[index],fNodeFeatures[index]);
    return node;
  }

  const std::vector<float> GCNGraph::GetNodePosition(const unsigned int index) const{
    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNode(): Can't access node with index " << index << std::endl;
      assert(0);
    }

    std::vector<float> node = fNodePositions[index];
    return node;
  }

  const std::vector<float> GCNGraph::GetNodeFeatures(const unsigned int index) const{
    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNode(): Can't access node with index " << index << std::endl;
      assert(0);
    }

    std::vector<float> node = fNodeFeatures[index];
    return node;
  }

  // Return minimum and maximum coordinate values ((xmin,xmax),(ymin,ymax),(zmin,zmax))
  const std::vector<std::pair<float,float>> GCNGraph::GetMinMaxPositions() const{

    std::pair<float,float> dummyPair = std::make_pair(1.e6,-1.e6);
    std::vector<std::pair<float,float>> minMaxVals;

    for(unsigned int i = 0; i < 3; ++i) minMaxVals.push_back(dummyPair);

    for(std::vector<float> node : fNodePositions){
      if(node[0] < minMaxVals[0].first) minMaxVals[0].first = node[0];
      if(node[1] < minMaxVals[1].first) minMaxVals[1].first = node[1];
      if(node[2] < minMaxVals[2].first) minMaxVals[2].first = node[2];
      if(node[0] > minMaxVals[0].second) minMaxVals[0].second = node[0];
      if(node[1] > minMaxVals[1].second) minMaxVals[1].second = node[1];
      if(node[2] > minMaxVals[2].second) minMaxVals[2].second = node[2];
    }
 
    return minMaxVals;

  }

  const std::pair<float,float> GCNGraph::GetMinMaxX() const{
    return this->GetMinMaxPositions()[0];
  }

  const std::pair<float,float> GCNGraph::GetMinMaxY() const{
    return this->GetMinMaxPositions()[1];
  }

  const std::pair<float,float> GCNGraph::GetMinMaxZ() const{
    return this->GetMinMaxPositions()[2];
  }

  // Return spacial extent of the graph in (x,y,z)
  const std::vector<float> GCNGraph::GetSpacialExtent() const{

    std::vector<std::pair<float,float>> minMaxVals = this->GetMinMaxPositions();

    std::vector<float> extent;
    for(std::pair<float,float> pair : minMaxVals){
      extent.push_back(pair.second - pair.first);
    }

    return extent;
  }

  const float GCNGraph::GetSpacialExtentX() const{
    return this->GetSpacialExtent()[0];
  }

  const float GCNGraph::GetSpacialExtentY() const{
    return this->GetSpacialExtent()[1];
  }

  const float GCNGraph::GetSpacialExtentZ() const{
    return this->GetSpacialExtent()[2];
  }

  std::ostream& operator<<(std::ostream& os, const GCNGraph& m)
  {
    os << "GCNGraph with " << m.GetNumberOfNodes() << " nodes, ";
    return os;
  }
}
