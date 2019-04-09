////////////////////////////////////////////////////////////////////////
/// \file    GCNGraph.cxx
/// \brief   GCNGraph for GCN
/// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNGraphNode.h"

namespace cvn
{

  GCNGraph::GCNGraph()
  {}

  GCNGraph::GCNGraph(std::vector<GCNGraphNode> nodes):
  fNodes(nodes)
  {

  }

  GCNGraph::GCNGraph(std::vector<std::vector<float>> positions,std::vector<std::vector<float>> features)
  {
    if(positions.size() != features.size()){
      std::cerr << "The number of nodes must be the same for the position and feature vectors" << std::endl;
      assert(0);
    }
    for(unsigned int n = 0; n < positions.size(); ++n){
      this->AddNode(positions.at(n),features.at(n));
    }
  }

  // Add a new node
  void GCNGraph::AddNode(std::vector<float> position, std::vector<float> features){
    GCNGraphNode newNode(position,features);
    fNodes.push_back(newNode);
  }

  // Get the number of nodes
  const unsigned int GCNGraph::GetNumberOfNodes() const{
    return fNodes.size();
  }

  // Access nodes
  const GCNGraphNode GCNGraph::GetNode(const unsigned int index) const{

    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNode(): Can't access node with index " << index << std::endl;
      assert(0);
    }

    return fNodes.at(index);
  }


  // Return minimum and maximum coordinate values ((xmin,xmax),(ymin,ymax),(zmin,zmax))
  const std::vector<std::pair<float,float>> GCNGraph::GetMinMaxPositions() const{

    std::pair<float,float> dummyPair = std::make_pair(1.e6,-1.e6);
    std::vector<std::pair<float,float>> minMaxVals;

    for(unsigned int i = 0; i < 3; ++i) minMaxVals.push_back(dummyPair);

    for(GCNGraphNode node : fNodes){
      std::vector<float> nodePos = node.GetPosition();
      if(nodePos[0] < minMaxVals[0].first) minMaxVals[0].first = nodePos[0];
      if(nodePos[1] < minMaxVals[1].first) minMaxVals[1].first = nodePos[1];
      if(nodePos[2] < minMaxVals[2].first) minMaxVals[2].first = nodePos[2];
      if(nodePos[0] > minMaxVals[0].second) minMaxVals[0].second = nodePos[0];
      if(nodePos[1] > minMaxVals[1].second) minMaxVals[1].second = nodePos[1];
      if(nodePos[2] > minMaxVals[2].second) minMaxVals[2].second = nodePos[2];
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


  // This function returns a vector of the format for a graph 
  // with N nodes and M features per node
  // (node0_posx, node0_posy, node0_posz, node0_feature0, ... ,node0_featureM,...
  // (nodeN_posx, nodeN_posy, nodeN_posz, nodeN_feature0, ... ,nodeN_featureM)
  const std::vector<float> GCNGraph::ConvertGraphToVector() const{

    std::vector<float> nodeVector;

    for(const GCNGraphNode &node : fNodes){
      // First add the position components
      for(const float pos : node.GetPosition()){
        nodeVector.push_back(pos);
      }
      // Now add the features
      for(const float feat : node.GetFeatures()){
        nodeVector.push_back(feat);
      }
    }

    return nodeVector; 
  }

 std::ostream& operator<<(std::ostream& os, const GCNGraph& m)
  {
    os << "GCNGraph with " << m.GetNumberOfNodes() << " nodes, ";
    return os;
  }
}
