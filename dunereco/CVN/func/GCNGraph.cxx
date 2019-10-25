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
    AddNode(newNode);
  }

  // Add a new node
  void GCNGraph::AddNode(std::vector<float> position, std::vector<float> features,
    std::vector<float> groundTruth){
    GCNGraphNode newNode(position,features,groundTruth);
    AddNode(newNode);
  }

  void GCNGraph::AddNode(cvn::GCNGraphNode node){
    fNodes.push_back(node);
  }

  // Get the number of nodes
  const unsigned int GCNGraph::GetNumberOfNodes() const{
    return fNodes.size();
  }

  // Access nodes
  const GCNGraphNode& GCNGraph::GetNode(const unsigned int index) const{
    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNode(): Can't access node with index " << index << std::endl;
      assert(0);
    }

    return fNodes.at(index);
  }

  GCNGraphNode& GCNGraph::GetNodeEditable(const unsigned int index){
    if(this->GetNumberOfNodes() == 0){
      std::cerr << "GCNGraph::GetNode(): Can't access node with index " << index << std::endl;
      assert(0);
    }

    return fNodes.at(index);
  }

  // Return minimum and maximum coordinate values
  const std::vector<std::pair<float,float>> GCNGraph::GetMinMaxPositions() const{

    std::pair<float,float> dummyPair = std::make_pair(1.e6,-1.e6);
    std::vector<std::pair<float,float>> minMaxVals;

    if(fNodes.size() == 0){
      std::cerr << "No nodes found in the graph, returning empty vector" << std::endl;
      return minMaxVals;
    }

    // Initialise the vector of pairs for the number of coordinates
    for(unsigned int i = 0; i < GetNumberOfNodeCoordinates(); ++i){
      minMaxVals.push_back(dummyPair);
    }

    for(GCNGraphNode node : fNodes){
      std::vector<float> nodePos = node.GetPosition();
      for(unsigned int p = 0; p < nodePos.size(); ++p){
        if(nodePos[p] < minMaxVals[p].first) minMaxVals[p].first = nodePos[p];
        if(nodePos[p] > minMaxVals[p].second) minMaxVals[p].second = nodePos[p];
      }
    }
 
    return minMaxVals;

  }

  const std::pair<float,float> GCNGraph::GetCoordinateMinMax(unsigned int coord) const{
    if(coord > fNodes.size()){
      std::cerr << "Node index is out of bounds" << std::endl;
      assert(0);
    }
    return this->GetMinMaxPositions()[coord];
  }

  // Return spacial extent of the graph in all coordinates
  const std::vector<float> GCNGraph::GetSpacialExtent() const{

    std::vector<std::pair<float,float>> minMaxVals = this->GetMinMaxPositions();

    std::vector<float> extent;
    for(std::pair<float,float> pair : minMaxVals){
      extent.push_back(pair.second - pair.first);
    }

    return extent;
  }

  const float GCNGraph::GetCoordinateSpacialExtent(unsigned int coord) const{
    if(coord > fNodes.size()){
      std::cerr << "Node index is out of bounds" << std::endl;
      assert(0);
    }
    return this->GetSpacialExtent()[coord];
  }

  // This function returns a vector of the format for a graph 
  // with N nodes, P positions and M features per node
  // (node0_pos0, node0_pos1, ... , node0_posP, node0_feature0, ... ,node0_featureM,...
  // (nodeN_pos0, nodeN_pos1, ... , nodeN_posP, nodeN_posz, nodeN_feature0, ... ,nodeN_featureM)
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
      // Now add the ground truth
      for (const float truth : node.GetGroundTruth()) {
        nodeVector.push_back(truth);
      }
    }

    return nodeVector; 
  }

  // Return the number of coordinates for each node
  const unsigned int GCNGraph::GetNumberOfNodeCoordinates() const{
    if(fNodes.size() == 0){
      std::cerr << "Graph has no nodes, returning 0" << std::endl;
      return 0;
    }
    else return fNodes[0].GetNumberOfCoordinates();
  }

  // Return the number of features for each node
  const unsigned int GCNGraph::GetNumberOfNodeFeatures() const{
    if(fNodes.size() == 0){
      std::cerr << "Graph has no nodes, returning 0" << std::endl;
      return 0;
    }
    else return fNodes[0].GetNumberOfFeatures();
  }

  std::ostream& operator<<(std::ostream& os, const GCNGraph& m)
  {
    os << "GCNGraph with " << m.GetNumberOfNodes() << " nodes, ";
    return os;
  }
}
