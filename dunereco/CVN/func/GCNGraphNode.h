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
		/// Constructor with position, feature and ground truth vectors
		GCNGraphNode(std::vector<float> position, std::vector<float> features,
			std::vector<float> groundTruth);
		/// Destructor
		~GCNGraphNode(){};
		
		/// Get the node position, features or ground truth
		const std::vector<float> GetPosition() const;
		const std::vector<float> GetFeatures() const;
		const std::vector<float> GetGroundTruth() const;

		/// Add a node position coordinate
		void AddPositionCoordinate(float pos);

		/// Add a node feature
		void AddFeature(float feature);

		/// Set true ID
		void AddGroundTruth(float truth);

		/// Get the number of features
		const unsigned int GetNumberOfFeatures() const;

		/// Get the number of position coordinates
		const unsigned int GetNumberOfCoordinates() const;

		/// Get feature - zero indexed - and returns -999. if feature doesn't exist
		const float GetFeature(const unsigned int feature) const;

	private:
		std::vector<float> fPosition;
		std::vector<float> fFeatures;
		std::vector<float> fGroundTruth;
	};

}

#endif  // CVN_GCNGRAPHNODE_H
