////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap3DProducer.h
/// \brief   RegPixelMap3DProducer for RegCNN modified from PixelMapProducer.h
/// \author  Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_PIXELMAP3DPRODUCER_H
#define REGCNN_PIXELMAP3DPRODUCER_H


#include <array>
#include <vector>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "dune/RegCNN/func/RegPixelMap3D.h"
#include "dune/RegCNN/func/RegCNNBoundary3D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace detinfo {
    class DetectorClocksData;
    class DetectorPropertiesData;
}

//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

namespace cnn
{
    /// Producer algorithm for RegPixelMap3D, input to CVN neural net
    class RegPixelMap3DProducer
    {
        public:
            RegPixelMap3DProducer(unsigned int nbinsX, float xmin, float xmax, 
                    unsigned int nbinsY, float ymin, float ymax,
                    unsigned int nbinsZ, float zmin, float zmax);

            /// Get boundaries for pixel map representation of cluster
            RegCNNBoundary3D Define3DBoundary(detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    const std::vector<float> &vtx);

            RegPixelMap3D Create3DMap(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    art::FindManyP<recob::SpacePoint> const& fmSPFromHits,
                    const std::vector<float> &vtx);

            RegPixelMap3D Create3DMapGivenBoundaryBySP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    const RegCNNBoundary3D& bound,
                    art::FindManyP<recob::SpacePoint> const& fmSPFromHits);

        private:

            unsigned int fNBinsX;
            unsigned int fNBinsY;
            unsigned int fNBinsZ;
            float fLengthX;
            float fLengthY;
            float fLengthZ;

            art::ServiceHandle<geo::Geometry> geom;
    };

}

#endif  // REGCNN_PIXELMAP3DPRODUCER_H
