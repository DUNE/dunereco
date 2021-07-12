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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

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
            RegPixelMap3DProducer(int nbinsX, int nbinsY, int nbinsZ,
                    double XResolution, double YResolution, double ZResolution,
                    bool Cropped, bool ProngOnly);

            /// Get boundaries for pixel map representation of cluster
            RegCNNBoundary3D Define3DBoundary(detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    const std::vector<float> &vtx);

            // prong tag by track
            RegPixelMap3D Create3DMap(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    art::FindManyP<recob::SpacePoint> const& fmSPFromHits,
                    art::FindManyP<recob::Track> const& fmtrkhit,
                    const std::vector<float> &vtx);

            RegPixelMap3D Create3DMapGivenBoundaryBySP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    const RegCNNBoundary3D& bound,
                    art::FindManyP<recob::SpacePoint> const& fmSPFromHits,
                    art::FindManyP<recob::Track> const& fmtrkhit,
                    const bool& Cropped,
                    const bool& ProngOnly);

            // prong tag by shower
            RegPixelMap3D Create3DMap(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    art::FindManyP<recob::SpacePoint> const& fmSPFromHits,
                    art::FindManyP<recob::Shower> const& fmshwkhit,
                    const std::vector<float> &vtx);

            RegPixelMap3D Create3DMapGivenBoundaryBySP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    std::vector< art::Ptr< recob::Hit > > const& cluster,
                    const RegCNNBoundary3D& bound,
                    art::FindManyP<recob::SpacePoint> const& fmSPFromHits,
                    art::FindManyP<recob::Shower> const& fmshwhit,
                    const bool& Cropped,
                    const bool& ProngOnly);

        private:

            int fNBinsX;
            int fNBinsY;
            int fNBinsZ;
            double fLengthX;
            double fLengthY;
            double fLengthZ;
            bool fCropped;
            bool fProngOnly;

            art::ServiceHandle<geo::Geometry> geom;
    };

}

#endif  // REGCNN_PIXELMAP3DPRODUCER_H
