////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap.h
/// \brief   RegPixelMap for RegCNN modified from PixelMap.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_REGPIXELMAP_H
#define REGCNN_REGPIXELMAP_H

#include <ostream>
#include <vector>

#include "dune/RegCNN/func/RegCNNBoundary.h"
#include "dune/RegCNN/func/HitType.h"
#include "TH2F.h"

namespace cnn
{


    /// RegPixelMap, basic input to CNN neural net
    class RegPixelMap
    {
        public:
            RegPixelMap(unsigned int nWire, unsigned int nWRes, unsigned int nTdc, unsigned int nTRes, const RegCNNBoundary& bound, const bool& prongOnly);
            RegPixelMap(){};

            /// Length in wires
            unsigned int NWire() const {return fNWire;};

            /// Number of Merged wires
            unsigned int NWRes() const {return fNWRes;};

            /// Width in tdcs
            unsigned int NTdc() const {return fNTdc;};

            /// Number of Merged tdcs
            unsigned int NTRes() const {return fNTRes;};

            /// Total number of pixels in map
            unsigned int NPixel() const {return fPE.size();};

            /// Map boundary
            RegCNNBoundary Bound() const {return fBound;};

            /// Number of inputs for the neural net
            unsigned int NInput() const {return NPixel();};

            void FillInputVector(float* input) const;

            /// Add a hit to the map if it is contained within the wire, tdc rcnne
            /// Could be expanded later to add to overflow accordingly.
            void Add(const int& wire, const int& tdc,  const unsigned int& view, const double& pe, const unsigned int& tpc, int hit_prong_tag);

            void GetTPC(const int& wire, const int& tdc, const unsigned int& view,  const unsigned int& tpc);
            /// Take global wire, tdc (detector) and return index in fPE vector
            unsigned int GlobalToIndex(const int& wire,
                    const int& tdc,
                    const unsigned int& view) ;

            /// Take local wire, tdc (within map) and return index in fPE vector
            unsigned int LocalToIndex(const unsigned int& wire,
                    const unsigned int& tdc)
                const;

            /// Take global wire, tdc (detector) and return index in fPE vector
            unsigned int GlobalToIndexSingle(const int& wire,
                    const int& tdc,
                    const unsigned int& view);

            /// Draw pixel map to the screen.  This is pretty hokey and the aspect ratio
            /// is totally unrealistic.
            void Print();

            /// Return the pixel map as a 2D histogram for visualization.
            TH2F* ToTH2() const;
            TH2F* ToLabTH2() const;
            TH2F* SingleViewToTH2(const unsigned int& view) const;

            // Add Finish method in order to determine whether to produce prong only
            // pixel maps or the full event pixel maps
            void Finish();

            unsigned int      fNWire;  ///< Number of wires, length of pixel map
            unsigned int      fNWRes;  
            unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
            unsigned int      fNTRes;  
            unsigned int      fInPM;   // check Empty Pixel Map
            unsigned int      fTPC; 
            double            fdist;
            std::vector<float>   fPE;   ///< Vector of PE measurements for pixels
            std::vector<float>   fPEX;  ///< Vector of X PE measurements for pixels
            std::vector<float>   fPEY;  ///< Vector of Y PE measurements for pixels
            std::vector<float>   fPEZ;  ///< Vector of Y PE measurements for pixels
            std::vector<double>  fPur;  ///< Vector of purity for pixels
            std::vector<double>  fPurX; ///< Vector of X purity for pixels
            std::vector<double>  fPurY; ///< Vector of Y purity for pixels
            std::vector<double>  fPurZ; ///< Vector of Y purity for pixels
            std::vector<HitType> fLab;  ///< Vector of Truth labels for pixels
            std::vector<HitType> fLabX; ///< Vector of X Truth labels for pixels
            std::vector<HitType> fLabY; ///< Vector of Y Truth labels for pixels
            std::vector<HitType> fLabZ; ///< Vector of Y Truth labels for pixels

            RegCNNBoundary          fBound;    //< RegCNNBoundary of pixel map
            bool fProngOnly;                   //< whether to use prong only pixel map
            std::vector<int> fProngTagX;
            std::vector<int> fProngTagY;
            std::vector<int> fProngTagZ;

    };

    std::ostream& operator<<(std::ostream& os, const RegPixelMap& m);

}

#endif  // CNN_PIXELMAP_H
