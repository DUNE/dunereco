#ifndef PAND_SELECT_PARAMS_H
#define PAND_SELECT_PARAMS_H

namespace pandoranusel
{
  class PandoraNuSelection
  {
  public:
    double selTrackPandizzleScore;
    double selTrackMichelNHits;
    double selTrackMichelElectronMVA;
    double selTrackMichelRecoEnergyPlane2;
    double selTrackDeflecAngleSD;
    double selTrackLength;
    double selTrackEvalRatio;
    double selTrackConcentration;
    double selTrackCoreHaloRatio;
    double selTrackConicalness;
    double selTrackdEdxStart;
    double selTrackdEdxEnd;
    double selTrackdEdxEndRatio;
    double selShowerBackupPandrizzleScore;
    double selShowerEnhancedPandrizzleScore;
    double selShowerEvalRatio;
    double selShowerConcentration;
    double selShowerCoreHaloRatio;
    double selShowerConicalness;
    double selShowerdEdxBestPlane;
    double selShowerDisplacement;
    double selShowerDCA;
    double selShowerWideness;
    double selShowerEnergyDensity;
    double selShowerPathwayLengthMin;
    double selShowerMaxShowerStartPathwayScatteringAngle2D;
    double selShowerMaxNPostShowerStartHits;
    double selShowerMaxPostShowerStartScatterAngle;
    double selShowerMaxPostShowerStartNuVertexEnergyAsymmetry;
    double selShowerMaxPostShowerStartShowerStartEnergyAsymmetry;
    double selShowerMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    double selShowerMinPostShowerStartShowerStartMoliereRadius;
    double selShowerMaxPostShowerStartOpeningAngle;
    double selShowerMaxFoundHitRatio;
    double selShowerMaxInitialGapSize;
    double selShowerMinLargestProjectedGapSize;
    double selShowerNViewsWithAmbiguousHits;
    double selShowereAmbiguousHitMaxUnaccountedEnergy;
  };
}

#endif
