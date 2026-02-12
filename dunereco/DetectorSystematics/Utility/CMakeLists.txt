cet_make_library(
  LIBRARY_NAME dunereco_WireMod_Utility

  SOURCE
    WireModUtility.cc

  LIBRARIES
    # LArSoft geometry
    larcorealg::Geometry
    larcorealg::Geometry::WireReadoutGeom

    # Detector properties
    lardataalg::DetectorInfo::DetectorPropertiesData

    # RecoBase data objects
    lardataobj::RecoBase::Track
    lardataobj::RecoBase::TrackHitMeta
    lardataobj::RecoBase::Hit
    lardataobj::RecoBase::Wire

    # Simulation objects
    lardataobj::Simulation::SimEnergyDeposit
    nusimdata::SimulationBase::MCParticle
    nusimdata::SimulationBase::MCTruth

    # Utilities
    lardata::Utilities

    # Message logging
    messagefacility::MF_MessageLogger

    # FHiCL
    fhiclcpp::fhiclcpp

    # ROOT libraries
    ROOT::Core
    ROOT::Hist
    ROOT::Graf
    ROOT::Graf2d
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::Tree
  BASENAME_ONLY
)

install_headers()
install_source()
