cet_make_library(
  LIBRARY_NAME WireModUtility

  SOURCE
    WireModUtility.cc

  LIBRARIES
    # LArSoft geometry
    larcorealg::Geometry

    # Detector properties
    lardataalg::DetectorInfo

    # RecoBase data objects
    lardataobj::RecoBase

    # Simulation objects
    lardataobj::Simulation
    nusimdata::SimulationBase

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
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::Tree
)

install_headers()
install_source()
