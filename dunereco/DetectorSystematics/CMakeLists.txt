find_package(ROOT REQUIRED)
message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")

add_subdirectory(Utility)

art_make( BASENAME_ONLY
          LIBRARY_NAME  WireModifier


          MODULE_LIBRARIES 
        
        WireModUtility
        
        # ART Core
        art::Framework_Core
        art::Framework_Principal
        art::Persistency_Common
        art::Persistency_Provenance

        # ROOT 
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        ROOT::Hist
        ROOT::Graf
        ROOT::Graf3d
        ROOT::MathCore


        # FHICL
        fhiclcpp::fhiclcpp

        # Message facility
        messagefacility::MF_MessageLogger

        # cetlib
        cetlib::cetlib
        cetlib_except::cetlib_except

        # ART Root I/O
        art_root_io::TFileService_service

        # LArSoft support libs you might need:
        lardataobj::RawData
        lardataobj::RecoBase
        lardataalg::DetectorInfo
        larcoreobj::SummaryData
        larcorealg::Geometry
        larcore::Geometry_Geometry_service

        # Pandora Interface (as you already listed)
        larpandora::LArPandoraInterface
)

install_headers()
install_fhicl()
install_source()
