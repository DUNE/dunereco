#find_package(libtorch REQUIRED)
find_package(Torch REQUIRED)

# using SYSTEM to prevent multiple "error: extra ‘;’ [-Werror=pedantic]" errors 
include_directories(SYSTEM ${TORCH_INCLUDE_DIRS})

art_make(
        MODULE_LIBRARIES
                larcore::headers
                lardataobj::RecoBase
                lardataobj::RawData
                larcorealg::Geometry
                lardataalg::DetectorInfo
                larevt::CalibrationDBI_Providers
                lardata::headers
                dunepdlegacy::Overlays
                art_root_io::tfile_support
                art::Persistency_Common
                art::Persistency_Provenance
                art::Utilities canvas::canvas
                cetlib::cetlib cetlib_except::cetlib_except
                ROOT::Core
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                CLHEP::CLHEP
                torch
                torch_cpu
                c10
        )

install_headers()
install_fhicl()
install_source()
