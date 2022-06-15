find_ups_product(libtorch)
find_package(Torch REQUIRED)

# using SYSTEM to prevent multiple "error: extra ‘;’ [-Werror=pedantic]" errors 
include_directories(SYSTEM ${TORCH_INCLUDE_DIRS})

art_make(
        MODULE_LIBRARIES

                lardataobj_RecoBase
                lardataobj_RawData
                larcorealg_Geometry
                lardataalg_DetectorInfo
                dunepdlegacy_Overlays

                art::Persistency_Common canvas
                art::Persistency_Provenance canvas
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                ART_FRAMEWORK_CORE
                ART_FRAMEWORK_PRINCIPAL
                ART_FRAMEWORK_SERVICES_REGISTRY
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                CLHEP

                torch
                torch_cpu
                c10

                ROOT_GEOM
                ROOT_XMLIO
                ROOT_GDML
                ROOT_BASIC_LIB_LIST
                ROOT_SPECTRUM
                ROOT_ROOFIT
                ROOT_ROOFITCORE
        )

install_headers()
install_fhicl()
install_source()
