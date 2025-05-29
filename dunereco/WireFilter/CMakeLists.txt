#include_directories($ENV{LARDATAOBJ_INC})
#include_directories($ENV{LARCOREOBJ_INC})

cet_build_plugin(WireFilter art::module LIBRARIES
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
			            lardataobj_RecoBase
                        nug4::ParticleNavigation
                        ROOT::Core
                        art::Persistency_Common 
                        art::Persistency_Provenance 
                        art::Utilities
                        #canvas_Utilities
                        messagefacility::MF_MessageLogger
                        BASENAME_ONLY
)
target_include_directories(
    WireFilter_module
    PUBLIC
    $ENV{LARDATAOBJ_INC}
    $ENV{LARCOREOBJ_INC}
    $ENV{LARCOREALG_INC}
    $ENV{LARCORE_INC}
    $ENV{RANGE_INC}
)
target_link_directories(
    WireFilter_module
    PUBLIC
    $ENV{LARCORE_LIB}
    $ENV{LARCOREALG_LIB}
    $ENV{LARDATAOBJ_LIB}
    $ENV{CANVAS_LIB}
)

install_fhicl()