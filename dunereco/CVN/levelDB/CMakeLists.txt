art_make_exec( cvnCreateZlibImages       
               SOURCE    cvnCreateZlibImages.cc
               LIBRARIES art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
                         art_root_io::tfile_support
                         art_root_io::TFileService_service
                         ROOT::Core
                         art::Persistency_Common
                         canvas::canvas
                         art::Persistency_Provenance
                         z
                         art::Utilities
                         Boost::program_options
                         fhiclcpp::fhiclcpp
                         CETLIB
                         ROOT_BASIC_LIB_LIST
                         dunereco::CVN_func
               )


install_source()
install_fhicl()
