# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

cet_build_plugin(SNSlicer art::module LIBRARIES
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_Simulation
                        lardataobj_RawData
			lardataobj_RecoBase
			nusimdata::SimulationBase
                        nug4::ParticleNavigation
                        lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
                        ROOT::Core
                        art::Persistency_Common canvas_Persistency_Common
                        art::Persistency_Provenance canvas_Persistency_Provenance
                        art::Utilities canvas_Utilities
                        messagefacility::MF_MessageLogger
)

cet_make( LIBRARIES lardataobj_RawData
                    canvas_Persistency_Common
                    canvas_Persistency_Provenance
                    canvas_Utilities
                    messagefacility::MF_MessageLogger
                    cetlib::cetlib
                    NO_DICTIONARY )


art_dictionary( DICTIONARY_LIBRARIES ) # dunereco_SNSlicer )

install_headers()
install_fhicl()
install_source()
install_scripts()
