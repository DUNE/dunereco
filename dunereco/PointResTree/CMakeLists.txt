
simple_plugin(PointResTree "module"
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		larsim::Simulation
		lardataobj::RawData
		lardataobj::RecoBase
		lardataalg::DetectorInfo
		larreco::Calorimetry
		larsim::MCCheater_BackTrackerService_service
		larsim::MCCheater_ParticleInventoryService_service
		nusimdata::SimulationBase
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		art_root_io::tfile_support
		ROOT::Tree
		dunereco_AnaUtils

)

install_headers()
install_fhicl()
install_source()
 
