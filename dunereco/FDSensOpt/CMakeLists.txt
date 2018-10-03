add_subdirectory(IniSegAlg)
add_subdirectory(MVAAlg)
add_subdirectory(FDSensOptData)

include_directories ( $ENV{GENIE_INC}/GENIE )

# genie
cet_find_library( GALGORITHM NAMES GAlgorithm PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GBARYONRESONANCE NAMES GBaryonResonance PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GBASE NAMES GBase PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GBODEKYANG NAMES GBodekYang PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GCHARM NAMES GCharm PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GCOH NAMES GCoh PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GDFRC NAMES GDfrc PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GDIS NAMES GDIS PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GCROSSSECTIONS NAMES GCrossSections PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GDECAY NAMES GDecay PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GELAS NAMES GElas PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GELFF NAMES GElFF PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GHEP NAMES GHEP PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GEVGCORE NAMES GEVGCore  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GEVGMODULES NAMES GEVGModules PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GEVGDRIVERS NAMES GEVGDrivers PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GGIBUU NAMES GGiBUU PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GHADRONTRANSP NAMES GHadronTransp PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFRAGMENTATION NAMES GFragmentation PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GINTERACTION NAMES GInteraction PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GLLEWELLYNSMITH NAMES GLlewellynSmith  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GMEC NAMES GMEC PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GMESSENGER NAMES GMessenger PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUGAMMA NAMES GNuGamma PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUE NAMES GNuE PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNTUPLE NAMES GNtuple PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUCLEAR NAMES GNuclear PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUMERICAL NAMES GNumerical PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GQPM NAMES GQPM PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GPDG NAMES GPDG PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GPDF NAMES GPDF PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GQEL NAMES GQEL PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GRES NAMES GRES PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GREGISTRY NAMES GRegistry PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
#cet_find_library( GREINSEGHAL NAMES GReinSeghal PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GUTILS NAMES GUtils PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GGEO NAMES GGeo PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GFLUXDRIVERS NAMES GFluxDrivers PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GMUELOSS NAMES GMuELoss PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GREWEIGHT NAMES GReWeight PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
cet_find_library( GNUCLEONDECAY NAMES GNucleonDecay PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
# prior to GENIE R-2_9_0 Sehgal's name was spelled incorrectly
# so if minor version # is 8 or less use GReinSeghal, otherwise GReinSehgal
# parse_ups_version returns VMAJ, VMIN, and VPRJ
parse_ups_version( ${GENIE_VERSION} )
set ( GRSSPELLING GReinSehgal )
if ( ${VMAJ} EQUAL 2 AND ${VMIN} LESS 9 )
  set ( GRSSPELLING GReinSeghal )
endif ( ${VMAJ} EQUAL 2 AND ${VMIN} LESS 9 )
message(STATUS "GENIE version ${GENIE_VERSION} major ${VMAJ} minor ${VMIN} ${GRSSPELLING}")
cet_find_library( GREINSEHGAL NAMES ${GRSSPELLING} PATHS ENV GENIE_LIB NO_DEFAULT_PATH )


simple_plugin(NueAna   "module"
		        lardataobj_RawData
		        lardataobj_RecoBase
		        lardata_ArtDataHelper
			larreco_Calorimetry
		        larreco_RecoAlg_PMAlg
		        IniSegAlg
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(MVASelect   "module"
                        MVAAlg
                        nutools_NuReweight
                        nutools_NuReweight_art
			nusimdata_SimulationBase
                        ${GHEP}
                        ${GEVGCORE}
			${GPDG}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			${GALGORITHM}		
			${GBARYONRESONANCE}
			${GBASE}
			${GBODEKYANG}
			${GCHARM}
			${GCOH}
			${GDFRC}
			${GDIS}
			${GCROSSSECTIONS}
			${GDECAY}
			${GELAS}
			${GELFF}
			${GHEP}
			${GEVGCORE}
			${GEVGMODULES}
			${GEVGDRIVERS}
			${GGIBUU}
			${GHADRONTRANSP}
			${GFRAGMENTATION}
			${GINTERACTION}
			${GLLEWELLYNSMITH}
			${GMEC}
			${GNUGAMMA}
			${GNUE}
			${GNTUPLE}
			${GNUCLEAR}
			${GNUMERICAL}
			${GQPM}
			${GPDG}
			${GPDF}
			${GQEL}
			${GRES}
			${GREGISTRY}
			${GREINSEHGAL}
			${GUTILS}
			${GGEO}
			${GFLUXDRIVERS}
			${GMUELOSS}
			${GREWEIGHT}
			${GNUCLEONDECAY}
			${GMESSENGER}
              BASENAME_ONLY
)

simple_plugin(IniSegReco "module"
		      lardataobj_RawData
		      lardataobj_RecoBase
		      larreco_RecoAlg
		      larreco_RecoAlg_PMAlg
		      IniSegAlg
		      larcorealg_Geometry
		      larcore_Geometry_Geometry_service
		      lardata_Utilities
		      larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
		      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
		      nusimdata_SimulationBase
      		      ${ART_FRAMEWORK_CORE}
		      ${ART_FRAMEWORK_PRINCIPAL}
		      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      		      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      		      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 		      art_Persistency_Common canvas
		      art_Persistency_Provenance canvas
		      art_Utilities canvas
		      ${MF_MESSAGELOGGER}
		      
		      cetlib cetlib_except
		      ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(ShSeg   "module"
			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			IniSegAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			
			
			larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)

simple_plugin(EnergyReco   "module"
                        MVAAlg
                        nutools_NuReweight
                        nutools_NuReweight_art
                        nusimdata_SimulationBase
                        larcore_Geometry_Geometry_service
                        larreco_Calorimetry
                        larreco_RecoAlg
                        lardataobj_RecoBase
                        larsim_MCCheater_BackTrackerService_service
                        ${GHEP}
                        ${GEVGCORE}
                        ${GPDG}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${GALGORITHM}
                        ${GBARYONRESONANCE}
                        ${GBASE}
                        ${GBODEKYANG}
                        ${GCHARM}
                        ${GCOH}
                        ${GDFRC}
                        ${GDIS}
                        ${GCROSSSECTIONS}
                        ${GDECAY}
                        ${GELAS}
                        ${GELFF}
                        ${GHEP}
                        ${GEVGCORE}
                        ${GEVGMODULES}
                        ${GEVGDRIVERS}
                        ${GGIBUU}
                        ${GHADRONTRANSP}
                        ${GFRAGMENTATION}
                        ${GINTERACTION}
                        ${GLLEWELLYNSMITH}
                        ${GMEC}
                        ${GNUGAMMA}
                        ${GNUE}
                        ${GNTUPLE}
                        ${GNUCLEAR}
                        ${GNUMERICAL}
                        ${GQPM}
                        ${GPDG}
                        ${GPDF}
                        ${GQEL}
                        ${GRES}
                        ${GREGISTRY}
                        ${GREINSEHGAL}
                        ${GUTILS}
                        ${GGEO}
                        ${GFLUXDRIVERS}
                        ${GMUELOSS}
                        ${GREWEIGHT}
                        ${GNUCLEONDECAY}
                        ${GMESSENGER}
              BASENAME_ONLY
)

simple_plugin(CheckRecoEnergy   "module"
                        nutools_NuReweight
                        nutools_NuReweight_art
                        nusimdata_SimulationBase
                        ${GHEP}
                        ${GEVGCORE}
                        ${GPDG}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ${GALGORITHM}           
                        ${GBARYONRESONANCE}
                        ${GBASE}
                        ${GBODEKYANG}
                        ${GCHARM}
                        ${GCOH}
                        ${GDFRC}
                        ${GDIS}
                        ${GCROSSSECTIONS}
                        ${GDECAY}
                        ${GELAS}
                        ${GELFF}
                        ${GHEP}
                        ${GEVGCORE}
                        ${GEVGMODULES}
                        ${GEVGDRIVERS}
                        ${GGIBUU}
                        ${GHADRONTRANSP}
                        ${GFRAGMENTATION}
                        ${GINTERACTION}
                        ${GLLEWELLYNSMITH}
                        ${GMEC}
                        ${GNUGAMMA}
                        ${GNUE}
                        ${GNTUPLE}
                        ${GNUCLEAR}
                        ${GNUMERICAL}
                        ${GQPM}
                        ${GPDG}
                        ${GPDF}
                        ${GQEL}
                        ${GRES}
                        ${GREGISTRY}
                        ${GREINSEHGAL}
                        ${GUTILS}
                        ${GGEO}
                        ${GFLUXDRIVERS}
                        ${GMUELOSS}
                        ${GREWEIGHT}
                        ${GNUCLEONDECAY}
                        ${GMESSENGER}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
