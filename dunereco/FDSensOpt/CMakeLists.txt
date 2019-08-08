add_subdirectory(IniSegAlg)
add_subdirectory(MVAAlg)
add_subdirectory(FDSensOptData)

include_directories ( $ENV{GENIE_INC}/GENIE )


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
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
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
                        nugen_NuReweight
                        nugen_NuReweight_art
			nusimdata_SimulationBase
                        ${GFWGHEP}
                        ${GFWEG}
			${GFWPARDAT}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			${GFWALG}
			${GFWPARDAT}
			${GPHDISXS}
			${GPHDISXS}
			${GPHCHMXS}
			${GPHCOHEG} ${GPHCOHXS}
			${GPHDFRCEG}
			${GPHDISEG}
			${GPHXSIG}
			${GPHDCY}
			${GPHQELXS}
			${GPHQELXS}
			${GFWGHEP}
			${GFWEG}
			${GPHCMN}
			${GFWEG}
			${GPHRESXS}
			${GPHHADTRANSP}
			${GPHHADNZ}
			${GFWINT}
			${GPHQELXS}
			${GPHMNUCEG} ${GPHMNUCXS}
			${GPHAMNGEG} ${GPHAMNGXS}
			${GPHNUELEG} ${GPHNUELXS}
			${GFWNTP}
			${GPHNUCLST}
			${GFWNUM}
			${PHDISXS}
			${GFWPARDAT}
			${GPHPDF}
			${GPHQELEG}
			${GPHRESEG} ${GPHRESXS}
			${GFWREG}
			${GPHRESXS}
			${GFWUTL}
			${GTLGEO}
			${GTLFLX}
			${GPHMEL}
			${GTLREW}
			${GPHNDCY}
			${GFWMSG}
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
		      larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
		      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
		      nusimdata_SimulationBase
      		      ${ART_FRAMEWORK_CORE}
		      ${ART_FRAMEWORK_PRINCIPAL}
		      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      		      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
      		      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
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
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			
			
			larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
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
                        nugen_NuReweight
                        nugen_NuReweight_art
                        nusimdata_SimulationBase
                        larcore_Geometry_Geometry_service
                        larreco_Calorimetry
                        larreco_RecoAlg
                        lardataobj_RecoBase
                        larsim_MCCheater_BackTrackerService_service
                        ${GFWGHEP}
                        ${GFWEG}
                        ${GFWPARDAT}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${GFWALG}
                        ${GFWPARDAT}
                        ${GPHDISXS}
                        ${GPHDISXS}
                        ${GPHCHMXS}
                        ${GPHCOHEG} ${GPHCOHXS}
                        ${GPHDFRCEG}
                        ${GPHDISEG}
                        ${GPHXSIG}
                        ${GPHDCY}
                        ${GPHQELXS}
                        ${GPHQELXS}
                        ${GFWGHEP}
                        ${GFWEG}
                        ${GPHCMN}
                        ${GFWEG}
                        ${GPHRESXS}
                        ${GPHHADTRANSP}
                        ${GPHHADNZ}
                        ${GFWINT}
                        ${GPHQELXS}
                        ${GPHMNUCEG} ${GPHMNUCXS}
                        ${GPHAMNGEG} ${GPHAMNGXS}
                        ${GPHNUELEG} ${GPHNUELXS}
                        ${GFWNTP}
                        ${GPHNUCLST}
                        ${GFWNUM}
                        ${PHDISXS}
                        ${GFWPARDAT}
                        ${GPHPDF}
                        ${GPHQELEG}
                        ${GPHRESEG} ${GPHRESXS}
                        ${GFWREG}
                        ${GPHRESXS}
                        ${GFWUTL}
                        ${GTLGEO}
                        ${GTLFLX}
                        ${GPHMEL}
                        ${GTLREW}
                        ${GPHNDCY}
                        ${GFWMSG}
              BASENAME_ONLY
)

simple_plugin(CheckRecoEnergy   "module"
                        nugen_NuReweight
                        nugen_NuReweight_art
                        nusimdata_SimulationBase
                        ${GFWGHEP}
                        ${GFWEG}
                        ${GFWPARDAT}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ${GFWALG}
                        ${GFWPARDAT}
                        ${GPHDISXS}
                        ${GPHDISXS}
                        ${GPHCHMXS}
                        ${GPHCOHEG} ${GPHCOHXS}
                        ${GPHDFRCEG}
                        ${GPHDISEG}
                        ${GPHXSIG}
                        ${GPHDCY}
                        ${GPHQELXS}
                        ${GPHQELXS}
                        ${GFWGHEP}
                        ${GFWEG}
                        ${GPHCMN}
                        ${GFWEG}
                        ${GPHRESXS}
                        ${GPHHADTRANSP}
                        ${GPHHADNZ}
                        ${GFWINT}
                        ${GPHQELXS}
                        ${GPHMNUCEG} ${GPHMNUCXS}
                        ${GPHAMNGEG} ${GPHAMNGXS}
                        ${GPHNUELEG} ${GPHNUELXS}
                        ${GFWNTP}
                        ${GPHNUCLST}
                        ${GFWNUM}
                        ${PHDISXS}
                        ${GFWPARDAT}
                        ${GPHPDF}
                        ${GPHQELEG}
                        ${GPHRESEG} ${GPHRESXS}
                        ${GFWREG}
                        ${GPHRESXS}
                        ${GFWUTL}
                        ${GTLGEO}
                        ${GTLFLX}
                        ${GPHMEL}
                        ${GTLREW}
                        ${GPHNDCY}
                        ${GFWMSG}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
