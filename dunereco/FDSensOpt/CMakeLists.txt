add_subdirectory(IniSegAlg)
add_subdirectory(MVAAlg)
add_subdirectory(FDSensOptData)
add_subdirectory(NeutrinoEnergyRecoAlg)

include_directories ( $ENV{GENIE_INC}/GENIE )

# genie
  cet_find_library( GFWMSG       NAMES GFwMsg       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWREG       NAMES GFwReg       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWALG       NAMES GFwAlg       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWINT       NAMES GFwInt       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWGHEP      NAMES GFwGHEP      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWNUM       NAMES GFwNum       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWUTL       NAMES GFwUtl       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWPARDAT    NAMES GFwParDat    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWEG        NAMES GFwEG        PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWNTP       NAMES GFwNtp       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHXSIG      NAMES GPhXSIg      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHPDF       NAMES GPhPDF       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNUCLST    NAMES GPhNuclSt    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCMN       NAMES GPhCmn       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDCY       NAMES GPhDcy       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHHADTRANSP NAMES GPhHadTransp PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHHADNZ     NAMES GPhHadnz     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDEEX      NAMES GPhDeEx      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHAMNGXS    NAMES GPhAMNGXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHAMNGEG    NAMES GPhAMNGEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCHMXS     NAMES GPhChmXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCOHXS     NAMES GPhCohXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCOHEG     NAMES GPhCohEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDISXS     NAMES GPhDISXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDISEG     NAMES GPhDISEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDFRCXS    NAMES GPhDfrcXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDFRCEG    NAMES GPhDfrcEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHGLWRESXS  NAMES GPhGlwResXS  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHGLWRESEG  NAMES GPhGlwResEG  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHIBDXS     NAMES GPhIBDXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHIBDEG     NAMES GPhIBDEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHMNUCXS    NAMES GPhMNucXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHMNUCEG    NAMES GPhMNucEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHMEL       NAMES GPhMEL       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNUELXS    NAMES GPhNuElXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNUELEG    NAMES GPhNuElEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHQELXS     NAMES GPhQELXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHQELEG     NAMES GPhQELEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHRESXS     NAMES GPhResXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHRESEG     NAMES GPhResEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHSTRXS     NAMES GPhStrXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHSTREG     NAMES GPhStrEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNDCY      NAMES GPhNDcy      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GTLGEO       NAMES GTlGeo       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GTLFLX       NAMES GTlFlx       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GRWFWK       NAMES GRwFwk       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GRWIO        NAMES GRwIO        PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GRWCLC       NAMES GRwClc       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )

  # GENIE pre-R-3 used GSL via ROOT's Math libs, now makes direct call not available via that route
  cet_find_library( GSL NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )

  set (GENIE_LIB_LIST GFWMSG
                      GFWREG
                      GFWALG
                      GFWINT
                      GFWGHEP
                      GFWNUM
                      GSL  # FWNUM relies on GSL
                      GFWUTL
                      GFWPARDAT
                      GFWEG
                      GFWNTP
                      GPHXSIG
                      GPHPDF
                      GPHNUCLST
                      GPHCMN
                      GPHDCY
                      GPHHADTRANSP
                      GPHHADNZ
                      GPHDEEX
                      GPHAMNGXS
                      GPHAMNGEG
                      GPHCHMXS
                      GPHCOHXS
                      GPHCOHEG
                      GPHDISXS
                      GPHDISEG
                      GPHDFRCXS
                      GPHDFRCEG
                      GPHGLWRESXS
                      GPHGLWRESEG
                      GPHIBDXS
                      GPHIBDEG
                      GPHMNUCXS
                      GPHMNUCEG
                      GPHMEL
                      GPHNUELXS
                      GPHNUELEG
                      GPHQELXS
                      GPHQELEG
                      GPHRESXS
                      GPHRESEG
                      GPHSTRXS
                      GPHSTREG
                      GPHNDCY
                      GTLGEO
                      GTLFLX
                      GRWFWK
                      GRWIO
                      GRWCLC
                      )


cet_build_plugin(NueAna   art::module
                        lardataalg_DetectorInfo
		        lardataobj_RawData
		        lardataobj_RecoBase
		        lardata_ArtDataHelper
			larreco_Calorimetry
		        larreco_RecoAlg_PMAlg
		        IniSegAlg
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
 			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			
			cetlib::cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

cet_build_plugin(MVASelect   art::module
                        MVAAlg
                        nugen::NuReweight
                        nugen::NuReweight_art
			nusimdata::SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
			${GENIE_LIB_LIST}
              BASENAME_ONLY
)

cet_build_plugin(IniSegReco art::module
                      lardataalg_DetectorInfo
		      lardataobj_RawData
		      lardataobj_RecoBase
		      larreco_RecoAlg
		      larreco_RecoAlg_PMAlg
		      IniSegAlg
		      larcorealg_Geometry
		      larcore_Geometry_Geometry_service
		      lardata_Utilities
		      larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
		      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
		      nusimdata::SimulationBase
      		      ART_FRAMEWORK_CORE
		      ART_FRAMEWORK_PRINCIPAL
		      ART_FRAMEWORK_SERVICES_REGISTRY
      		      ART_ROOT_IO_TFILE_SUPPORT
      		      ROOT::Core
      		      ART_ROOT_IO_TFILESERVICE_SERVICE
 		      art::Persistency_Common canvas
		      art::Persistency_Provenance canvas
		      art::Utilities canvas
		      messagefacility::MF_MessageLogger
		      
		      cetlib::cetlib cetlib_except
		      ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

cet_build_plugin(ShSeg   art::module
                        lardataalg_DetectorInfo
			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			IniSegAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			
			
			larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
			ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
 			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			
			cetlib::cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(EnergyReco   art::module
                        MVAAlg
                        NeutrinoEnergyRecoAlg
                        nugen::NuReweight
                        nugen::NuReweight_art
                        nusimdata::SimulationBase
                        larcore_Geometry_Geometry_service
                        larreco_Calorimetry
                        larreco_RecoAlg
                        lardataalg_DetectorInfo
                        lardataobj_RecoBase
                        larsim_MCCheater_BackTrackerService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                         ${GENIE_LIB_LIST}
              BASENAME_ONLY
)

cet_build_plugin(CheckRecoEnergy   art::module
                        nugen::NuReweight
                        nugen::NuReweight_art
                        nusimdata::SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
                        ${GENIE_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
