simple_plugin(LongestRecoVertexTrackSelector "tool"
dunereco_AnaUtils
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
		       	   ${ART_FRAMEWORK_SERVICES_REGISTRY}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas
             )

simple_plugin(HighestPandrizzleScoreRecoVertexShowerSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
		       	   ${ART_FRAMEWORK_SERVICES_REGISTRY}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandrizzleAlg    
                #dune_FDSelections_pandizzle
             )


simple_plugin(HighestPandizzleScoreRecoVertexTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
		       	   ${ART_FRAMEWORK_SERVICES_REGISTRY}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas
                   PandizzleAlg
             )

simple_plugin(HighestEnergyRecoVertexShowerSelector "tool"
    	  	        lardataobj_RecoBase
                larcorealg_Geometry
                larcore_Geometry_Geometry_service
                larreco_RecoAlg
    	  	        larpandora_LArPandoraInterface
dunereco_AnaUtils
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
		       	${ART_FRAMEWORK_SERVICES_REGISTRY}
	  	       	${ART_FRAMEWORK_SERVICES_OPTIONAL}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas
             )
