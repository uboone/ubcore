

art_make(

		MODULE_LIBRARIES
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${MF_MESSAGELOGGER}
		
		nusimdata_SimulationBase

                art_Persistency_Common canvas
                art_Persistency_Provenance canvas
                art_Utilities canvas

		${ROOT_BASIC_LIB_LIST}


	)




install_headers()
install_fhicl()
install_source()

