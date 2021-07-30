
art_make( BASENAME_ONLY
          LIBRARY_NAME uboonecode_DQMTools
          LIB_LIBRARIES
			lardataobj_RawData
			lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except

          MODULE_LIBRARIES
			uboonecode_DQMTools
                        lardataobj_RecoBase
			lardataobj_MCBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        nug4_ParticleNavigation lardataobj_Simulation nug4_ParticleNavigation lardataobj_Simulation
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}

	                ${FHICLCPP}
			cetlib cetlib_except

         )

install_headers()
install_fhicl()
install_source()
