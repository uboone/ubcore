art_make(
          LIBRARY_NAME ubcore_SampledInput
          # Needed for standalone c++ helper files
          LIB_LIBRARIES

          # Needed for the main module
          MODULE_LIBRARIES
                    larcoreobj_SummaryData
                    nusimdata::SimulationBase
                    ${ART_FRAMEWORK_CORE}
                    ${ART_FRAMEWORK_PRINCIPAL}
                    ${ART_FRAMEWORK_SERVICES_REGISTRY}
                    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                    #${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		    ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                    art_Framework_Services_Registry
                    art_Persistency_Common
                    art_Persistency_Provenance
                    art_Utilities
                    canvas
                    ${MF_MESSAGELOGGER}
                    ${FHICLCPP}
                    cetlib cetlib_except
                    ROOT::Tree
          )

install_headers()
install_source()
install_fhicl()

