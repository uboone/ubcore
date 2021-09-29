
art_make(
          LIB_LIBRARIES
                        ${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except

          MODULE_LIBRARIES
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}

	                ${FHICLCPP}
			cetlib cetlib_except

         )

install_headers()
install_fhicl()
install_source()
