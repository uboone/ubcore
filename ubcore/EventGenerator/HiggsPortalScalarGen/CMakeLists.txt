link_directories( $ENV{DK2NUDATA_LIB})

art_make(BASENAME_ONLY
          LIBRARY_NAME uboonecode_HiggsPortalScalarGen
          LIB_LIBRARIES
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom_RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            art_Persistency_Common canvas
            art_Persistency_Provenance canvas
            art_Utilities canvas
            cetlib cetlib_except
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            CLHEP::CLHEP
            ${ROOT_BASIC_LIB_LIST}
            ${SQLITE3}
            ${ROOT_EG}
            ${ROOT_GEOM}
	    ifdh
            ${Boost_SYSTEM_LIBRARY}
            dk2nuTree
            GFwMsg
	    GFwParDat
	    log4cpp
      MODULE_LIBRARIES
            uboonecode_HiggsPortalScalarGen
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom_RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
	    ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
            art_Persistency_Common canvas
            art_Persistency_Provenance canvas
            art_Utilities canvas
            cetlib cetlib_except
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            CLHEP::CLHEP
            ${ROOT_BASIC_LIB_LIST}
            ${SQLITE3}
            ${ROOT_EG}
            ${ROOT_GEOM}
	    ifdh
	    dk2nuTree
            GFwMsg
	    GFwParDat
	    log4cpp
            )

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

