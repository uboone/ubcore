include_directories($ENV{LARSIM_INC})
link_directories($ENV{LARSIM_LIB})

art_make(
  MODULE_LIBRARIES
    larcoreobj_SummaryData
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    nutools_RandomUtils_NuRandomService_service
    lardataobj_Simulation
    nusimdata_SimulationBase
    larsim_EventGenerator_MARLEY
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_EG}
    ${IFDH}
 )

install_headers()
install_fhicl()
install_source()

