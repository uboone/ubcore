cet_make_library(
  SOURCE
  WireReadoutUBooNEGeom.cxx
  LIBRARIES
  PUBLIC
  larcorealg::Geometry
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(
  UBOpReadoutMap art::service
)

cet_build_plugin(
  UBooNEWireReadout art::service
  LIBRARIES
  PUBLIC
  ubcore::Geometry
  larcore::Geometry_Geometry_service
)

add_subdirectory(gdml)

install_headers()
install_fhicl()
install_source()
