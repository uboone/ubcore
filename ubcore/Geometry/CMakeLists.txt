
art_make( BASENAME_ONLY
	  LIBRARY_NAME ubcore_Geometry
	  LIB_LIBRARIES larcorealg_Geometry
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib cetlib_except
                        fhiclcpp::fhiclcpp

	  SERVICE_LIBRARIES larcorealg_Geometry
	  		ubcore_Geometry

	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml

        )

add_subdirectory(gdml)

install_headers()
install_fhicl()
install_source()


