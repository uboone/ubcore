link_directories( $ENV{DK2NUDATA_LIB})
link_directories( $ENV{GENIE_LIB})
link_directories( $ENV{LOG4CPP_LIB})

art_make(BASENAME_ONLY
          LIBRARY_NAME uboonecode_HiggsPortalScalarGen
          LIB_LIBRARIES
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom::RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            
            art::Framework_Services_Optional_RandomNumberGenerator_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            cetlib::cetlib
            cetlib_except::cetlib_except
            messagefacility::MF_MessageLogger
            
            fhiclcpp::fhiclcpp
            CLHEP::CLHEP
            SQLite::SQLite3
            ROOT::EG
            ROOT::Tree
            ROOT::Geom
	    ifdh::ifdh
            
            dk2nuTree
            GFwMsg
	    GFwParDat
	    log4cpp
      MODULE_LIBRARIES
            uboonecode_HiggsPortalScalarGen
            larcoreobj_SummaryData
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            nurandom::RandomUtils_NuRandomService_service
            lardataobj_Simulation
            nusimdata::SimulationBase
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
	    art_root_io::tfile_support ROOT::Core
            art_root_io::TFileService_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            cetlib::cetlib
            cetlib_except::cetlib_except
            messagefacility::MF_MessageLogger
            
            fhiclcpp::fhiclcpp
            CLHEP::CLHEP
            SQLite::SQLite3
            ROOT::EG
            ROOT::Geom
	    ifdh::ifdh
	    dk2nuTree
            GFwMsg
	    GFwParDat
	    log4cpp
            )

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

