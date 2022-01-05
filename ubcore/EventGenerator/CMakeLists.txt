install_fhicl()
install_source()
add_subdirectory(GENIE)
add_subdirectory(HeavySterileNeutrino)
add_subdirectory(HiggsPortalScalarGen)

art_make( BASENAME_ONLY
          MODULE_LIBRARIES
		   larcorealg_Geometry
		   ubcore_Geometry
		   larcore_Geometry_Geometry_service
		   nug4::ParticleNavigation lardataobj_Simulation
		   nusimdata::SimulationBase
		   lardata_Utilities
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   nurandom::RandomUtils_NuRandomService_service
		   art::Framework_Core
		   art::Framework_Principal
		   art::Framework_Services_Registry
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
                   art::Framework_Services_Optional_RandomNumberGenerator_service
		   art::Persistency_Common
		   art::Persistency_Provenance
		   art::Utilities
		   canvas::canvas
		   messagefacility::MF_MessageLogger
		   fhiclcpp::fhiclcpp
		   CLHEP::CLHEP
		   cetlib::cetlib
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml
                   ROOT::EG
		)
