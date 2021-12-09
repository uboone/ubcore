
art_make( BASENAME_ONLY
          LIBRARY_NAME uboonecode_DQMTools
          LIB_LIBRARIES
			lardataobj_RawData
			lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			art::Persistency_Common canvas::canvas
			art::Persistency_Provenance canvas::canvas
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except

          MODULE_LIBRARIES
			uboonecode_DQMTools
                        lardataobj_RecoBase
			lardataobj_MCBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        nug4::ParticleNavigation
                        lardataobj_Simulation
                        nusimdata::SimulationBase
                        art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
			art::Persistency_Common canvas::canvas
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        ROOT::Tree
         )

install_headers()
install_fhicl()
install_source()
