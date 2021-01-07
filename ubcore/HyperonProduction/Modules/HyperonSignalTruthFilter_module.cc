////////////////////////////////////////////////////////////////////////
// Class:       HyperonSignalTruthFilter
// Plugin Type: filter (art v3_01_02)
// File:        HyperonSignalTruthFilter_module.cc
//
// Generated at Mon Nov 23 05:39:28 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"

//local includes

#include "../Alg/FV.h" //fiducial volume checker


namespace hyperon {
	class HyperonSignalTruthFilter;
}


class hyperon::HyperonSignalTruthFilter : public art::EDFilter {
	public:
		explicit HyperonSignalTruthFilter(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		HyperonSignalTruthFilter(HyperonSignalTruthFilter const&) = delete;
		HyperonSignalTruthFilter(HyperonSignalTruthFilter&&) = delete;
		HyperonSignalTruthFilter& operator=(HyperonSignalTruthFilter const&) = delete;
		HyperonSignalTruthFilter& operator=(HyperonSignalTruthFilter&&) = delete;

		// Required functions.
		bool filter(art::Event& e) override;

		// Selected optional functions.
		void beginJob() override;
		void endJob() override;

	private:

		// Declare member data here.

		//create map between particles and their IDs
		std::map<int,art::Ptr<simb::MCParticle>> partByID;
		std::pair<int,art::Ptr<simb::MCParticle>>  part_and_id;


		//used by G4 to track particles
		std::vector<int>daughter_IDs; //ids of semistable hyperon decay products
		//std::vector<int>Sigma0_daughter_IDs; //ids of sigma0 decay products
		std::vector<int>primary_IDs; //ids of particles produced at primary vertex


		//////////////////////////////
		//       Module Labels      //
		//////////////////////////////


		std::string fGenieGenModuleLabel;
		std::string fGeantModuleLabel;

		//////////////////////////////
		//      Misc Parameters     //
		//////////////////////////////

		bool fPrint;



};


hyperon::HyperonSignalTruthFilter::HyperonSignalTruthFilter(fhicl::ParameterSet const& p)
	: EDFilter{p}  // ,
	// More initializers here.
{
	// Call appropriate produces<>() functions here.
	// Call appropriate consumes<>() for any products to be retrieved by this module.

	fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");
	fGeantModuleLabel = p.get<std::string>("GeantLabel");
	fPrint = p.get<bool>("Print",false);
	


}

bool hyperon::HyperonSignalTruthFilter::filter(art::Event& e)
{


	if(fPrint) std::cout << std::endl << "New Event" << std::endl;

	bool IsLambda=false;
	bool IsNuMuBar=false;

	TVector3 PrimaryVertex;
	TVector3 DecayVertex;

	/*
	//Get event generator information

	art::Handle<std::vector<simb::MCTruth>> genMCTruthHandle;
	std::vector<art::Ptr<simb::MCTruth>> genMCTruthVect;
	genMCTruthVect.clear();
	
	std::vector<int> genie_primary_IDs;

	
	if(e.getByLabel(fGenieGenModuleLabel,genMCTruthHandle)) art::fill_ptr_vector(genMCTruthVect,genMCTruthHandle);
	else {	std::cout << "No MC Truth Object in event!" << std::endl; return false; }


	if( genMCTruthVect.size() != 1 ){ std::cout << "N MC Truths != 1" << std::endl; return false; }

	for(const art::Ptr<simb::MCTruth> &mct : genMCTruthVect){

	//get list of mc particles in MC Truth

                for(int k_particle=0;k_particle<mct->NParticles();k_particle++){


                        simb::MCParticle Part = mct->GetParticle(k_particle);

			if(Part.StatusCode() == 1){ genie_primary_IDs.push_back(Part.TrackId()); std::cout << Part.TrackId() << " " << Part.PdgCode() << std::endl; }
	
		}



	}
	
*/

	//Get Geant 4 information

	art::Handle<std::vector<simb::MCParticle>>g4particleHandle;
	std::vector< art::Ptr<simb::MCParticle>>g4particleVect;
	g4particleVect.clear();



	if(e.getByLabel(fGeantModuleLabel,g4particleHandle)) art::fill_ptr_vector(g4particleVect,g4particleHandle);
	else
		std::cout << "Geant Particles Missing!" << std::endl;

	//id numbers of hyperon daughters
	daughter_IDs.clear();

	//id numbers of particles produced at primary vertex
	primary_IDs.clear();

	//id nubmers of Sigma Zero decay products
	//Sigma0_daughter_IDs.clear();

	//map between particle ID numbers (g4p->TrackId()) and pointers to simb::MCParticle 
	partByID.clear();


	for(const art::Ptr<simb::MCParticle> &g4p : g4particleVect){

		if(g4p->Mother() == 0){

			primary_IDs.push_back(g4p->TrackId());

		

			if(g4p->PdgCode() == 3122 && g4p->EndProcess() == "Decay"){
				IsLambda = true;
				//get the hyperon's daughters
				for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

					daughter_IDs.push_back( g4p->Daughter(i_d) );

				}
				// Get decay vertex
				DecayVertex = { g4p->EndPosition().X() , g4p->EndPosition().Y() , g4p->EndPosition().Z() };
			

			}


		}

		part_and_id = std::make_pair(g4p->TrackId() , g4p);

		partByID.insert( part_and_id );

	}



	//if there is no lambda in final state, reject event
	if(!IsLambda){ if(fPrint) std::cout << "Event does not contain Lambda in final state" << std::endl; return false; }

	//search list of primary IDs for antimuon

	for(size_t i_p=0;i_p<primary_IDs.size();i_p++){

		if(partByID.find(primary_IDs[i_p]) == partByID.end()) continue;

		art::Ptr<simb::MCParticle> part = partByID[primary_IDs.at(i_p)];

		if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these


		if(part->PdgCode() == -13){
			IsNuMuBar = true;
			PrimaryVertex = { part->Position().X() , part->Position().Y() , part->Position().Z() };

		}


	}



	//check fiducial volume
	if(!inActiveTPC(PrimaryVertex)){ if(fPrint) std::cout << "Outside fiducial volume" << std::endl; return false; }


	//check event is muon antineutrino
	if(!IsNuMuBar){ std::cout << "Not a NuMuBar event" << std::endl; return false; }

	//check decay products
	bool GotProton=false;
	bool GotPion=false;
	
	int nProducts=0;
	
	for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){

         if(partByID.find(daughter_IDs[i_d]) == partByID.end()) continue;

                art::Ptr<simb::MCParticle> part = partByID[daughter_IDs[i_d]];

		if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

		TVector3 StartPosition( part->Position().X() , part->Position().Y() , part->Position().Z() ); 
		
		//ignore any particles not created at the decay vertex
		if( StartPosition != DecayVertex ) continue;

		nProducts++;		
				
		if( part->PdgCode() == 2212 ) GotProton = true;
		if( part->PdgCode() == -211 ) GotPion = true;


}

	if(nProducts != 2 && fPrint) std::cout << "Num decay products != 2" << std::endl;
	if(!GotProton  && fPrint) std::cout << "No true decay proton found" << std::endl;
	if(!GotPion  && fPrint) std::cout << "No true decay pi minus found" << std::endl;

	if(nProducts != 2 || !GotProton || !GotPion) return false;


	return true;
}

void hyperon::HyperonSignalTruthFilter::beginJob()
{
	// Implementation of optional member function here.
}

void hyperon::HyperonSignalTruthFilter::endJob()
{
	// Implementation of optional member function here.
}

DEFINE_ART_MODULE(hyperon::HyperonSignalTruthFilter)
