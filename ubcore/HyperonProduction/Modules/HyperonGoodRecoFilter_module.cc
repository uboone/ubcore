////////////////////////////////////////////////////////////////////////
// Class:       HyperonGoodRecoFilter
// Plugin Type: filter (art v3_01_02)
// File:        HyperonGoodRecoFilter_module.cc
//
// Generated at Mon Nov 23 12:12:59 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_05_01.
//
// Filter to find Lambda analysis signal events and check if hyperon decay was
// reconstructed
//
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
#include <vector>


#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"


#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"


#include "lardataobj/RecoBase/PFParticle.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"



//local includes
#include "../Alg/FV.h" //fiducial volume checker


namespace hyperon {
	class HyperonGoodRecoFilter;
}


class hyperon::HyperonGoodRecoFilter : public art::EDFilter {
	public:
		explicit HyperonGoodRecoFilter(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		HyperonGoodRecoFilter(HyperonGoodRecoFilter const&) = delete;
		HyperonGoodRecoFilter(HyperonGoodRecoFilter&&) = delete;
		HyperonGoodRecoFilter& operator=(HyperonGoodRecoFilter const&) = delete;
		HyperonGoodRecoFilter& operator=(HyperonGoodRecoFilter&&) = delete;

		// Required functions.
		bool filter(art::Event& e) override;

		// Selected optional functions.
		void beginJob() override;
		void endJob() override;

	private:

		// Declare member data here.


		/////////////////////////////
		//         G4 Info         //
		/////////////////////////////

		bool fIsSignal=false;
		bool fHyperonGoodReco=false;

		bool IsLambda=false;
		bool IsNuMuBar=false;

		TVector3 PrimaryVertex;
		TVector3 DecayVertex;

		//track IDs of decay products, used to find tracks truth matching to them later
		int Proton_TrackID;				
		int Pion_TrackID;				


		//create map between particles and their IDs
		std::map<int,art::Ptr<simb::MCParticle>> partByID;
		std::pair<int,art::Ptr<simb::MCParticle>>  part_and_id;


		//used by G4 to track particles
		std::vector<int>daughter_IDs; //ids of semistable hyperon decay products
		//std::vector<int>Sigma0_daughter_IDs; //ids of sigma0 decay products
		std::vector<int>primary_IDs; //ids of particles produced at primary vertex



		////////////////////////////
		//        Metadata        //
		////////////////////////////

		int nPassedTruth;
		int nPassedReco;	

		/////////////////////////////
		//      Module Labels      //
		/////////////////////////////

		std::string fGeantModuleLabel;
		std::string fTrackLabel;
		std::string fPFParticleLabel;
		std::string fCaloLabel;
		std::string fHitLabel;
		std::string fTrackHitAssnLabel;
		std::string fHitTruthAssnLabel;

		//////////////////////////////
		//      Misc Parameters     //
		//////////////////////////////

		bool fPrint;

};


hyperon::HyperonGoodRecoFilter::HyperonGoodRecoFilter(fhicl::ParameterSet const& p)
	: EDFilter{p}  // ,
	// More initializers here.
{

	fPFParticleLabel = p.get<std::string>("PFParticleLabel");
	fGeantModuleLabel = p.get<std::string>("GeantLabel");
	fTrackLabel = p.get<std::string>("TrackLabel");
	fCaloLabel = p.get<std::string>("CaloLabel");
	fHitLabel  = p.get<std::string>("HitLabel");
	fTrackHitAssnLabel = p.get<std::string>("TrackHitAssnLabel");
	fHitTruthAssnLabel = p.get<std::string>("HitTruthAssnLabel");

	fPrint = p.get<bool>("Print",false);

}

bool hyperon::HyperonGoodRecoFilter::filter(art::Event& e)
{


	if(fPrint) std::cout << "New Event" << std::endl;

	//reset everything 
	fIsSignal = false;
	fHyperonGoodReco = false;

	Proton_TrackID = -1000;	
	Pion_TrackID = -1000;	



	//first check if event is signal - discard it if it isn't


	bool IsLambda=false;
	bool IsNuMuBar=false;

	TVector3 PrimaryVertex;
	TVector3 DecayVertex;


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
	bool GotProtonTruth=false;
	bool GotPionTruth=false;

	int nProducts=0;

	for(size_t i_d=0;i_d<daughter_IDs.size();i_d++){

		if(partByID.find(daughter_IDs[i_d]) == partByID.end()) continue;

		art::Ptr<simb::MCParticle> part = partByID[daughter_IDs[i_d]];

		if(part->PdgCode() > 10000) continue; //anything with very large pdg code is a nucleus, skip these

		TVector3 StartPosition( part->Position().X() , part->Position().Y() , part->Position().Z() ); 

		//ignore any particles not created at the decay vertex
		if( StartPosition != DecayVertex ) continue;

		nProducts++;		

		if( part->PdgCode() == 2212 ){ GotProtonTruth = true; Proton_TrackID = part->TrackId(); }
		if( part->PdgCode() == -211 ){ GotPionTruth = true; Pion_TrackID = part->TrackId();  }

	}


	if(nProducts != 2 && fPrint) std::cout << "Num decay products != 2" << std::endl;
	if(!GotProtonTruth  && fPrint) std::cout << "No true decay proton found" << std::endl;
	if(!GotPionTruth  && fPrint) std::cout << "No true decay pi minus found" << std::endl;

	if(nProducts != 2 || !GotProtonTruth || !GotPionTruth) return false;
	
	nPassedTruth++;


	/////////////////////////////////////////////
	//       Analyse Reconstructed Info        //
	/////////////////////////////////////////////


	bool GotProtonReco=false;
	bool GotPionReco=false;

	//setup handles
	art::Handle<std::vector<recob::PFParticle>>pfparticleHandle; //PFParticles
	art::Handle<std::vector<recob::Track>  > trackHandle; //reconstructed tracks
	art::Handle<std::vector<anab::Calorimetry>>caloHandle;
	art::Handle<std::vector<recob::Hit> > hitHandle;


	std::vector< art::Ptr<recob::PFParticle>>pfparticleVect; //vector of PFparticles
	std::vector <art::Ptr <recob::Track>  > trackVect; //vector of tracks
	std::vector<art::Ptr<anab::Calorimetry>>caloVect; //calorimetry
	std::vector <art::Ptr <recob::Hit>  > hitVect; //hits


	//setup handles
	if(e.getByLabel(fPFParticleLabel,pfparticleHandle)) art::fill_ptr_vector(pfparticleVect,pfparticleHandle);
	else return false;

	if(!pfparticleVect.size()) return false;


	//fill track vector
	if(e.getByLabel(fTrackLabel,trackHandle)) art::fill_ptr_vector(trackVect,trackHandle);
	else std::cout << "Track handle not setup" << std::endl;




	e.getByLabel(fHitLabel,hitHandle);	

	//setup assns
	
	art::FindManyP<recob::Track> trackAssoc(pfparticleVect,e,fTrackLabel);
	art::FindManyP<anab::Calorimetry> caloTrackAssoc(trackVect,e,fCaloLabel);

	art::FindManyP<recob::Hit> trackHitAssoc(trackHandle,e,fTrackHitAssnLabel);

	art::FindMany< simb::MCParticle , anab::BackTrackerHitMatchingData> particlesPerHit(hitHandle,e,fHitTruthAssnLabel);

	//check if there are tracks truth matching to decay products
	for(const art::Ptr<recob::PFParticle> &pfp : pfparticleVect){


		//skip anything that isn't a track
		if(pfp->PdgCode() != 13) continue;

		std::vector< art::Ptr<recob::Track> > pfpTracks = trackAssoc.at(pfp.key());

		if(pfpTracks.size() == 1){

			for(const art::Ptr<recob::Track> &trk : pfpTracks){			


				//Truth match 
				std::vector< art::Ptr< recob::Hit> > hits = trackHitAssoc.at(trk.key());

				std::unordered_map<int , double>  trkide;
				//			double maxe = -1;
				double tote = 0;

				int maxhits=-1;

				simb::MCParticle const* matchedParticle = NULL;

				std::vector< simb::MCParticle const*> particleVec;
				std::vector< anab::BackTrackerHitMatchingData const*> matchVec;


				for (size_t i_hit = 0; i_hit < hits.size(); ++i_hit){

					//clear vectors
					particleVec.clear();
					matchVec.clear();
					particlesPerHit.get(hits[i_hit].key(), particleVec, matchVec);



					//loop over particles that deposit energy in this hit
					for (size_t i_particle = 0; i_particle < particleVec.size(); ++i_particle){

						//      trkide[ particleVec[i_particle]->TrackId() ] += matchVec[i_particle]->energy;
						trkide[ particleVec[i_particle]->TrackId() ] ++; //just increment the number of hits

						tote += matchVec[i_particle]->energy;

						/*
						   if ( trkide[ particleVec[i_particle]->TrackId() ] > maxe ){
						   maxe = trkide[ particleVec[i_particle]->TrackId() ];
						   matchedParticle = particleVec[i_particle];

						   }
						   */



						if ( trkide[ particleVec[i_particle]->TrackId() ] > maxhits ){
							maxhits = trkide[ particleVec[i_particle]->TrackId() ];
							matchedParticle = particleVec[i_particle];

						}



					}


				}


				if( matchedParticle == NULL ) continue;

				if( matchedParticle->TrackId() == Proton_TrackID ) GotProtonReco = true; 
				if( matchedParticle->TrackId() == Pion_TrackID ) GotPionReco = true; 

			}

		}		





	}



if(GotProtonReco && GotPionReco){ nPassedReco++; return true; }
	else return false;


}

void hyperon::HyperonGoodRecoFilter::beginJob()
{
	// Implementation of optional member function here.

	nPassedTruth=0;
	nPassedReco=0;


}


void hyperon::HyperonGoodRecoFilter::endJob()
{
	// Implementation of optional member function here.
	
	std::cout << "Events passing truth filter: " << nPassedTruth << std::endl;
	std::cout << "Events passing reco filter: " << nPassedReco << std::endl;
}

DEFINE_ART_MODULE(hyperon::HyperonGoodRecoFilter)
