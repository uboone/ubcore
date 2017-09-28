////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackProducer
// Module Type: producer
// File:        CRTTrackProducer_module.cc
// Description: Module for constructiong CRT tracks.
// Generated at Tue May 16 08:31:48 2017 by David Lorca Galindo using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
//#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "artdaq-core/Data/Fragments.hh"

#include "art/Framework/Services/Optional/TFileService.h"

#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "uboone/CRT/CRTProducts/CRTTrack.hh"
#include "uboone/CRT/CRTAuxFunctions.hh"

#include "TTree.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3S.h"
#include "TProfile.h"
#include "TF1.h"
#include "TDatime.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

namespace bernfebdaq {
  class CRTTrackProducer;
}

class bernfebdaq::CRTTrackProducer : public art::EDProducer {
public:
  explicit CRTTrackProducer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTrackProducer(CRTTrackProducer const &) = delete;
  CRTTrackProducer(CRTTrackProducer &&) = delete;
  CRTTrackProducer & operator = (CRTTrackProducer const &) = delete;
  CRTTrackProducer & operator = (CRTTrackProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_label_;
  double max_time_difference_ ;//max time for coincidence 
  int store_track_;



  //quallity plots
  TTree* my_tree_;     
  double track_time_ns = -1e18;
  double track_time_s = -1e18;
  double time_diff = 1e24;
  double length = -1e18;
  double theta = -1e18;
  double phi = -1e18;

  TH2F* hplavspla;
  TH1F* hTlength;
  TH2F* hTlengthvsTime;
  TH1F* htheta;
  TH1F* hphi;
  TH1F* hts0_ns;
  TH2F* hTvsH;

 //quallity plots                                                                                                                                                      

  //test                                                                                                                                         
  int verbose_ = 0;
  //test

};


bernfebdaq::CRTTrackProducer::CRTTrackProducer(fhicl::ParameterSet const & p)
  :
  // Initialize member data here.
  data_label_(p.get<std::string>("data_label")),
  max_time_difference_(p.get<double>("max_time_difference")),
  store_track_(p.get<int>("store_track")),
  verbose_(p.get<int>("verbose"))
  
{
  // Call appropriate produces<>() functions here.
  if(store_track_ == 1) 
    produces< std::vector<crt::CRTTrack>   >();
}

void bernfebdaq::CRTTrackProducer::produce(art::Event & evt)
{
  // Implementation of required member function here.
  
  art::Handle< std::vector<crt::CRTHit> > rawHandle;
  evt.getByLabel(data_label_, rawHandle); //what is the product instance name? no BernZMQ
  
  //check to make sure the data we asked for is valid                                                                                                      
  if(!rawHandle.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_label_ << std::endl;
    std::cout << std::endl;
    return;
  }
  
  //get better access to the data               
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle);
  
  //CRTTrack collection on this event                                                                                   
  std::unique_ptr<std::vector<crt::CRTTrack> > CRTTrackCol(new std::vector<crt::CRTTrack>);
  
  // std::cout<<"  CRTHitCollection.size()  "<<CRTHitCollection.size()<<std::endl;
  //getchar();
  
  int N_CRTHits = CRTHitCollection.size();
  int N_CRTTrack = 0;
  
  for(std::vector<int>::size_type i = 0; i != CRTHitCollection.size(); i++) {//A //comparar solo con el resto del vector
    
    crt::CRTHit CRTHiteventA = CRTHitCollection[i];	
    
    for(std::vector<int>::size_type j = 0; j != CRTHitCollection.size(); j++) {//B
      if(j>i){//C
	
	crt::CRTHit CRTHiteventB = CRTHitCollection[j];
	
	//look for coincidences
	double time_s_A = CRTHiteventA.ts0_s;
	double time_ns_A = CRTHiteventA.ts0_ns;
	double time_s_B = CRTHiteventB.ts0_s;
	double time_ns_B = CRTHiteventB.ts0_ns;
	int planeA = CRTHiteventA.plane; 
	int planeB = CRTHiteventB.plane; 
	double time_diffABS = abs(time_ns_A - time_ns_B);
	time_diff = time_ns_A - time_ns_B;
	
	
	if( (planeA != planeB) && (time_s_A == time_s_B) && (time_diffABS<max_time_difference_)  ){//D
	  
	  N_CRTTrack++;
	  
	  double xA = CRTHiteventA.x_pos;
	  double yA = CRTHiteventA.y_pos;
	  double zA = CRTHiteventA.z_pos;
	  double xB = CRTHiteventB.x_pos;
	  double yB = CRTHiteventB.y_pos;
	  double zB = CRTHiteventB.z_pos;
	  
	  length = sqrt(pow(xA-xB,2)+pow(yA-yB,2)+pow(zA-zB,2));
	  theta = 1/( std::abs(xA-xB) / std::abs(yA-yB));
	  phi = (1/( std::abs(zA-zB) / std::abs(yA-yB) ) );
	  track_time_ns = (CRTHiteventA.ts0_ns + CRTHiteventB.ts0_ns)/2;
	  track_time_s = CRTHiteventA.ts0_s;
	  
	  if(verbose_ == 1){
	    std::cout.precision(19);
	    std::cout<<"tAs: "<<time_s_A<< " s" <<std::endl;
	    std::cout<<"tAns: "<<time_ns_A<< " ns" <<std::endl;
	    std::cout<<"tBs: "<<time_s_B<< " s" <<std::endl;
	    std::cout<<"tBns: "<<time_ns_B<< " ns" <<std::endl;
	    std::cout<<"planeA: "<<planeA<<std::endl;
	    std::cout<<"planeB: "<<planeB<<std::endl;
	    std::cout<<"time_diff: "<<time_diff<< " ns" <<std::endl;
	    std::cout<<"length: "<<length<< " cm" <<std::endl;
	    std::cout<<"theta: "<<theta<< " degrees in xz" <<std::endl;
	    std::cout<<"phi: "<<phi<< " degrees in yz" <<std::endl;
	    
	    getchar();
	  }
	  
	  //do something
	  
	  crt::CRTTrack CRTcanTrack;
	  CRTcanTrack.x1_pos = xA;
	  CRTcanTrack.y1_pos = yA;
	  CRTcanTrack.z1_pos = zA;
	  CRTcanTrack.x2_pos = xB;
	  CRTcanTrack.y2_pos = yB;
	  CRTcanTrack.z2_pos = zB;
	  CRTcanTrack.ts0_s = time_s_A;
	  CRTcanTrack.ts0_ns = track_time_ns;
	  CRTcanTrack.ts0_ns_1 = time_ns_A; 
	  CRTcanTrack.ts0_ns_2 = time_ns_B; 

	  CRTcanTrack.length = length;
	  CRTcanTrack.thetaxy = theta;
	  CRTcanTrack.phizy = phi;
	
	  CRTTrackCol->emplace_back(CRTcanTrack);
	  
	  ////quality plot
	  my_tree_->Fill();
	  hplavspla->Fill(planeA,planeB);	
	  hTlength->Fill(length);
	  hTlengthvsTime->Fill(length,time_diff);
	  htheta->Fill(theta);
	  hphi->Fill(phi);
	  hts0_ns->Fill(track_time_ns);
	  //quallity plots                                                                                                                                                 		
	}//D
      }//C
    }//B
  }//A
  
  hTvsH->Fill(N_CRTHits,N_CRTTrack);
  
  //store track collection into event
  if(store_track_ == 1)
    evt.put(std::move(CRTTrackCol));
  
}

void bernfebdaq::CRTTrackProducer::beginJob()
{

  my_tree_ = tfs->make<TTree>("my_tree","CRT Tree");
  my_tree_->Branch("track_time_s", &track_time_s, "time (s)");
  my_tree_->Branch("track_time_ns", &track_time_ns, "time (ns)");
  my_tree_->Branch("track_length", &length, "Track lenght (cm)");
  my_tree_->Branch("track_timdiff", &time_diff, "Time_diff (ns)");
  my_tree_->Branch("track_theta", &theta, "Theta_xy (º)");
  my_tree_->Branch("track_phi", &phi, "Phi_xy (º)");

  hplavspla = tfs->make<TH2F>("hplavspla","PlanevsPlane",4,0,3,4,0,3);
  hplavspla->GetXaxis()->SetTitle("Plane (0=Bottom, 1=FT, 2=Pipe, 3=Top)");
  hplavspla->GetYaxis()->SetTitle("Plane (0=Bottom, 1=FT, 2=Pipe, 3=Top)");
  hplavspla->GetZaxis()->SetTitle("Entries/bin");
  hplavspla->SetOption("COLZ");

  hTvsH = tfs->make<TH2F>("hTvsH","Track_vs_Hits",500,0,500,500,0,500);
  hTvsH->GetXaxis()->SetTitle("Number of CRTHits per event");
  hTvsH->GetYaxis()->SetTitle("Number of CRTTracks per event");
  hTvsH->GetZaxis()->SetTitle("Entries/bin");
  hTvsH->SetOption("COLZ");

  hTlength = tfs->make<TH1F>("hTlength","Track_Length",1500,0,1500);
  hTlength->GetXaxis()->SetTitle("Thack_Length (cm)");
  hTlength->GetYaxis()->SetTitle("Entries/bin");

  hTlengthvsTime = tfs->make<TH2F>("hTlengthvsTime","Track_LengthvsTime",1500,0,1500,200,-100,100);
  hTlengthvsTime->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTime->GetYaxis()->SetTitle("Track_time (ns)");
  hTlengthvsTime->GetZaxis()->SetTitle("Entries/bin");
  hTlengthvsTime->SetOption("COLZ");

  htheta = tfs->make<TH1F>("htheta","Track_theta",900,0,180);
  htheta->GetXaxis()->SetTitle("Theta_xy (º)");
  htheta->GetYaxis()->SetTitle("Entries/bin");
 
  hphi = tfs->make<TH1F>("hphi","Track_phi",900,0,180);
  hphi->GetXaxis()->SetTitle("Phi_zy (º)");
  hphi->GetYaxis()->SetTitle("Entries/bin");

  hts0_ns = tfs->make<TH1F>("hts0_ns","Track_time_ns",100000,0,1e9);
  hts0_ns->GetXaxis()->SetTitle("Track time (ns)");
  hts0_ns->GetYaxis()->SetTitle("Entries/bin");


}

void bernfebdaq::CRTTrackProducer::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(bernfebdaq::CRTTrackProducer)
