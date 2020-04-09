////////////////////////////////////////////////////////////////////////
// Class:       SampledInputInfo
// Plugin Type: producer (art v3_05_00)
// File:        SampledInputInfo_module.cc
//
// Generated at Wed Apr  8 13:35:20 2020 by Wesley Ketchum using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/Sampled.h"
#include "canvas/Persistency/Provenance/SampledInfo.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"

#include <memory>


//include the truth objects
#include "larcoreobj/SummaryData/POTSummary.h"

namespace util {
  class SampledInputInfo;
}


class util::SampledInputInfo : public art::EDProducer {
public:
  explicit SampledInputInfo(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SampledInputInfo(SampledInputInfo const&) = delete;
  SampledInputInfo(SampledInputInfo&&) = delete;
  SampledInputInfo& operator=(SampledInputInfo const&) = delete;
  SampledInputInfo& operator=(SampledInputInfo&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;
  void endSubRun(art::SubRun &) override;

  void SetupEventIDTree();

private:

  // Declare member data here.

  art::InputTag fSamplingInputTag;

  bool          fCreatePOTInfo;
  std::vector<std::string>   fDatasetNamesForPOT;
  art::InputTag fPOTInfoTag;

  bool   fCreateEventIDTree;
  TTree* fEventIDTree;

  unsigned int fRun_new;
  unsigned int fSubrun_new;
  unsigned int fEvent_new;

  unsigned int fRun_old;
  unsigned int fSubrun_old;
  unsigned int fEvent_old;
  std::string  fDataset_old;

  bool fVerbose;

};



void util::SampledInputInfo::SetupEventIDTree()
{
  fEventIDTree->Branch("run_new",&fRun_new,"run_new/i");
  fEventIDTree->Branch("subrun_new",&fSubrun_new,"subrun_new/i");
  fEventIDTree->Branch("event_new",&fEvent_new,"event_new/i");

  fEventIDTree->Branch("run_old",&fRun_old,"run_old/i");
  fEventIDTree->Branch("subrun_old",&fSubrun_old,"subrun_old/i");
  fEventIDTree->Branch("event_old",&fEvent_old,"event_old/i");

  fEventIDTree->Branch("dataset_old",&fDataset_old);
}


util::SampledInputInfo::SampledInputInfo(fhicl::ParameterSet const& p)
  : EDProducer{p} ,
  fSamplingInputTag(p.get<art::InputTag>("SamplingInputTag")),
  fCreatePOTInfo(p.get<bool>("CreatePOTInfo")),
  fCreateEventIDTree(p.get<bool>("CreateEventIDTree")),
  fVerbose(p.get<bool>("Verbose",false))

  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  if(fCreatePOTInfo){
    this->produces<sumdata::POTSummary,art::InSubRun>();
    fDatasetNamesForPOT = p.get< std::vector<std::string> >("DatasetNamesForPOT");
    fPOTInfoTag = p.get<art::InputTag>("POTInfoTag");
  }
  
}

void util::SampledInputInfo::beginJob()
{

  if(fCreateEventIDTree) {
    art::ServiceHandle<art::TFileService> tfs;    
    fEventIDTree = tfs->make<TTree>("id_tree","Event ID Tree");
    SetupEventIDTree();
  }
}


void util::SampledInputInfo::produce(art::Event& e)
{
  
  if(fCreateEventIDTree){
    fRun_new = e.run();
    fSubrun_new = e.subRun();
    fEvent_new = e.event();

    auto const& sampled_info = *e.getValidHandle<art::SampledEventInfo>(fSamplingInputTag);
    auto const& oid = sampled_info.id;

    fRun_old = oid.run();
    fSubrun_old = oid.subRun();
    fEvent_old = oid.event();
    fDataset_old = sampled_info.dataset;

    if(fVerbose) std::cout << "sampled in " << fRun_old << " " << fEvent_old << " from " << fDataset_old << std::endl;
    
    fEventIDTree->Fill();
  }

}

void util::SampledInputInfo::endSubRun(art::SubRun & s)
{
  if(!fCreatePOTInfo) return;

  auto const id = s.id();  // SubRunID{1, 0};

  //get underlying subrun info
  auto const& sampled_info = *s.getValidHandle<art::SampledSubRunInfo>(fSamplingInputTag);

  art::InputTag const orig_potinfo_tag{fPOTInfoTag.label(),fPOTInfoTag.instance(), art::sampled_from(fPOTInfoTag.process())};
  auto const orig_potinfo = 
    s.getValidHandle< art::Sampled<sumdata::POTSummary> >(orig_potinfo_tag);

  //init new pot
  std::unique_ptr<sumdata::POTSummary> srpot_ptr(new sumdata::POTSummary());
  srpot_ptr->totpot = 0;
  srpot_ptr->totgoodpot = 0;
  srpot_ptr->totspills = 0;
  srpot_ptr->goodspills = 0;


  for(auto const& name : fDatasetNamesForPOT){
    auto const& pot_info = sampled_info.at(name);

    //sum over all the sampled subruns
    for(auto const& sid : pot_info.ids){
      auto const potsum_ptr = orig_potinfo->get(name,sid);
      srpot_ptr->totpot += potsum_ptr->totpot;
      srpot_ptr->totgoodpot += potsum_ptr->totgoodpot;
      srpot_ptr->totspills += potsum_ptr->totspills;
      srpot_ptr->goodspills += potsum_ptr->goodspills;
      
      if(fVerbose) std::cout << "\t\tPOT for subrun " << sid.run() << ":" << sid.subRun()
			     << " = " << potsum_ptr->totpot << std::endl;
    }
  }

  if(fVerbose) std::cout << "\t Total POT for subrun " << id.run() << ":" << id.subRun()
			 << " = " << srpot_ptr->totpot << std::endl;

  s.put(std::move(srpot_ptr));
  
}

DEFINE_ART_MODULE(util::SampledInputInfo)
