////////////////////////////////////////////////////////////////////////
// Class:       AnodeCathodeTrackFilter
// Module Type: filter
// File:        AnodeCathodeTrackFilter_module.cc
//
// Generated at Sun Aug 21 22:09:38 2016 by Wesley Ketchum using artmod
// from cetpkgsupport v1_10_02.
//
// wketchum@fnal.gov
//
// Algorithm for quickly identifying anode-to-cathode tracks based on a
// simple pattern-matching approach.
//
// Input:  vector<recob::Hit>
// Output: true if anode-cathode track candidate in event (false if not)
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

#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>
#include <iostream>

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "uboone/PatternFilter/PMAlgs/AnodeCathodePMAlg.h"

namespace pm {
  class AnodeCathodeTrackFilter;
}

class pm::AnodeCathodeTrackFilter : public art::EDFilter {
public:
  explicit AnodeCathodeTrackFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnodeCathodeTrackFilter(AnodeCathodeTrackFilter const &) = delete;
  AnodeCathodeTrackFilter(AnodeCathodeTrackFilter &&) = delete;
  AnodeCathodeTrackFilter & operator = (AnodeCathodeTrackFilter const &) = delete;
  AnodeCathodeTrackFilter & operator = (AnodeCathodeTrackFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) ;
  void beginJob() override;

private:

  art::InputTag     fHitLabel;
  float             fFractionMatchingThreshold;
  AnodeCathodePMAlg fAlg;
  bool              fVerbose;

};


pm::AnodeCathodeTrackFilter::AnodeCathodeTrackFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool pm::AnodeCathodeTrackFilter::filter(art::Event & e)
{

  auto const& hitVector = *e.getValidHandle< std::vector<recob::Hit> >(fHitLabel);
  float result;
  fAlg.RunPatternMatching(hitVector,result);

  if(fVerbose) std::cout << "\tDone matching. Result = " << result << std::endl;
  
  if(result > fFractionMatchingThreshold) return true;

  return false;
}

void pm::AnodeCathodeTrackFilter::reconfigure(fhicl::ParameterSet const & p)
{

  auto const* geo     = lar::providerFrom<geo::Geometry>();  
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();  
  fHitLabel = art::InputTag( p.get<std::string>("HitLabel") );
  fFractionMatchingThreshold = p.get<float>("FractionMatchingThreshold");

  float frac = p.get<fhicl::ParameterSet>("AnodeCathodPMAlg").get<float>("MatchFractionThreshold",10.0);
  if(fFractionMatchingThreshold > frac)
    {
      mf::LogWarning("AnodeCathodeTrackFilter")
	<< "FractionMatchingThreshold misconfigured:"
	<< "  it's greater than MatchFractionThreshold in AnodeCathodePMAlg."
	<< "\nSetting them equal or FractionMatchingThreshold to 0.0, so all events may pass.";
      if(frac>0.0 && frac<1.0)
	fFractionMatchingThreshold = frac;
      else
	fFractionMatchingThreshold = 0.0;
    }
  
  fAlg.Configure(p.get<fhicl::ParameterSet>("AnodeCathodPMAlg"),*geo,*detprop);

  fVerbose = p.get<bool>("Verbose",false);
}

void pm::AnodeCathodeTrackFilter::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pm::AnodeCathodeTrackFilter)
