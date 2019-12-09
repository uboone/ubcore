////////////////////////////////////////////////////////////////////////
// Class:       BurstNoiseMetricsFilter
// Module Type: filter
// File:        BurstNoiseMetricsFilter_module.cc
//
// Based on a producer module BurstNoiseMetrics_module.cc written by Ivan Caro
// Generated at Tue Nov 28 11:16:47 2019 by Elena Gramellini using artmod
// from cetpkgsupport v1_13_00.
// Filter is set by the nf and fft filter threshold values. If event has
// less than these values then keep the event.
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
#include "art/Framework/Services/Optional/TFileService.h"

// Data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
//#include "BNMS.h"

// ROOT includes
#include "TVector3.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtuple.h"

// Standard Library
#include <memory>
#include <iostream>
#include <utility>

class BurstNoiseMetricsFilter;
class BurstNoiseMetricsFilter : public art::EDFilter {
public:
 
  explicit BurstNoiseMetricsFilter(fhicl::ParameterSet const & p);
  BurstNoiseMetricsFilter(BurstNoiseMetricsFilter const &) = delete;
  BurstNoiseMetricsFilter(BurstNoiseMetricsFilter &&) = delete;
  BurstNoiseMetricsFilter & operator = (BurstNoiseMetricsFilter const &) = delete;
  BurstNoiseMetricsFilter & operator = (BurstNoiseMetricsFilter &&) = delete;
  // Required functions.
  bool filter(art::Event &e) override;
  void beginJob() override; 

private:
  int median_func(std::vector<double> Chans); 
  // Declare member data here.
  std::string fFlashLabel;
  std::string fRawDigitLabel;
  int    fWin_start;
  int    fWin_end;
  int    fIntWinFFTsum;
  int    fIntWinNF;
  bool   fuseAsFilter;
  double fFFT_FiltThreshold;
  double fnf_FiltThreshold;
  TH1D *h_time;
  TH1D *U_uberwf;
  TH1D *FFT_U_uberwf;
  TH1D *h_nf;
  TH1D *h_fftsum;

  bool _debug{false};
};
void BurstNoiseMetricsFilter::BurstNoiseMetricsFilter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  
  // Histogram storing the event flashes
  h_time       = tfs->make<TH1D>("h_time","Histogram;Counts;Time", fWin_end - fWin_start, fWin_start, fWin_end);
  
  // Histogram storing the U plane Uberwaveform
  U_uberwf     = tfs->make<TH1D>("U_uberwf", "Event UberWF for U plane;Time Tick; ADC Value",9594, -0.5, 9593.5);
  
  // Histogram used to store the FFT of U plane uber waveform
  FFT_U_uberwf = tfs->make<TH1D>("FFT_U_uberwf", "FFT of U_UberWF; ADCvals; Frequencies",9594, -0.5, 9593.5);

  // Purity Burst Metric Histogram
  h_nf = tfs->make<TH1D>("h_nf", "Purtity Burst Metric; nf; Entries",100, 0, 50);

  // Cathode Burst Metric Histogram
  h_fftsum = tfs->make<TH1D>("h_fftsum", "Cathode Burst Metric; fftsum; Entries",100, 0, 150000);

}

BurstNoiseMetricsFilter::BurstNoiseMetricsFilter(fhicl::ParameterSet const & p) {

  produces< double >("CBmetric" ); // Cathode Burst Metric
  produces< int >   ("PMBmetric"); // Purity Monitor Burst Metric
  
  fFlashLabel        = p.get<std::string>("FlashLabel"          );
  fRawDigitLabel     = p.get<std::string>("RawDigitLabel"       );
  fWin_start         = p.get<int>        ("Win_start"           );
  fWin_end           = p.get<int>        ("Win_end"             );
  fIntWinFFTsum      = p.get<int>        ("IntWinFFTsum"        );
  fIntWinNF          = p.get<int>        ("IntWinNF"            );
  fuseAsFilter       = p.get<bool>       ("useAsFilter"         );
  fFFT_FiltThreshold = p.get<double>     ("FFT_FiltThreshold"   );
  fnf_FiltThreshold  = p.get<double>     ("nf_FiltThreshold"   );
  _debug             = p.get<bool>       ("Debug"               );
}

bool BurstNoiseMetricsFilter::filter(art::Event & e) {

  // Initial parameters. Some have to be parametrized in the .fcl file 
  double fftsum;
  int nf;
  std::vector<double> U_UberADCvals;
  U_UberADCvals.clear();
  double U_median = 0;
  std::unique_ptr< double > unqptr_fftsum;
  std::unique_ptr< int > unqptr_nf;
  

  // Locate the flash objects within artroot file
  auto const& opflash_handle = e.getValidHandle< std::vector < recob::OpFlash > > (fFlashLabel);
  auto const& opflash_vec(*opflash_handle);
 
  // Locate the raw digit objects within the artroot file
  auto const& rawdigit_handle = e.getValidHandle< std::vector < raw::RawDigit > > (fRawDigitLabel);
  auto const& allrawdigits_vec(*rawdigit_handle);
  
  // ------------------ unqptr_nf00 Calculations -------------------------------
  // Store any flash within the time window -1600 to 3200 ms
  for( auto const& flash : opflash_vec){
    if (flash.Time() < fWin_start || flash.Time() > fWin_end) continue;
    h_time -> Fill(flash.Time());      
  }
  
  std::vector<int> vec_getmax;
  for (int k = fIntWinNF - 1; k < fWin_end - fWin_start; k++ ){
     vec_getmax.push_back(h_time -> Integral(k - fIntWinNF - 1, k));
  }
  nf = *std::max_element(vec_getmax.begin(), vec_getmax.end());
  
  //------------------------ FFTSUM2 Calculations ------------------------------
  double UChanADCval_temp = 0;
  U_median = 0.;
  
  
  for (int k = 0; k < allrawdigits_vec.at(1).Samples(); k++){
    
    // UberWaveform Calculations
    for (size_t i_ar = 0, size_allrawdigits = rawdigit_handle->size(); i_ar != size_allrawdigits; ++i_ar){
      
      // Collection plane only
      if (allrawdigits_vec.at(i_ar).Channel() < 2400){
        UChanADCval_temp += allrawdigits_vec.at(i_ar).ADC(k);
      } // if channel less than 2400
    
    } // i_are loop

    U_UberADCvals.push_back(UChanADCval_temp);
    UChanADCval_temp = 0;
  
  } // k loop

  U_median = median_func(U_UberADCvals);
  if (_debug) std::cout<<"U_median "<<U_median<<"\n";
  
  for (unsigned j = 0; j < U_UberADCvals.size(); j++){
    U_uberwf -> SetBinContent(j+1,( U_UberADCvals[j] - U_median)/pow(10.,3.));
  }

  
  //------------------------- UWF Calculations ---------------------------------
  FFT_U_uberwf = (TH1D*) U_uberwf -> FFT(FFT_U_uberwf, "MAG");
  fftsum = FFT_U_uberwf -> Integral(0, fIntWinFFTsum); 
  
  if (_debug) std::cout<<"******** nf, fftsum "<<nf<<","<<fftsum<<"\n";
  if ( fftsum > fFFT_FiltThreshold || nf > fnf_FiltThreshold) std::cout<<"BURST Event Detected: "<<e.run()<<" "<<e.subRun()<<" "<<e.event()<<"\n";
 
  h_nf ->Fill(nf);
  h_fftsum->Fill(fftsum); 
 
  unqptr_nf = std::unique_ptr< int >(new int (nf));
  unqptr_fftsum = std::unique_ptr< double >(new double (fftsum));
   
  e.put(std::move(unqptr_nf),      "PMBmetric"  );
  e.put(std::move(unqptr_fftsum),  "CBmetric"   );
  
  // Reset and clear
  h_time->Reset("ICESM");
  U_uberwf->Reset("ICESM");
  FFT_U_uberwf->Reset("ICESM");
  U_UberADCvals.clear();
  vec_getmax.clear();
  
  // Apply the filter condition
  if ( !fuseAsFilter )           return true;
  if ( fftsum < fFFT_FiltThreshold && nf < fnf_FiltThreshold) return true;
  
  return false;

}

//---------------------------- Functions ---------------------------------------
int BurstNoiseMetricsFilter::median_func(std::vector<double> Chans){
  size_t size = Chans.size();
  int median;
  sort(Chans.begin(), Chans.end());
  if (size  % 2 == 0){
     median = (Chans[size / 2 - 1] + Chans[size / 2]) / 2;
  }else{
    median = Chans[size / 2];
  }
  return median;
}

DEFINE_ART_MODULE(BurstNoiseMetricsFilter)
