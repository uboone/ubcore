////////////////////////////////////////////////////////////////////////
// Class:       HiggsPortalGen
// Plugin Type: producer (art v3_01_02)
// File:        HiggsPortalGen_module.cc
//
// Generated at Thu Mar 26 20:58:26 2020 by Andrew Mastbaum using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

/**
 * E1 px1 py1 pz1 E2 px2 py2 pz2 weight x_K y_K z_K px_K py_K pz_K L_enter L_exit
 *
 * (E1 , px1 , py1 , pz1): Four vector for the first scalar decay product (GeV)
 * (E2 , px2 , py2 , pz2): Four vector for the second scalar decay product (GeV)
 * weight: Event weight at listed mixing angle
 * (x_K , y_K , z_K): Position of kaon at time of decay to scalar (cm)
 * (px_K , py_K , pz_K): Momentum of kaon at time of decay to scalar (GeV)
 * L_enter, L_exit: Distance from kaon decay location to entry and exit from detector (cm)
 *
 * The first 9 columns are the most crucial for now. The last ones are useful
 * for rescaling the benchmark model to different mixing angles and for
 * debugging, which will come at a later stage in this process. The
 * coordinate system is such that the z-axis is the beamline. The scalar
 * decay positions should be to good-enough-for-now approximation evenly
 * distributed throughout the detector volume.
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larsim/EventGenerator/MARLEY/ActiveVolumeVertexSampler.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

namespace evgen {

class HiggsPortalGen : public art::EDProducer {
public:
  explicit HiggsPortalGen(fhicl::ParameterSet const& p);
  virtual ~HiggsPortalGen();

  HiggsPortalGen(HiggsPortalGen const&) = delete;
  HiggsPortalGen(HiggsPortalGen&&) = delete;
  HiggsPortalGen& operator=(HiggsPortalGen const&) = delete;
  HiggsPortalGen& operator=(HiggsPortalGen&&) = delete;

  void beginRun(art::Run & run) override;
  void produce(art::Event& e) override;

private:
  std::ifstream* fInputFile;  //!< Content of input file
  CLHEP::HepRandomEngine& fRandomEngine;  //!< RNG engine
  std::unique_ptr<evgen::ActiveVolumeVertexSampler> fVtxSampler;  //!< Vertex sampler
  int fDecayPDG;  //!< Decay particle PDG
  int fKPDG;  //!< Parent kaon PDG
};

HiggsPortalGen::HiggsPortalGen(fhicl::ParameterSet const& p)
    : art::EDProducer{p},
      fInputFile(nullptr),
      fRandomEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(
                    *this, "HepJamesRandom", "hpgen", p, "Seed")) {
  // Input file
  std::string inputFileName = p.get<std::string>("InputFileName");
  fInputFile = new std::ifstream(inputFileName.c_str());

  // FHiCL parameters
  fDecayPDG = p.get<int>("DecayPDG");
  fKPDG = p.get<int>("KaonPDG");
  auto const& vtxConfig = p.get<fhicl::ParameterSet>("VertexConfig", {});

  // Vertex sampler
  const auto& seedService = art::ServiceHandle<rndm::NuRandomService>();
  const auto& geomService = art::ServiceHandle<geo::Geometry const>();
  fVtxSampler = std::make_unique<evgen::ActiveVolumeVertexSampler>(
    vtxConfig, *seedService, *geomService, "HPGVtxSampler");

  produces<std::vector<simb::MCTruth> >();
  produces<sumdata::RunData, art::InRun>();
}


HiggsPortalGen::~HiggsPortalGen() {}


void HiggsPortalGen::beginRun(art::Run& run) {
  art::ServiceHandle<geo::Geometry> geo;
  std::unique_ptr<sumdata::RunData> runcol(
    new sumdata::RunData(geo->DetectorName()));
  run.put(std::move(runcol));
}


void HiggsPortalGen::produce(art::Event& e) {
  if (!fInputFile->good()) {
    throw cet::exception("HiggsPortalGen")
      << "input text file cannot be read.\n";
  }

  std::unique_ptr<std::vector<simb::MCTruth> > truths(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;
  truth.SetOrigin(simb::kSingleParticle);

  float E1, px1, py1, pz1,
        E2, px2, py2, pz2,
        weight,
        xK, yK, zK, pxK, pyK, pzK,
        lEnter, lExit;

  std::string line;
  std::getline(*fInputFile, line);

  if (line.empty()) {
    e.put(std::move(truths));
    return;
  }

  std::istringstream iline(line);
  iline >> E1 >> px1 >> py1 >> pz1
        >> E2 >> px2 >> py2 >> pz2
        >> weight
        >> xK >> yK >> zK
        >> pxK >> pyK >> pzK
        >> lEnter >>  lExit;

  // Parent kaon
  /*
  TLorentzVector posK(xK, yK, zK, 0);
  TVector3 pK(pxK, pyK, pzK);
  const TParticlePDG* kdef = TDatabasePDG::Instance()->GetParticle(fKPDG);
  float mK = kdef ? kdef->Mass() : 0;
  float EK = sqrt(pK.Mag() * pK.Mag() + mK * mK);
  TLorentzVector momK(pK, EK);

  simb::MCParticle partK(-1, fDecayPDG, "", -1, mK, 0);
  partK.AddTrajectoryPoint(posK, momK);
  truth.Add(partK);
*/
  // Decay products
  art::ServiceHandle<geo::Geometry> geo;
  TLorentzVector pos = fVtxSampler->sample_vertex_pos(*geo);

  TLorentzVector mom1(px1, py1, pz1, E1);
  simb::MCParticle part1(-1, fDecayPDG, "primary", -1);
  part1.AddTrajectoryPoint(pos, mom1);
  truth.Add(part1);

  TLorentzVector mom2(px2, py2, pz2, E2);
  simb::MCParticle part2(-2, -fDecayPDG, "primary", -2);
  part2.AddTrajectoryPoint(pos, mom2);
  truth.Add(part2);

  std::cout << truth;

  truths->push_back(truth);

  e.put(std::move(truths));
}

DEFINE_ART_MODULE(HiggsPortalGen)

}  // namespace evgen

