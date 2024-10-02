////////////////////////////////////////////////////////////////////////////////
/// \file UBOONEGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#include "ubcore/Geometry/UBooNEWireReadout.h"

#include "larcore/Geometry/Geometry.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace uboone {

  UBooNEWireReadout::UBooNEWireReadout(fhicl::ParameterSet const& pset)
    : fWireReadout{ art::ServiceHandle<geo::Geometry>().get(), pset}
  {
    if ( art::ServiceHandle<geo::Geometry>()->DetectorName().find("microboone") == std::string::npos ) {
      std::cout << __PRETTY_FUNCTION__ << ": WARNING USING WIRE READOUT WITH NON-MICROBOONE GEO!" << std::endl;
    }
  }

}
