////////////////////////////////////////////////////////////////////////////////
/// \file UBOONEGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> ubcore/Geometry
#include "ubcore/Geometry/UBooNEGeometryHelper.h"

#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h" // larcore. geo::GeometryData_t

// Migration note:
// Geometry --> ubcore/Geometry for the two below
#include "ubcore/Geometry/ChannelMapUBooNEAlg.h"

namespace uboone {

  UBooNEGeometryHelper::UBooNEGeometryHelper(fhicl::ParameterSet const& pset)
    : fPset(pset)
  {}

  std::unique_ptr<geo::ChannelMapAlg>
  UBooNEGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                                                 std::string const& detectorName) const
  {
    if ( detectorName.find("microboone") == std::string::npos ) {
      std::cout << __PRETTY_FUNCTION__ << ": WARNING USING CHANNEL MAP ALG WITH NON-MICROBOONE GEO!" << std::endl;
    }
    return std::make_unique<geo::ChannelMapUBooNEAlg>(fPset, sortingParameters);
  }

}
