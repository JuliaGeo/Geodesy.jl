module Geodesy

using CoordinateTransformations
using FixedSizeArrays
using Compat

import CoordinateTransformations.transform,
       CoordinateTransformations.transform,
       CoordinateTransformations.transform_deriv_params,
       CoordinateTransformations.compose,
       CoordinateTransformations.∘

export
    # Points
    ECEF,
    ENU,
    LLA,
    LatLon,
    UTM,
    UTMZ,

    # Other types
    Ellipsoid,

    # Constants
    WGS84, wgs84,
    OSGB36, osgb36,
    NAD27, nad27,
    GRS80, grs80,

    # Methods
    distance,

    # transformation methods
    transform, transform, transform_deriv_params, compose, ∘,
    ECEFfromLLA, LLAfromECEF, ENUfromECEF, ECEFfromENU, ENUfromLLA, LLAfromENU,
    UTMfromLLA, LLAfromUTM, UTMfromECEF, ECEFfromUTM,
    UTMZfromLLA, LLAfromUTMZ, UTMZfromECEF, ECEFfromUTMZ,

    # UTM helpers
    utm_zone



include("points.jl")
include("datums.jl")
include("transformations.jl")
include("conversion.jl")
include("distances.jl")
include("transverse_mercator.jl")
include("utm.jl")

end # module Geodesy
