module Geodesy

using CoordinateTransformations
using FixedSizeArrays
using Compat

import CoordinateTransformations.transform_deriv,
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
    Ellipsoid, ellipsoid,

    # Constants
    WGS84, wgs84,
    OSGB36, osgb36,
    NAD27, nad27,
    GRS80, grs80,

    # Methods
    distance,

    # transformation methods
    transform_deriv, transform_deriv_params, compose, ∘,
    ECEFfromLLA, LLAfromECEF, ENUfromECEF, ECEFfromENU, ENUfromLLA, LLAfromENU,
    UTMfromLLA, LLAfromUTM, UTMfromECEF, ECEFfromUTM, ENUfromUTM, UTMfromENU,
    UTMZfromLLA, LLAfromUTMZ, UTMZfromECEF, ECEFfromUTMZ, ENUfromUTMZ, UTMZfromENU,
    UTMZfromUTM, UTMfromUTMZ,

    # UTM helpers
    utm_zone

include("points.jl")
include("datums.jl")
include("transverse_mercator.jl")
include("polar_stereographic.jl")
include("transformations.jl")
include("conversion.jl")
include("distances.jl")
include("utm.jl")

end # module Geodesy
