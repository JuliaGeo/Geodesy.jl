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
    ECEFfromLLA, LLAfromECEF, ENUfromECEF, ECEFfromENU, ENUfromLLA, LLAfromENU



include("points.jl")
include("datums.jl")
include("transformations.jl")
include("conversion.jl")
include("distances.jl")


end # module Geodesy
