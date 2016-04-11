module Geodesy

using FixedSizeArrays
using Compat

export
    # Points
    ECEF,
    ENU,
    LLA,

    # Other types
    Ellipsoid,

    # Constants
    WGS84, wgs84,
    OSGB36, osgb36,
    NAD27, nad27,

    # Methods
    distance


include("points.jl")
include("datums.jl")
include("conversion.jl")

end # module Geodesy
