module Geodesy

export
    # Points
    ECEF,
    ENU,
    LLA,

    # Other types
    Bounds,
    Ellipsoid,

    # Constants
    WGS84,
    OSGB36,
    NAD27,

    # Methods
    center,
    distance,
    getX,
    getY,
    getZ,
    inBounds

include("points.jl")
include("bounds.jl")
include("conversion.jl")

end # module Geodesy
