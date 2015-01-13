
###################
### Point Types ###
###################

### Point in Latitude-Longitude-Altitude (LLA) coordinates
# Used to store node data in OpenStreetMap XML files
immutable LLA
    lat::Float64
    lon::Float64
    alt::Float64
end
LLA(lat, lon) = LLA(lat, lon, 0.0)

### Point in Earth-Centered-Earth-Fixed (ECEF) coordinates
# Global cartesian coordinate system rotating with the Earth
immutable ECEF
    x::Float64
    y::Float64
    z::Float64
end

### Point in East-North-Up (ENU) coordinates
# Local cartesian coordinate system
# Linearized about a reference point
immutable ENU
    east::Float64
    north::Float64
    up::Float64
end
ENU(x, y) = ENU(x, y, 0.0)

### World Geodetic Coordinate System of 1984 (WGS 84)
# Standardized coordinate system for Earth
# Global ellipsoidal reference surface
immutable WGS
    a::Float64        # Semi-major axis
    b::Float64        # Semi-minor axis
    e::Float64        # Eccentricity
    e_prime::Float64  # Second eccentricity
    #N # Was undef & not used anywhere in OpenStreetMap.jl
end

const WGS84 = let a = 6378137.0,
                  b = 6356752.31424518,
                  e = sqrt((a*a - b*b) / (a*a)),
                  e_prime = sqrt((a*a - b*b) / (b*b))
    WGS(a, b, e, e_prime)
end

### Helper for creating other point types
type XYZ
    x::Float64
    y::Float64
    z::Float64
end
XY(x, y) = XYZ(x, y, 0.0)

LLA(xyz::XYZ) = LLA(xyz.y, xyz.x, xyz.z)
ENU(xyz::XYZ) = ENU(xyz.x, xyz.y, xyz.z)

### Point translators
getX(lla::LLA) = lla.lon
getY(lla::LLA) = lla.lat
getZ(lla::LLA) = lla.alt

getX(enu::ENU) = enu.east
getY(enu::ENU) = enu.north
getZ(enu::ENU) = enu.up

####################
### Bounds Types ###
####################

type Bounds{T <: Union(LLA, ENU)}
    min_y::Float64
    max_y::Float64
    min_x::Float64
    max_x::Float64
end
function Bounds(min_lat, max_lat, min_lon, max_lon)
    if !(-90 <= min_lat <= max_lat <= 90 && -180 <= min_lon <= max_lon <= 180)
        throw(ArgumentError("Bounds out of range of LLA coordinate system. " *
                            "Perhaps you're looking for Bounds{ENU}(...)"))
    end
    Bounds{LLA}(min_lat, max_lat, min_lon, max_lon)
end
