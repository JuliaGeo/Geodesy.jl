
###################
### Point Types ###
###################

### Point in Latitude-Longitude-Altitude (LLA) coordinates
immutable LLA{T <: Datum}
    lat::Float64
    lon::Float64
    alt::Float64
end
Base.call{T}(::Type{LLA{T}}, lat::Real, lon::Real) = LLA{T}(lat, lon, 0.0)
LLA(args...) = LLA{WGS84}(args...)

ellipsoid{T}(::Type{LLA{T}}) = ellipsoid(T)

### Latitude-Longitude (LL) coordinates
immutable LL{T <: Datum}
    lat::Float64
    lon::Float64
end
LL(args...) = LL{WGS84}(args...)

ellipsoid{T}(::Type{LL{T}}) = ellipsoid(T)

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

### XYZ
# Helper for creating other point types in generic code
# e.g. myfunc{T <: Union(ENU, LLA)}(...) = (x, y = ...; T(XY(x, y)))
type XYZ
    x::Float64
    y::Float64
    z::Float64
end
XY(x, y) = XYZ(x, y, 0.0)

Base.call{T}(::Type{LL{T}}, xyz::XYZ) = LL{T}(xyz.y, xyz.x)
Base.call{T}(::Type{LLA{T}}, xyz::XYZ) = LLA{T}(xyz.y, xyz.x, xyz.z)
ENU(xyz::XYZ) = ENU(xyz.x, xyz.y, xyz.z)

### get*
# Point translators
getX(ll::LL) = ll.lon
getY(ll::LL) = ll.lat

getX(lla::LLA) = lla.lon
getY(lla::LLA) = lla.lat
getZ(lla::LLA) = lla.alt

getX(enu::ENU) = enu.east
getY(enu::ENU) = enu.north
getZ(enu::ENU) = enu.up

### distance
# Point translators
distance(a::ENU, b::ENU) = distance(a.east, a.north, a.up,
                                    b.east, b.north, b.up)

distance(a::ECEF, b::ECEF) = distance(a.x, a.y, a.z,
                                      b.x, b.y, b.z)

function distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
end
