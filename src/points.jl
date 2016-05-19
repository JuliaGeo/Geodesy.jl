###################
### Point Types ###
###################

"""
    LLA(lat, lon, alt = 0.0)
    LLA(lat = ϕ, lon = Θ, alt = h)

Latitude, longitude, and alititude co-ordinates. *Note:* assumes degrees not radians
"""
immutable LLA{T}
    lat::T
    lon::T
    alt::T
end
LLA(;lat=NaN,lon=NaN,alt=0.0) = LLA(lat,lon,alt) # Constructor that is idependent of storage order
LLA(lat, lon) = LLA(lat, lon, 0.0)
Base.show(io::IO, lla::LLA) = print(io, "LLA(lat=$(lla.lat)°, lon=$(lla.lon)°, alt=$(lla.alt))") # Maybe show to nearest mm, or 6 decimel places?
Base.isapprox(lla1::LLA, lla2::LLA; atol = 1e-6, kwargs...) = isapprox(lla1.lat, lla2.lat; atol = 180*atol/6.371e6, kwargs...) & isapprox(lla1.lon, lla2.lon; atol = 180*atol/6.371e6, kwargs...) & isapprox(lla1.alt, lla2.alt; atol = atol, kwargs...) # atol in metres (1μm)

# For radians, we could export ° = pi/180, then output in degrees... ??

"""
    LatLon(lat, lon)
    LatLon(lat = ϕ, lon = Θ)

Latitude and longitude co-ordinates. *Note:* assumes degrees not radians
"""
immutable LatLon{T}
    lat::T
    lon::T
end
LatLon(;lat=NaN,lon=NaN) = LatLon(lat,lon) # Constructor that is idependent of storage order
LatLon(lla::LLA) = LatLon(lla.lat, lla.lon)
function LatLon(x, datum) #?
    lla = LLA(x, datum)
    return LatLon(lla.lat, lla.lon)
end
LLA(ll::LatLon) = LLA(ll.lat, ll.lon)
Base.show(io::IO, ll::LatLon) = print(io, "LatLon(lat=$(ll.lat)°, lon=$(ll.lon)°)") # Maybe show to nearest mm, or 6 decimel places?
Base.isapprox(ll1::LatLon, ll2::LatLon; atol = 1e-6, kwargs...) = isapprox(ll1.lat, ll2.lat; atol = 180*atol/6.371e6, kwargs...) & isapprox(ll1.lon, ll2.lon; atol = 180*atol/6.371e6, kwargs...) # atol in metres (1μm)


"""
    ECEF(x, y, z)

Earth-Centered-Earth-Fixed (ECEF) coordinates. A global Cartesian coordinate
system rotating with the Earth.
"""
immutable ECEF{T} <: FixedVectorNoTuple{3,Float64}
    x::T
    y::T
    z::T
end
Base.show(io::IO, ecef::ECEF) = print(io, "ECEF($(ecef.x), $(ecef.y), $(ecef.z))")


"""
    ENU(e, n, u = 0.0)

East-North-Up (ENU) coordinates. A local Cartesian coordinate system, linearized about a reference point.
"""
immutable ENU{T} <: FixedVectorNoTuple{3,Float64}
    e::T
    n::T
    u::T
end
ENU(x, y) = ENU(x, y, 0.0)
Base.show(io::IO, enu::ENU) = print(io, "ENU($(enu.e), $(enu.n), $(enu.u))")


"""
    UTM(x, y, z = 0.0)

Universal transverse Mercator (UTM) coordinates. Common projection type for
world points. Zone not included in coordinates - it is a parameterized in the
relavant transformations `UTMfromLLA` and `LLAfromUTM` (see also the `UTMZ` type).

This type may be used to parameterize UPS coordinates (Universal Polar
Stereographic) to accurately represent the polar regions, in zone "0".
"""
immutable UTM{T}
    x::T
    y::T
    z::T
end
UTM(x, y) = UTM(x, y, 0.0)
Base.isapprox(utm1::UTM, utm2::UTM; atol = 1e-6, kwargs...) = isapprox(utm1.x, utm2.x; atol = atol, kwargs...) & isapprox(utm1.y, utm2.y; atol = atol, kwargs...) & isapprox(utm1.z, utm2.z; atol = atol, kwargs...) # atol in metres (1μm)
Base.show(io::IO, utm::UTM) = print(io, "UTM($(utm.x), $(utm.y), $(utm.z))")


"""
    UTMZ(x, y, z = 0.0, zone::Integer, hemisphere::Bool)

Universal transverse Mercator (UTM) coordinates with zone number. Common
projection type for world points. The UTM zone is included in coordinates
(see also the `UTM` type).

This type may be used to parameterize UPS coordinates (Universal Polar
Stereographic) to accurately represent the polar regions, in zone "0".
"""
immutable UTMZ{T}
    x::T
    y::T
    z::T
    zone::UInt8
    hemisphere::Bool # true = north, false = south
end
UTMZ{T}(x::T, y::T, zone::Integer, hemisphere::Bool) = UTMZ(x, y, 0.0, UInt8(zone), hemisphere)
UTMZ(x, y, z, zone::Integer, hemisphere::Bool) = UTMZ(x, y, z, UInt8(zone), hemisphere)
UTMZ(utm::UTM, zone::Integer, hemisphere::Bool) = UTMZ(utm.x, utm.y, utm.z, UInt8(zone), hemisphere)
UTM(utmz::UTMZ) = UTM(utmz.x, utmz.y, utmz.z)
Base.isapprox(utm1::UTMZ, utm2::UTMZ; atol = 1e-6, kwargs...) = (utm1.zone == utm2.zone) & (utm1.hemisphere == utm2.hemisphere) & isapprox(utm1.x, utm2.x; atol = atol, kwargs...) & isapprox(utm1.y, utm2.y; atol = atol, kwargs...) & isapprox(utm1.z, utm2.z; atol = atol, kwargs...) # atol in metres (1μm)
Base.show(io::IO, utm::UTMZ) = print(io, "UTMZ($(utm.x), $(utm.y), $(utm.z), zone=$(utm.zone == 0 ? "polar" : Int(utm.zone)) ($(utm.hemisphere ? "north" : "south")))")
