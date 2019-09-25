using LinearAlgebra: norm

"""
    distance(a, b, [datum = wgs84])

The Cartesian distance between points `a` and `b`. Uses `datum` to perform
transformations as necessary, and unlike other *Geodesy* functions, this defaults
to use the common WGS-84 datum for convenience.
"""
distance(a::T, b::T) where {T <: AbstractVector} = norm(a-b)

# Automatic transformations to ECEF
# Uses wgs84 datum as a default. In most cases, the datum choice will only make
# a small difference to the answer. Nevertheless, is this acceptable?
distance(a::T, b::T, datum) where {T <: AbstractVector} = distance(a, b)
distance(a::LLA, b, datum = wgs84) = distance(ECEFfromLLA(datum)(a), b, datum)
distance(a::ECEF, b::LLA, datum = wgs84) = distance(a, ECEFfromLLA(datum)(b), datum)
distance(a::LatLon, b, datum = wgs84) = distance(ECEFfromLLA(datum)(LLA(a)), b, datum)
distance(a::ECEF, b::LatLon, datum = wgs84) = distance(a, ECEFfromLLA(datum)(LLA(b)), datum)
distance(a::UTMZ, b, datum = wgs84) = distance(ECEFfromUTMZ(datum)(a), b, datum)
distance(a::ECEF, b::UTMZ, datum = wgs84) = distance(a, ECEFfromUTMZ(datum)(b), datum)

"""
    distance(utm1, utm2, zone, isnorth, [datum = wgs84])
    distance(a, utm2, zone, isnorth, [datum = wgs84])
    distance(utm1, b, zone, isnorth, [datum = wgs84])

If one or both points are UTM, we need the zone (and particularly the hemisphere,
isnorth = true/false) to determine the Cartesian distance.
"""
distance(a::UTM, b::UTM, zone::Integer, isnorth::Bool, datum = wgs84) = distance(ECEFfromUTM(zone, isnorth, datum)(a), ECEFfromUTM(zone, isnorth, datum)(b), datum)
distance(a, b::UTM, zone::Integer, isnorth::Bool, datum = wgs84) = distance(a, ECEFfromUTM(zone, isnorth, datum)(b), datum)
distance(a::UTM, b, zone::Integer, isnorth::Bool, datum = wgs84) = distance(ECEFfromUTM(zone, isnorth, datum)(a), b, datum)

import GeographicLib

export GreatCircle

"""
    GreatCircle

A `GreatCircle` object stores a cache of values needed to quickly compute
great circles.
"""
struct GreatCircle
    geod::GeographicLib.Geodesic
    datum::Datum
end

"""
    GreatCircle(datum) -> gc

Construct a `GreatCircle` `gc` from a datum.

The object `gc` can then be used as a function to calculate great circle properties
"""
GreatCircle(datum::Datum) =
    (el = ellipsoid(datum); GreatCircle(GeographicLib.Geodesic(el.a, el.f), datum))


(g::GreatCircle)(a::Union{LatLon,LLA}, b::Union{LatLon,LLA}) =
    GeographicLib.inverse(g.geod, a.lon, a.lat, b.lon, b.lat)

function (g::GreatCircle)(p::Union{LatLon,LLA}; azi, dist=nothing, angle=nothing)
    sum(x -> x === nothing, (dist, angle)) == 1 ||
        throw(ArgumentError("must specify one and only one of dist or angle"))
    if dist !== nothing
        GeographicLib.forward(g.geod, p.lon, p.lat, azi, dist)
    elseif angle !== nothing
        GeographicLib.forward_deg(g.geod, p.lon, p.lat, azi, angle)
    end
end

(g::GreatCircle)(a; kwargs...) = g(LLA(a, g.datum); kwargs...)
