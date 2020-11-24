using LinearAlgebra: norm

"""
    euclidean_distance(a, b, [datum = wgs84])

The straight-line distance between points `a` and `b`. Uses `datum` to perform
transformations as necessary, and unlike other *Geodesy* functions, this
defaults to use the common WGS-84 datum for convenience.
"""
euclidean_distance(a::T, b::T) where {T <: AbstractVector} = norm(a-b)

# Automatic transformations to ECEF
# Uses wgs84 datum as a default. In most cases, the datum choice will only make
# a small difference to the answer. Nevertheless, is this acceptable?
euclidean_distance(a::T, b::T, datum) where {T <: AbstractVector} = euclidean_distance(a, b)
euclidean_distance(a::LLA, b, datum = wgs84) = euclidean_distance(ECEFfromLLA(datum)(a), b, datum)
euclidean_distance(a::ECEF, b::LLA, datum = wgs84) = euclidean_distance(a, ECEFfromLLA(datum)(b), datum)
euclidean_distance(a::LatLon, b, datum = wgs84) = euclidean_distance(ECEFfromLLA(datum)(LLA(a)), b, datum)
euclidean_distance(a::ECEF, b::LatLon, datum = wgs84) = euclidean_distance(a, ECEFfromLLA(datum)(LLA(b)), datum)
euclidean_distance(a::UTMZ, b, datum = wgs84) = euclidean_distance(ECEFfromUTMZ(datum)(a), b, datum)
euclidean_distance(a::ECEF, b::UTMZ, datum = wgs84) = euclidean_distance(a, ECEFfromUTMZ(datum)(b), datum)

"""
    euclidean_distance(utm1, utm2, zone, isnorth, [datum = wgs84])
    euclidean_distance(a, utm2, zone, isnorth, [datum = wgs84])
    euclidean_distance(utm1, b, zone, isnorth, [datum = wgs84])

If one or both points are UTM, we need the zone (and particularly the hemisphere,
isnorth = true/false) to determine the Euclidean distance.
"""
euclidean_distance(a::UTM, b::UTM, zone::Integer, isnorth::Bool, datum = wgs84) = euclidean_distance(ECEFfromUTM(zone, isnorth, datum)(a), ECEFfromUTM(zone, isnorth, datum)(b), datum)
euclidean_distance(a, b::UTM, zone::Integer, isnorth::Bool, datum = wgs84) = euclidean_distance(a, ECEFfromUTM(zone, isnorth, datum)(b), datum)
euclidean_distance(a::UTM, b, zone::Integer, isnorth::Bool, datum = wgs84) = euclidean_distance(ECEFfromUTM(zone, isnorth, datum)(a), b, datum)

