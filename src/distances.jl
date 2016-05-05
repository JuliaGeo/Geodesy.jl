"""
    distance(a, b, [datum = wgs84])

The Cartesian distance between points `a` and `b`. Uses `datum` to perform
transformations as necessary.
"""
distance{T <: FixedVector}(a::T, b::T) = norm(a-b)

# Automatic transformations to ECEF
# Uses wgs84 datum as a default. In most cases, the datum choice will only make
# a small difference to the answer. Nevertheless, is this acceptable?
distance{T <: FixedVector}(a::T, b::T, datum) = distance(a, b)
distance(a::LLA, b, datum = wgs84) = distance(transform(ECEFfromLLA(datum), a), b, datum)
distance(a::ECEF, b::LLA, datum = wgs84) = distance(a, transform(ECEFfromLLA(datum), b), datum)
distance(a::LatLon, b, datum = wgs84) = distance(transform(ECEFfromLLA(datum), LLA(a)), b, datum)
distance(a::ECEF, b::LatLon, datum = wgs84) = distance(a, transform(ECEFfromLLA(datum), LLA(b)), datum)


# Also add geodesic distances here