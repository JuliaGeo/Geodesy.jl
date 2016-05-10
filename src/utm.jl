# Utilities to convert between UTM+zone and Mercator projections

# bound theta in degs, so on output -pi <= theta < pi
bound_thetad(theta) = theta - floor((theta+180) / (360)) * 360

"""
    utm_zone(lat, lon)
    utm_zone(ll::LatLon)
    utm_zone(lla::LLA)
    utm_zone(ecef::ECEF, datum)

Find the UTM zone and hemisphere for the given latitude and longitude (or world
point). Zone 0 corresponds to the UPS regions.
"""
function utm_zone{T}(lat::T, lon::T)
    if lat > 84
        return (0, true)
    elseif lat < -80
        return (0, false)
    end

    # int versions
    ilat = floor(Int64, bound_thetad(lat))
    ilon = floor(Int64, bound_thetad(lon))

    # get the latitude band
    band = max(-10, min(9,  fld((ilat + 80), 8) - 10))

    # and check for weird ones
    zone = fld((ilon + 186), 6);
    if ((band == 7) && (zone == 31) && (ilon >= 3))
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42))
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    # TODO: Svalbard (31X-37X)

    return (zone, lat >= 0)

end

utm_zone(ll::LatLon) = utm_zone(ll.lat, ll.lon)
utm_zone(lla::LLA) = utm_zone(lla.lat, lla.lon)
utm_zone(ecef::ECEF, datum) = utm_zone(transform(LLAfromECEF(datum), ecef))

function utm_meridian(zone::Integer)
    if zone < 1 || zone > 60
        error("UTM zone $zone not defined")
    end
    return 6 * zone - 183
end
