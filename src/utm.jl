# Utilities to convert between UTM+zone and Mercator projections

# bound theta in degs, so on output -pi <= theta < pi
bound_thetad(theta) = theta - floor((theta+180) / (360)) * 360

"""
    (zone, isnorth) = utm_zone(lat, lon)
    (zone, isnorth) = utm_zone(ll::LatLon)
    (zone, isnorth) = utm_zone(lla::LLA)
    (zone, isnorth) = utm_zone(ecef::ECEF, datum)

Find the UTM zone and hemisphere (`isnorth = true` or `false`) for the given
latitude and longitude (or world point), including the special rules for Norway
and Svalbard. Zone 0 corresponds to the poles, using the UPS regions.
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
    zone = fld((ilon + 186), 6)
    if ((band == 7) && (zone == 31) && (ilon >= 3)) # Norway
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42)) # Svalbard
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    return (zone, lat >= 0)
end

utm_zone(ll::LatLon) = utm_zone(ll.lat, ll.lon)
utm_zone(lla::LLA) = utm_zone(lla.lat, lla.lon)
utm_zone(ecef::ECEF, datum) = utm_zone(LLAfromECEF(datum)(ecef))

"""
    utm_meridian(zone)

Central meridian of the given UTM `zone` (note - does not include the
conventiional 500km false easting offset).
"""
function utm_meridian(zone::Integer)
    if zone < 1 || zone > 60
        error("UTM zone $zone not defined")
    end
    return 6 * zone - 183
end
