"""
    WebMercatorfromLLA(datum=wgs84)

Convert from LLA to Web Mercator / Pseudo Mercator, following the convention of
proj, which uses a scaling factor of the semi-major axis of the ellipsoid
https://proj.org/operations/projections/webmerc.html.

!!! warning
    For web mapping applications this projection is ubiquitous due to its
    simplicity of implementation, but this simplicity gives rise to poor
    mathematical properties: it doesn't preserve angles away from the equator.
    Other projections should be preferred if possible, and especially when
    measurement is important. See
    https://earth-info.nga.mil/GandG/wgs84/web_mercator/(U)%20NGA_SIG_0011_1.0.0_WEBMERC.pdf
    for an extended discussion.
"""
struct WebMercatorfromLLA <: Transformation
    el::Ellipsoid
end

WebMercatorfromLLA(d::Datum=wgs84) = WebMercatorfromLLA(ellipsoid(d))

"""
    LLAfromWebMercator

Inverse of WebMercatorfromLLA — see the docs for that transformation.
"""
struct LLAfromWebMercator <: Transformation
    el::Ellipsoid
end

LLAfromWebMercator(d::Datum=wgs84) = LLAfromWebMercator(ellipsoid(d))


function (trans::WebMercatorfromLLA)(lla::LLA)
    xy = trans(LatLon(lla.lat, lla.lon))
    SA[xy[1], xy[2], lla.alt]
end

function (trans::WebMercatorfromLLA)(ll::LatLon)
    lat_lim = 85.6 # according to https://epsg.io/3857
    if abs(ll.lat) > lat_lim
        throw(ArgumentError("Exceeded Web Mercator latitude bounds"))
    end
    x, y = web_mercator_forward(ll.lat, ll.lon, trans.el.a)
    SA[x, y]
end

function (trans::LLAfromWebMercator)(point::AbstractVector)
    x = point[1]
    y = point[2]
    lat, lon = web_mercator_reverse(x, y, trans.el.a)
    if length(point) == 2
        LatLon(lat, lon)
    elseif length(point) == 3
        LLA(lat, lon, point[3])
    else
        throw(ArgumentError("Expected input vector of length 2 or 3"))
    end
end

Base.inv(trans::WebMercatorfromLLA) = LLAfromWebMercator(trans.el)
Base.inv(trans::LLAfromWebMercator) = WebMercatorfromLLA(trans.el)

# Web / Psuedo Mercator. Following the scaling convention of proj, the scaling
# will be set to the ellipsoid semi-major axis.
#   https://proj.org/operations/projections/webmerc.html
function web_mercator_forward(lat, lon, scaling)
    x = scaling * deg2rad(lon)
    y = scaling * log(tand((90 + lat)/2))
    (x,y)
end

function web_mercator_reverse(x, y, scaling)
    lon = rad2deg(x / scaling)
    lat = 2 * atand(exp(y / scaling)) - 90
    (lat,lon)
end
