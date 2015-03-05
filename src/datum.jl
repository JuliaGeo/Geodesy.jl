### Ellipsoid
# Specify datum for translation between LLA and other coordinate systems
immutable Ellipsoid
    a::Float64        # Semi-major axis
    b::Float64        # Semi-minor axis
    e²::Float64       # Eccentricity squared
    e′²::Float64      # Second eccentricity squared
end

function Ellipsoid(; a::String="", b::String="", f_inv::String="")
    if isempty(a) || isempty(b) == isempty(f_inv)
        throw(ArgumentError("Specify parameter 'a' and either 'b' or 'f_inv'"))
    end
    if isempty(b)
        _ellipsoid_af(BigFloat(a), BigFloat(f_inv))
    else
        _ellipsoid_ab(BigFloat(a), BigFloat(b))
    end
end

function _ellipsoid_ab(a::BigFloat, b::BigFloat)
    e² = (a^2 - b^2) / a^2
    e′² = (a^2 - b^2) / b^2

    Ellipsoid(a, b, e², e′²)
end
function _ellipsoid_af(a::BigFloat, f_inv::BigFloat)
    b = a * (1 - inv(f_inv))

    _ellipsoid_ab(a, b)
end

### World Geodetic Coordinate System of 1984 (WGS 84)
# Standardized coordinate system for Earth
# Global ellipsoidal reference surface
const WGS84  = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")

const OSGB36 = Ellipsoid(a = "6377563.396", b = "6356256.909")
const NAD27  = Ellipsoid(a = "6378206.4",   b = "6356583.8")
