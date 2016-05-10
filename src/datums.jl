# Built-in geodetic datums (singletons types and instances)
immutable WGS84; end
const wgs84 = WGS84()
Base.show(io::IO, ::WGS84) = print(io,"wgs84")

immutable OSGB36; end
const osgb36 = OSGB36()
Base.show(io::IO, ::OSGB36) = print(io,"osbg84")

immutable NAD27; end
const nad27 = NAD27()
Base.show(io::IO, ::NAD27) = print(io,"nad27")

immutable GRS80; end
const grs80 = GRS80()
Base.show(io::IO, ::GRS80) = print(io,"grs80")

"""
An ellipsoidal representation of the earth, for converting between LLA and
other co-ordinate systems such as ECEF.
"""
immutable Ellipsoid
    a::Float64        # Semi-major axis
    b::Float64        # Semi-minor axis
    e²::Float64       # Eccentricity squared
    e′²::Float64      # Second eccentricity squared
end

function Ellipsoid(; a::@compat(AbstractString)="", b::@compat(AbstractString)="", f_inv::@compat(AbstractString)="")
    if isempty(a) || isempty(b) == isempty(f_inv)
        throw(ArgumentError("Specify parameter 'a' and either 'b' or 'f_inv'"))
    end
    if isempty(b)
        _ellipsoid_af(@compat(parse(BigFloat,a)), @compat(parse(BigFloat,f_inv)))
    else
        _ellipsoid_ab(@compat(parse(BigFloat,a)), @compat(parse(BigFloat,b)))
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
const WGS84_el  = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const OSGB36_el = Ellipsoid(a = "6377563.396", b = "6356256.909")
const NAD27_el  = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const GRS80_el  = Ellipsoid(a = "6378137.0", f_inv = "298.2572221008827112431628366")

@inline ellipsoid(::Union{WGS84,Type{WGS84}})   = WGS84_el
@inline ellipsoid(::Union{OSGB36,Type{OSGB36}}) = OSGB36_el
@inline ellipsoid(::Union{NAD27,Type{NAD27}})   = NAD27_el
@inline ellipsoid(::Union{GRS80,Type{GRS80}})   = GRS80_el
@inline ellispoid(x::Ellipsoid) = x
ellipsoid(x) = error("No ellipsoid defined for datum $x")
