# Built-in geodetic datums (singletons types and instances)

# FIXME: WGS84 should probably have the realization in the type parameters.
# Though this may confuse people who just want to use "WGS84" naively...
immutable WGS84; end
const wgs84 = WGS84()
Base.show(io::IO, ::WGS84) = print(io,"wgs84")

"""
`OSGB36` - Datum for Ordinance Survey of Great Britian, 1936
"""
immutable OSGB36; end
const osgb36 = OSGB36()
Base.show(io::IO, ::OSGB36) = print(io,"osbg84")

"""
`NAD27` - North American Datum of 1927
"""
immutable NAD27; end
const nad27 = NAD27()
Base.show(io::IO, ::NAD27) = print(io,"nad27")

# FIXME: Remove - this is not a datum!
immutable GRS80; end
const grs80 = GRS80()
Base.show(io::IO, ::GRS80) = print(io,"grs80")
Base.deprecate(:grs80)

"""
`GDA94` - Geocentric Datum of Australia, 1994
"""
immutable GDA94; end


"""
    ITRF{Year}([epoch])

Construct an object representing the International Terrestrial Reference Frame
for the given `Year` of realization.  ITRF is the standard high accuarcy
terrestrial reference frame for worldwide use.  An optional `epoch` parameter
defines the time of interest - typically a date at which coordinates were
measured, using, eg a GPS device.  Without the `epoch` parameter, the resulting
`ITRF{Year}` object represents the full dynamic datum.

A *realization* is created every few years by computing the position of a large
set of ground control stations from satellite and celestial measurements.  The
`Year` parameter represents the last year from which data was used in the
frame processing regression problem.  A full list of realizations is available
at http://itrf.ensg.ign.fr/ITRF_solutions; as of 2016-07 this included
ITRF2014 ITRF2008 ITRF2005 ITRF2000 ITRF1997 ITRF1996 ITRF1994 ITRF1993
ITRF1992.

See http://itrf.ensg.ign.fr/general.php for a technical overview.  Useful
technical papers:

1. "ITRF2008: an improved solution of the international terrestrial reference
    frame", Altamimi et al., J. Geodesy (2011) 85: 457,
    http://dx.doi.org/10.1007/s00190-011-0444-4
3. "IGS08: the IGS realization of ITRF2008", Rebischung et al., GPS Solutions
    (2012) 16: 483, http://dx.doi.org/10.1007/s10291-011-0248-2,
    ftp://igs.org/pub/resource/pubs/IGS08_The_IGS_Realization_of_ITRF2008.pdf
2. "The IGS contribution to ITRF2014", Rebischung et al., J. Geodesy (2016) 90:
    611, http://dx.doi.org/10.1007/s00190-016-0897-6
"""
immutable ITRF{Year, EpochT}
    epoch::EpochT

    function ITRF(epoch::EpochT)
        check_itrf_year(new(epoch))
    end
end

@generated function check_itrf_year{Y}(itrf::ITRF{Y})
    if Y < 100
        :(throw(ErrorException("Two digit year $Y for ITRF found - this library expects a full four-digit year.")))
    elseif Y <= 2016 && Y ∉ [2014, 2008, 2005, 2000, 1997, 1996, 1994, 1993, 1992]
        :(throw(ErrorException("No ITRF realization exists for year $Y")))
    else
        :(itrf)
    end
end

@compat (::Type{ITRF{Year}}){Year}() = ITRF{Year,Void}(nothing)
@compat (::Type{ITRF{Year}}){Year}(epoch) = ITRF{Year,typeof(epoch)}(epoch)

Base.show{Y}(io::IO, ::ITRF{Y,Void}) = print(io,"ITRF{$Y}")
Base.show(io::IO, itrf::ITRF) = print(io,"ITRF{$Y}($(itrf.epoch))")


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
const wgs84_el  = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const osgb36_el = Ellipsoid(a = "6377563.396", b = "6356256.909")
const nad27_el  = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const grs80_el  = Ellipsoid(a = "6378137.0", f_inv = "298.2572221008827112431628366")

@inline ellipsoid(::Union{WGS84,Type{WGS84}})   = wgs84_el
@inline ellipsoid(::Union{OSGB36,Type{OSGB36}}) = osgb36_el
@inline ellipsoid(::Union{NAD27,Type{NAD27}})   = nad27_el
@inline ellipsoid(::Union{GRS80,Type{GRS80}})   = grs80_el
@inline ellipsoid(::Union{GDA94,Type{GDA94}})   = grs80_el
@inline ellipsoid(::ITRF)                       = grs80_el
@inline ellipsoid{D<:ITRF}(::Type{D})           = grs80_el
@inline ellispoid(x::Ellipsoid) = x
ellipsoid(x) = error("No ellipsoid defined for datum $x")
