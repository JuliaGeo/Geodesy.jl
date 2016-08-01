"""
Abstract type for geodetic datums

A datum is a set of reference objects and assigned coordinates, relative to
which other objects may be positioned.  We model these in code with subtypes of
`Datum`.  Each geodetic datum has an associated ellipsoid model of the Earth
which is required when transforming between coordinate systems.  The ellipsoid
can be accessed with the `ellipsoid()` function.
"""
abstract Datum

#------------------------------------------------------------------------------
# Worldwide geodetic datums

const valid_WGS84_frame_weeks = [0, 730,  873, 1150, 1674, 1762]
const WGS84_frames_last_updated = 1907

"""
    WGS84()

Construct an object representing the World Geodetic System of 1984, as a dynamic
datum - see below for gory details about specific WGS84 frames).

If you're getting positions from a consumer GPS device, you're probably going to
have WGS84 by default because it's the datum in which GPS satellites
broadcast their position ("broadcast ephemerides").  Note however that many
devices can also provide position in a national datum, so you should check your
device settings to be sure.

As a special case for low accuracy work (worse than a meter or so), `WGS84` will
assume that coordinates supplied without a capture time are in the *latest*
frame realization known to Geodesy.jl, WGS84 (G$(valid_WGS84_frame_weeks[end])).
Note that this may not be correct if you're processing historical data, or
Geodesy.jl itself is out of date.  For higher accuracy, you should supply a date
of capture - either with each coordinate, or explicitly using the `GpsWeek`
parameter to the WGS84 type:

    WGS84{GpsWeek}()

Construct an object representing the static WGS84 datum computed using data
gathered prior to the given `GpsWeek`.  WGS84 is maintained and updated by the
US National Geospatial-Intelligence-Agency (NGA) at irregular intervals to align
with the ITRF to within 0.1m (see [1]); if you care about accuracy at that
level, you probably want to be solving for position in a different datum, for
example, ITRF.  As of 2016, `GpsWeek` should be one out of
$valid_WGS84_frame_weeks.

Note that the dates of implementation of these frames as broadcast by the
satellites are not the same as the associated GPS week - see Ref. [1], table
2.1.  (TODO: Perhaps Geodesy should have a table to figure out which frame to
use at a given date?  Does anybody care?)

1. "World Geodetic System 1984", NGA standard NGA.STND.0036_1.0.0_WGS84, 2014-07-08,
   http://earth-info.nga.mil/GandG/publications/NGA_STND_0036_1_0_0_WGS84/NGA.STND.0036_1.0.0_WGS84.pdf
"""
immutable WGS84{GpsWeek} <: Datum
    WGS84() = check_wgs84_params(new())
end

WGS84() = WGS84{Void}()

check_wgs84_params(wgs84::WGS84{Void}) = wgs84

@generated function check_wgs84_params{GpsWeek}(wgs84::WGS84{GpsWeek})
    if GpsWeek <= WGS84_frames_last_updated && GpsWeek ∉ valid_WGS84_frame_weeks
        :(throw(ErrorException("No WGS84 realization exists for week $GpsWeek")))
    else
        :(wgs84)
    end
end

Base.show(io::IO, ::WGS84{Void}) = print(io,"WGS84")
Base.show{W}(io::IO, ::WGS84{W}) = print(io,"WGS84 (G$W)")


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
frame processing regression problem.  A list of recent realizations is available
at http://itrf.ensg.ign.fr/ITRF_solutions; as of 2016-07 valid `Year`s were
2014 2008 2005 2000 1997 1996 1994 1993 1992 1991 1990 1989 1988.

See http://itrf.ensg.ign.fr/general.php for a technical overview.  Useful
technical papers:

* "IERS Conventions (2010)", Petit and Luzum (eds.), IERS Technical note No.36,
   Chapter 4, https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
* "ITRF2008: an improved solution of the international terrestrial reference
   frame", Altamimi et al., J. Geodesy (2011) 85: 457,
   http://dx.doi.org/10.1007/s00190-011-0444-4
* "IGS08: the IGS realization of ITRF2008", Rebischung et al., GPS Solutions
   (2012) 16: 483, http://dx.doi.org/10.1007/s10291-011-0248-2,
   ftp://igs.org/pub/resource/pubs/IGS08_The_IGS_Realization_of_ITRF2008.pdf
* "The IGS contribution to ITRF2014", Rebischung et al., J. Geodesy (2016) 90:
   611, http://dx.doi.org/10.1007/s00190-016-0897-6
"""
immutable ITRF{Year, EpochT} <: Datum
    epoch::EpochT

    function ITRF(epoch::EpochT)
        check_itrf_year(new(epoch))
    end
end

@generated function check_itrf_year{Y}(itrf::ITRF{Y})
    if Y < 100
        :(throw(ErrorException("Two digit year $Y for ITRF found - this library expects a full four-digit year.")))
    elseif Y <= 2016 && Y ∉ [2014, 2008, 2005, 2000, 1997, 1996,
                             1994, 1993, 1992, 1991, 1990, 1989, 1988]
        :(throw(ErrorException("No ITRF realization exists for year $Y")))
    else
        :(itrf)
    end
end

@compat (::Type{ITRF{Year}}){Year}() = ITRF{Year,Void}(nothing)
@compat (::Type{ITRF{Year}}){Year}(epoch) = ITRF{Year,typeof(epoch)}(epoch)

Base.show{Y}(io::IO, ::ITRF{Y,Void}) = print(io,"ITRF{$Y}")
Base.show(io::IO, itrf::ITRF) = print(io,"ITRF{$Y}($(itrf.epoch))")

#------------------------------------------------------------------------------
# National geodetic datums
"""
`OSGB36` - Datum for Ordinance Survey of Great Britian, 1936
"""
immutable OSGB36 <: Datum; end
Base.show(io::IO, ::OSGB36) = print(io,"osbg84")

"""
`NAD27` - North American Datum of 1927
"""
immutable NAD27 <: Datum; end
Base.show(io::IO, ::NAD27) = print(io,"nad27")

"""
`GDA94` - Geocentric Datum of Australia, 1994
"""
immutable GDA94 <: Datum; end


#-------------------------------------------------------------------------------
# FIXME: Remove - this is not a datum!!
immutable GRS80 <: Datum; end
Base.show(io::IO, ::GRS80) = print(io,"grs80")

const wgs84 = WGS84()
const osgb36 = OSGB36()
const grs80 = GRS80()
const nad27 = NAD27()

Base.deprecate(:grs80)

# TODO: Figure out how to deprecate these - will we be using types or instances
# for passing datums around.
#Base.deprecate(:wgs84)
#Base.deprecate(:osgb36)
#Base.deprecate(:nad27)

#-------------------------------------------------------------------------------
# Ellipsoids

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

const wgs84_el  = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const airy1830_el = Ellipsoid(a = "6377563.396", b = "6356256.909")
const clarke1866_el  = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const grs80_el  = Ellipsoid(a = "6378137.0", f_inv = "298.2572221008827112431628366") #NB: not the definition - f_inv is derived.

# TODO: Also deprecate other incorrectly named ellipsoids
Base.@deprecate_binding nad27_el clarke1866_el
Base.@deprecate_binding osgb36_el airy1830_el

@inline ellipsoid(::Type{WGS84}) = wgs84_el
@inline ellipsoid(::WGS84)    = wgs84_el
@inline ellipsoid{GpsWeek}(::Type{WGS84{GpsWeek}}) = wgs84_el

@inline ellipsoid(::Union{OSGB36,Type{OSGB36}}) = airy1830_el
@inline ellipsoid(::Union{NAD27,Type{NAD27}})   = clarke1866_el
@inline ellipsoid(::Union{GRS80,Type{GRS80}})   = grs80_el
@inline ellipsoid(::Union{GDA94,Type{GDA94}})   = grs80_el

@inline ellipsoid(::ITRF)                       = grs80_el
@inline ellipsoid{D<:ITRF}(::Type{D})           = grs80_el

ellipsoid(x) = error("No ellipsoid defined for datum $x")
