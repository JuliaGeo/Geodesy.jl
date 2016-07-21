# Generic datum shifting function in ECEF coordinates

"""
    datum_shift_ECEF(dest_datum, source_datum)

Return a transformation object which transforms ECEF points from `source_datum`
to `dest_datum`, provided there's a well-accepted transformation available.

Note that any computed datum shift is inherently prone to measurement error,
with different physical measurements resulting in different transformations.
This is in contrast to coordinate system transformations which are precisely
mathematically defined.  Because of this, `datum_shift_ECEF`

Suppose you've got two points with different datums.  These were measured with
respect to different physical reference frames.

Examples of measurement systems:
* Datum 1: Measure the angles to the top of the church steeple and the big
  fence post on the hill.  Knowing the distance between those accuratly, your
  *relative* position may be calculated using triangulation.
* Datum 2: Measure the time and phase of GPS satellite signals .. FIXME


This means that the transformation parameters between datums
are measured quantities relating two different physical measurement procedures.
As such, different organizations will not agree on the exact way to transform
between datums if they have a different set of observations with which to fit
the transformation.

FIXME DOCS!
"""
datum_shift_ECEF(dest_datum, source_datum) = error("No datum shift implemented for $dest_datum from $source_datum")

# Convenience function for datums types without parameters
datum_shift_ECEF{D1,D2}(::Type{D1}, ::Type{D2}) = datum_shift_ECEF(D1(), D2())
datum_shift_ECEF{D1}(::Type{D1}, d2) = datum_shift_ECEF(D1(), d2)
datum_shift_ECEF{D2}(d1, ::Type{D2}) = datum_shift_ECEF(d1, D2())


#-------------------------------------------------------------------------------
# Functions to compute transformations between known datums

"""
    GDA94_from_ITRF_Dawson2010(ITRF_year, epoch)

Return a `Transformation` converting ECEF points from a given ITRF to GDA94.
Datum shift parameters are taken from [1], supporting `ITRF_year` 2008, 2005,
2000, 1997, 1996.  `epoch` is the `Date` (or `DateTime`) of interest at which
the input `ECEF` coordinates to the transformation were measured in ITRF.

[1] J. Dawson and A. Woods, "ITRF to GDA94 coordinate transforms",
    Journal of Applied Geodesy, 4, p. 189 (2010).

TODO: We don't yet support `epoch` varying per input point, but there should be
a `Transformation` object for this at some stage.
"""
function GDA94_from_ITRF_Dawson2010(ITRF_year, epoch)
    # ITRF transformation parameters from
    # J. Dawson and A. Woods, "ITRF to GDA94 coordinate transforms",
    # Journal of Applied Geodesy, 4, p. 189 (2010):
    #
    # These are intentionally kept exactly the same as the paper for
    # consistency and ease of checking.
    #
    # (Unit key: ppb = parts per billion; mas = milli arc seconds; yr = year)
    #
    #          Values at 1994 reference epoch                                 Derivatives
    #  Name    tx      ty      tz     sc      rx       ry      rz             tx      ty    tz     sc     rx      ry      rz
    #  Units   mm      mm      mm     ppb     mas      mas     mas            mm/yr  mm/yr  mm/yr  ppb/yr mas/yr  mas/yr  mas/yr
    #
    table = Dict(
        2008 => ([-84.68, -19.42,  32.01, 9.710, -0.4254,  2.2578, 2.4015], [  1.42,  1.34,  0.90, 0.109, 1.5461, 1.1820, 1.1551]),
        2005 => ([-79.73,  -6.86,  38.03, 6.636, -0.0351,  2.1211, 2.1411], [  2.25, -0.62, -0.56, 0.294, 1.4707, 1.1443, 1.1701]),
        2000 => ([-45.91, -29.85, -20.37, 7.070, -1.6705,  0.4594, 1.9356], [ -4.66,  3.55, 11.24, 0.249, 1.7454, 1.4868, 1.2240]),
        1997 => ([-14.63, -27.62, -25.32, 6.695, -1.7893, -0.6047, 0.9962], [ -8.60,  0.36, 11.25, 0.007, 1.6394, 1.5198, 1.3801]),
        1996 => ([ 24.54, -36.43, -68.12, 6.901, -2.7359, -2.0431, 0.3731], [-21.80,  4.71, 26.27, 0.388, 2.0203, 2.1735, 1.6290]),
    )
    reference_epoch = DateTime(1994,1,1)

    if !haskey(table, ITRF_year)
        throw(ErrorException("No ITRF to GDA94 parameters for ITRF year $ITRF_year"))
    end
    params, rates = table[ITRF_year]

    # Fractional years since reference epoch.  Seems a bit unnecessarily cumbersome!
    millisecs_per_year = 365.25 * Dates.value(Dates.Millisecond(Dates.Day(1)))
    dt = Dates.value(Dates.Millisecond(convert(DateTime, epoch) - reference_epoch)) / millisecs_per_year

    # Convert units to meters and radians
    mas2rad = deg2rad(1e-3/(60*60))
    unitconv = [1e-3, 1e-3, 1e-3, 1e-9, mas2rad, mas2rad, mas2rad]
    Tx,Ty,Tz, Sc, Rx,Ry,Rz = unitconv .* (params .+ rates*dt)

    M = (1+Sc) * @fsa [1.0  Rz  -Ry;
                      -Rz  1.0   Rx;
                       Ry  -Rx  1.0]

    AffineTransformation(M, Vec(Tx,Ty,Tz))
end


datum_shift_ECEF{Y}(::GDA94, itrf::ITRF{Y}) = GDA94_from_ITRF_Dawson2010(Y, itrf.epoch)
datum_shift_ECEF{Y}(itrf::ITRF{Y}, ::GDA94) = inv(GDA94_from_ITRF_Dawson2010(Y, itrf.epoch))

# TODO - time-based transformation!
# make_datum_transform_ECEF{Y}(::ITRF{Y,Void}, ::GDA94) =

