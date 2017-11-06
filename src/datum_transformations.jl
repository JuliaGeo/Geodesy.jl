# Generic datum shifting function in ECEF coordinates

"""
    datum_shift_ECEF(dest_datum, source_datum)

Return a transformation object which transforms ECEF points from `source_datum`
to `dest_datum`.  The function should attempt to supply the current publicly
accepted best estimate.

Note that the best known version of a datum transformation will inherently
improve with time (see below), so we *cannot simultaneously* guarantee that:

1. We return the best publicly accepted version of a datum shift.
2. We return the same thing in future versions of Geodesy.jl.

This function will attempt to satisfy condition 1. rather than 2.; if you want a
stable version of a transformation you should use one of the lower level
functions, for example, `GDA94_from_ITRF_Dawson2010()`.


### Important note about accuracy

If you care about accuracy, you should *always* store long term archival data in
the source datum where possible, along with metadata defining the full datum and
coordinate system in use.  The time to do a datum shift is when you want to
compare information from two different datums.

Why all this bother about inaccuracy? The errors are twofold: First, the
parameters of a datum shift come from a physical measurement process. Typically
this involves measuring the coordinates of some physical locations in *both*
datums, and a physical measurement procedure is always subject to some inaccuracy.
Second, these **tie points** are used to infer a compact representation of the
datum shift, with as few numerical parameters as possible. A small number of
parameters will result in an overly smooth representation; this modelling error
is a second source of inaccuracy.  Both these errors can be reduced if you
choose to measure more tie points or improve the complexity of the numerical
model at a future date.
"""
datum_shift_ECEF(dest_datum, source_datum) = error("No datum shift implemented for $dest_datum from $source_datum")

# Convenience function for datums types without parameters
datum_shift_ECEF(::Type{D1}, ::Type{D2}) where {D1,D2} = datum_shift_ECEF(D1(), D2())
datum_shift_ECEF(::Type{D1}, d2) where {D1} = datum_shift_ECEF(D1(), d2)
datum_shift_ECEF(d1, ::Type{D2}) where {D2} = datum_shift_ECEF(d1, D2())


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

    M = (1+Sc) * @SMatrix [1.0  Rz  -Ry;
                          -Rz  1.0   Rx;
                           Ry  -Rx  1.0]

    AffineMap(M, SVector(Tx,Ty,Tz))
end


datum_shift_ECEF(::GDA94, itrf::ITRF{Y}) where {Y} = GDA94_from_ITRF_Dawson2010(Y, itrf.epoch)
datum_shift_ECEF(itrf::ITRF{Y}, ::GDA94) where {Y} = inv(GDA94_from_ITRF_Dawson2010(Y, itrf.epoch))

# TODO - time-based transformation!
# make_datum_transform_ECEF{Y}(::ITRF{Y,Void}, ::GDA94) =
