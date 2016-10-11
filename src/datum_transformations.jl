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
datum_shift_ECEF{D1,D2}(::Type{D1}, ::Type{D2}) = datum_shift_ECEF(D1(), D2())
datum_shift_ECEF{D1}(::Type{D1}, d2) = datum_shift_ECEF(D1(), d2)
datum_shift_ECEF{D2}(d1, ::Type{D2}) = datum_shift_ECEF(d1, D2())


#-------------------------------------------------------------------------------
# Functions to compute transformations between known datums

#
# ITRF -> ITRF transforms
#
"""
    ITRF_from_ITRF(Ydest, Ysrc, epoch)

Return a `Transformation` converting ECEF points from one ITRF realization to another to GDA94.  Datum
shift parameters are taken from [1], supporting `ITRF_year` 2014, 2008, 2005, 2000, 1997, 1996
1997, 1996, 1994, 1993, 1992, 1991, 1990, 1989, 1988.  `epoch` is the `Date` (or `DateTime`) of interest at which the
input `ECEF` coordinates were measured in ITRF.

[1] http://itrf.ign.fr/trans_para.php

"""

function ITRF_from_ITRF{Ydest, Ysrc}(dest::ITRF{Ydest}, src::ITRF{Ysrc}, epoch)

    # Values comes from http://itrf.ign.fr/trans_para.php
    # (there's a scraping function parse_itrf_src(...) below)
    #
    # (Unit key: ppb = parts per billion; mas = milli arc seconds; yr = year)
    #
    #          Values at 2000 reference epoch                                 Derivatives                                          EPOCH
    #  Name    tx      ty      tz     D       rx       ry      rz             tx      ty    tz     D      rx      ry      rz
    #  Units   mm      mm      mm     ppb     mas      mas     mas            mm/yr  mm/yr  mm/yr  ppb/yr mas/yr  mas/yr  mas/yr
    #
    # dict key is (src, dest)
    table = Dict(
        # 2014 ->
        (2014, 2005) => ([ 2.6000,  1.0000, -2.3000,  0.9200,  0.0000,  0.0000,  0.0000], [ 0.3000,  0.0000, -0.1000,  0.0300,  0.0000,  0.0000,  0.0000], 2010),
        (2014, 2000) => ([ 0.7000,  1.2000, -26.1000,  2.1200,  0.0000,  0.0000,  0.0000], [ 0.1000,  0.1000, -1.9000,  0.1100,  0.0000,  0.0000,  0.0000], 2010),
        (2014, 1997) => ([ 7.4000, -0.5000, -62.8000,  3.8000,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1996) => ([ 7.4000, -0.5000, -62.8000,  3.8000,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1994) => ([ 7.4000, -0.5000, -62.8000,  3.8000,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1993) => ([-50.4000,  3.3000, -60.2000,  4.2900, -2.8100, -3.3800,  0.4000], [-2.8000, -0.1000, -2.5000,  0.1200, -0.1100, -0.1900,  0.0700], 2010),
        (2014, 1992) => ([ 15.4000,  1.5000, -70.8000,  3.0900,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1991) => ([ 27.4000,  15.5000, -76.8000,  4.4900,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1990) => ([ 25.4000,  11.5000, -92.8000,  4.7900,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1989) => ([ 30.4000,  35.5000, -130.8000,  8.1900,  0.0000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        (2014, 1988) => ([ 25.4000, -0.5000, -154.8000,  11.2900,  0.1000,  0.0000,  0.2600], [ 0.1000, -0.5000, -3.3000,  0.1200,  0.0000,  0.0000,  0.0200], 2010),
        # 2008 ->
        (2008, 2005) => ([-2.0000, -0.9000, -4.7000,  0.9400,  0.0000,  0.0000,  0.0000], [ 0.3000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000], 2000),
        (2008, 2000) => ([-1.9000, -1.7000, -10.5000,  1.3400,  0.0000,  0.0000,  0.0000], [ 0.1000,  0.1000, -1.8000,  0.0800,  0.0000,  0.0000,  0.0000], 2000),
        (2008, 1997) => ([ 4.8000,  2.6000, -33.2000,  2.9200,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1996) => ([ 4.8000,  2.6000, -33.2000,  2.9200,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1994) => ([ 4.8000,  2.6000, -33.2000,  2.9200,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1993) => ([-24.0000,  2.4000, -38.6000,  3.4100, -1.7100, -1.4800, -0.3000], [-2.8000, -0.1000, -2.4000,  0.0900, -0.1100, -0.1900,  0.0700], 2000),
        (2008, 1992) => ([ 12.8000,  4.6000, -41.2000,  2.2100,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1991) => ([ 24.8000,  18.6000, -47.2000,  3.6100,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1990) => ([ 22.8000,  14.6000, -63.2000,  3.9100,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1989) => ([ 27.8000,  38.6000, -101.2000,  7.3100,  0.0000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        (2008, 1988) => ([ 22.8000,  2.6000, -125.2000,  10.4100,  0.1000,  0.0000,  0.0600], [ 0.1000, -0.5000, -3.2000,  0.0900,  0.0000,  0.0000,  0.0200], 2000),
        # 2005 ->  N.B. there are more available from alternate sources if needed
        (2005, 2000) => ([0.1 -0.8 -5.8 0.40 0.000 0.000 0.000], [-0.2 0.1 -1.8 0.08 0.000 0.000 0.000], 2000),
        # 2000 ->
        (2000, 1997) => ([ 0.6700,  0.6100, -1.8500,  1.5500,  0.0000,  0.0000,  0.0000], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1997),
        (2000, 1996) => ([ 0.6700,  0.6100, -1.8500,  1.5500,  0.0000,  0.0000,  0.0000], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1997),
        (2000, 1994) => ([ 0.6700,  0.6100, -1.8500,  1.5500,  0.0000,  0.0000,  0.0000], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1997),
        (2000, 1993) => ([ 1.2700,  0.6500, -2.0900,  1.9500, -0.3900,  0.8000, -1.1400], [-0.2900, -0.0200, -0.0600,  0.0100, -0.1100, -0.1900,  0.0700], 1988),
        (2000, 1992) => ([ 1.4700,  1.3500, -1.3900,  0.7500,  0.0000,  0.0000, -0.1800], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1988),
        (2000, 1991) => ([ 2.6700,  2.7500, -1.9900,  2.1500,  0.0000,  0.0000, -0.1800], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1988),
        (2000, 1990) => ([ 2.4700,  2.3500, -3.5900,  2.4500,  0.0000,  0.0000, -0.1800], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1988),
        (2000, 1989) => ([ 2.9700,  4.7500, -7.3900,  5.8500,  0.0000,  0.0000, -0.1800], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1988),
        (2000, 1988) => ([ 2.4700,  1.1500, -9.7900,  8.9500,  0.1000,  0.0000, -0.1800], [ 0.0000, -0.0600, -0.1400,  0.0100,  0.0000,  0.0000,  0.0200], 1988)
    )
    key = (max(Ydest, Ysrc), min(Ydest, Ysrc))
    if !haskey(table, key)
        throw(ErrorException("No $(src) to $(dest) transformation available"))
    end
    params, rates, reference_epoch = table[key]
    reference_epoch = DateTime(reference_epoch,1,1)

    # Fractional years since reference epoch.  Seems a bit unnecessarily cumbersome!
    millisecs_per_year = 365.25 * Dates.value(Dates.Millisecond(Dates.Day(1)))
    dt = Dates.value(Dates.Millisecond(convert(DateTime, epoch) - reference_epoch)) / millisecs_per_year

    # Convert units to meters and radians
    mas2rad = deg2rad(1e-3/(60*60))
    unitconv = [1e-3, 1e-3, 1e-3, 1e-9, mas2rad, mas2rad, mas2rad]
    Tx,Ty,Tz, D, Rx,Ry,Rz = unitconv .* (params .+ rates*dt)

    # assemble the backwards in time version
    M = @SMatrix [1+D  Rz  -Ry;
                 -Rz   1+D  Rx;
                  Ry  -Rx   1+D]
    trans = AffineMap(M, SVector(Tx,Ty,Tz))

    # want the inverse to go forward in time
    if (Ydest > Ysrc)
        trans = inv(trans)
    end
    return trans
end
ITRF_from_ITRF{Y}(::ITRF{Y}, ::ITRF{Y}, epoch) = AffineMap(eye(SMatrix{3}), zero(SVector{3}))

function datum_shift_ECEF{Ydest, Ysrc}(dest::ITRF{Ydest}, src::ITRF{Ysrc})
    @assert dest.epoch == src.epoch "mismatched epochs are not supported when transforming $(src) -> $(dest)"
    ITRF_from_ITRF(dest, src, dest.epoch)
end



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


datum_shift_ECEF{Y}(::GDA94, itrf::ITRF{Y}) = GDA94_from_ITRF_Dawson2010(Y, itrf.epoch)
datum_shift_ECEF{Y}(itrf::ITRF{Y}, ::GDA94) = inv(GDA94_from_ITRF_Dawson2010(Y, itrf.epoch))

# TODO - time-based transformation!
# make_datum_transform_ECEF{Y}(::ITRF{Y,Void}, ::GDA94) =


#
# ETRS related transforms
#

"""
    ETRF_from_ITRF(ITRF_year, epoch)

Return a `Transformation` converting ECEF points from ITRF to ETRF2000.  Datum
shift parameters are taken from [1], supporting `ITRF_year` 2008, 2005, 2000,
1997, 1996, 1994, 1993, 1992, 1991, 1990, 1989.  `epoch` is the `Date` (or `DateTime`) of interest at which the
input `ECEF` coordinates were measured in ITRF.

[1] C. Boucher and Z. Altamimi, "Memo : Specifications for reference frame fixing in the analysis of a
    EUREF GPS campaign" Version 8

"""
function ETRF_from_ITRF{ETRF_year, ITRF_year}(dest::ETRF{ETRF_year}, src::ITRF{ITRF_year}, epoch)
    inter = ITRF{ETRF_year}() # go via this datum
    itrf_trans = ITRF_from_ITRF(inter, src, epoch)
    etrf_trans = ETRF_from_ITRF(dest, inter, epoch)
    compose(etrf_trans, itrf_trans)
end

function ETRF_from_ITRF{Year}(dest::ETRF{Year}, src::ITRF{Year}, epoch)

    # ITRF transformation parameters from
    # C. Boucher and Z. Altamimi, "Memo : Specifications for reference frame fixing in the analysis of a
    # EUREF GPS campaign" Version 8
    #
    # These are intentionally kept exactly the same as the paper for
    # consistency and ease of checking.
    #
    # (Unit key: ppb = parts per billion; mas = milli arc seconds; yr = year)
    #
    #          Values at 1989 reference epoch                                 Derivatives
    #  Name    tx      ty      tz     D       rx       ry      rz             tx      ty    tz     D      rx      ry      rz
    #  Units   cm      cm      cm     1e-8    mas      mas     mas            cm/yr   cm/yr cm/yr  1e-8/yr mas/yr  mas/yr  mas/yr
    #
    table = Dict(
        2005 => ([ 5.6000,  4.8000, -3.7000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.0540,  0.5180, -0.7810]),
        2000 => ([ 5.4000,  5.1000, -4.8000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.0810,  0.4900, -0.7920]),
        1997 => ([ 4.1000,  4.1000, -4.9000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.2000,  0.5000, -0.6500]),
        1996 => ([ 4.1000,  4.1000, -4.9000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.2000,  0.5000, -0.6500]),
        1994 => ([ 4.1000,  4.1000, -4.9000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.2000,  0.5000, -0.6500]),
        1993 => ([ 1.9000,  5.3000, -2.1000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.3200,  0.7800, -0.6700]),
        1992 => ([ 3.8000,  4.0000, -3.7000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.2100,  0.5200, -0.6800]),
        1991 => ([ 2.1000,  2.5000, -3.7000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.2100,  0.5200, -0.6800]),
        1990 => ([ 1.9000,  2.8000, -2.3000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.1100,  0.5700, -0.7100]),
        1989 => ([ 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000], [ 0.0000,  0.0000,  0.0000,  0.0000,  0.1100,  0.5700, -0.7100])
    )
    reference_epoch = DateTime(1989,1,1)

    if !haskey(table, Year)
        throw(ErrorException("No $(src) to $(dest) parameters available"))
    end
    params, rates = table[Year]

    # Fractional years since reference epoch.  Seems a bit unnecessarily cumbersome!
    millisecs_per_year = 365.25 * Dates.value(Dates.Millisecond(Dates.Day(1)))
    dt = Dates.value(Dates.Millisecond(convert(DateTime, epoch) - reference_epoch)) / millisecs_per_year

    # Convert units to meters and radians
    mas2rad = deg2rad(1e-3/(60*60))
    unitconv = [1e-2, 1e-2, 1e-2, 1e-8, mas2rad, mas2rad, mas2rad]
    Tx,Ty,Tz, D, Rx,Ry,Rz = unitconv .* (params .+ rates*dt)
    M =  @SMatrix [1+D  Rz  -Ry;
                  -Rz   1+D  Rx;
                   Ry  -Rx   1+D]
    AffineMap(M, SVector(Tx,Ty,Tz))
end


# special case for ETRF2000, we know the transformation to it
# TODO: is this worth keeping?
function ETRF_from_ITRF{ITRF_year}(::ETRF{2000}, ::ITRF{ITRF_year}, epoch)
    # ITRF transformation parameters from
    # C. Boucher and Z. Altamimi, "Memo : Specifications for reference frame fixing in the analysis of a
    # EUREF GPS campaign" Version 8
    #
    # These are intentionally kept exactly the same as the paper for
    # consistency and ease of checking.
    #
    # (Unit key: ppb = parts per billion; mas = milli arc seconds; yr = year)
    #
    #          Values at 2000 reference epoch                                 Derivatives
    #  Name    tx      ty      tz     D       rx       ry      rz             tx      ty    tz     D      rx      ry      rz
    #  Units   mm      mm      mm     ppb     mas      mas     mas            mm/yr  mm/yr  mm/yr  ppb/yr mas/yr  mas/yr  mas/yr
    #
    table = Dict(
        2008 => ([52.1, 49.3, -58.5,  1.34, 0.891, 5.390, -8.712], [0.1,  0.1, -1.8,  0.08, 0.081, 0.490, -0.792]),
        2005 => ([54.1, 50.2, -53.8,  0.40, 0.891, 5.390, -8.712], [-0.2, 0.1, -1.8,  0.08, 0.081, 0.490, -0.792]),
        2000 => ([54.0, 51.0, -48.0,  0.00, 0.891, 5.390, -8.712], [0.0,  0.0,  0.0,  0.00, 0.081, 0.490, -0.792]),
        1997 => ([47.3, 46.7, -25.3, -1.58, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812]),
        1996 => ([47.3, 46.7, -25.3, -1.58, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812]),
        1994 => ([47.3, 46.7, -25.3, -1.58, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812]),
        1993 => ([76.1, 46.9, -19.9, -2.07, 2.601, 6.870, -8.412], [2.9,  0.2,  0.6, -0.01, 0.191, 0.680, -0.862]),
        1992 => ([39.3, 44.7, -17.3, -0.87, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812]),
        1991 => ([27.3, 30.7, -11.3, -2.27, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812]),
        1990 => ([29.3, 34.7,  4.7,  -2.57, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812]),
        1989 => ([24.3, 10.7,  42.7, -5.97, 0.891, 5.390, -8.772], [0.0,  0.6,  1.4, -0.01, 0.081, 0.490, -0.812])
    )
    reference_epoch = DateTime(2000,1,1)

    if !haskey(table, ITRF_year)
        throw(ErrorException("No ITRF to ETRF2000 parameters for ITRF$(ITRF_year)"))
    end
    params, rates = table[ITRF_year]

    # Fractional years since reference epoch.  Seems a bit unnecessarily cumbersome!
    millisecs_per_year = 365.25 * Dates.value(Dates.Millisecond(Dates.Day(1)))
    dt = Dates.value(Dates.Millisecond(convert(DateTime, epoch) - reference_epoch)) / millisecs_per_year

    # Convert units to meters and radians
    mas2rad = deg2rad(1e-3/(60*60))
    unitconv = [1e-3, 1e-3, 1e-3, 1e-9, mas2rad, mas2rad, mas2rad]
    Tx,Ty,Tz, D, Rx,Ry,Rz = unitconv .* (params .+ rates*dt)
    M = @SMatrix [1+D  Rz  -Ry;
                 -Rz   1+D  Rx;
                  Ry  -Rx   1+D]

    AffineMap(M, SVector(Tx,Ty,Tz))
end

datum_shift_ECEF{Ydest, Ysrc}(dest::ETRF{Ydest}, src::ITRF{Ysrc}) = ETRF_from_ITRF(dest, src, src.epoch)
datum_shift_ECEF{Ydest, Ysrc}(dest::ITRF{Ydest}, src::ETRF{Ysrc}) = inv(ETRF_from_ITRF(src, dest, dest.epoch))




#
# Helper functions
#

# parse the itrf conversion files and return a dict
# example http://itrf.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt
# hope they keep using the same format
function parse_itrf_src(file)

    # (year_out, year_in) => (params, rates, epeoch)
    dict = Dict{NTuple{2,Int}, Tuple{Vector{Float64}, Vector{Float64}, Int}}()

    # read it all
    local lines
    open(file) do fid

        # grab the realization year
        header = readline(fid)
        src_year = match(r"(?<=ITRF)\d+", header).match
        src_year = (length(src_year) < 4 ? 1900 : 0) + parse(Int, src_year)

        # now read till we get into the parameters part
        line = chomp(readline(fid))
        done = false
        while !ismatch(r"(?<=ITRF)\d+", line)
            line = chomp(readline(fid))
        end

        # and loop till thge end of the table
        while ismatch(r"(?<=ITRF)\d+", line)

            # read the year
            dest_year = match(r"(?<=ITRF)\d+", line).match
            dest_year = (length(dest_year) < 4 ? 1900 : 0) + parse(Int, dest_year)

            # read params
            params = Vector{Float64}(0)
            for str in eachmatch(r"(?<=\s)[-+]*\d+\.\d+", line); push!(params, parse(Float64, str.match)); end
            epoch = Int(pop!(params))
            @assert length(params) == 7

            # read rates
            line = chomp(readline(fid))
            rates = Vector{Float64}(0)
            for str in eachmatch(r"(?<=\s)[-+]*\d+\.\d+", line); push!(rates, parse(Float64, str.match)); end
            @assert length(rates) == 7

            # add it to the dict
            dict[(src_year, dest_year)] = (params, rates, epoch)
            @printf("(%i, %i) => ([% 0.4f, % 0.4f, % 0.4f, % 0.4f, % 0.4f, % 0.4f, % 0.4f], [% 0.4f, % 0.4f, % 0.4f, % 0.4f, % 0.4f, % 0.4f, % 0.4f], %i),\n",
                     src_year, dest_year, params..., rates..., epoch)

            # move along
            line = chomp(readline(fid))
        end
    end
    return dict
end
