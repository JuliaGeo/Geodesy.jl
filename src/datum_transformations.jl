# Functions to compute transformations between known datums

"""
    GDA94_from_ITRF(ITRF_realization, epoch)

Compute ITRF to GDA94 datum shift using parameters from [1].  The year of the
desired ITRF realization should be given by `ITRF_realization`.  `epoch` is the
Date (or DateTime) of interest at which the input `ECEF` coordinates were
measured in ITRF.

[1] J. Dawson and A. Woods, "ITRF to GDA94 coordinate transforms",
    Journal of Applied Geodesy, 4, p. 189 (2010)
"""
function GDA94_from_ITRF(ITRF_realization, epoch)
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

    if !haskey(table, ITRF_realization)
        throw(ErrorException("No ITRF to GDA94 parameters for ITRF realization $ITRF_realization"))
    end
    params, rates = table[ITRF_realization]

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

