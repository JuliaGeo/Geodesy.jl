"""Define a GeodesicLine object in terms of the direct geodesic
problem specified in terms of spherical arc length

:param lat1: latitude of the first point in degrees
:param lon1: longitude of the first point in degrees
:param azi1: azimuth at the first point in degrees
:param s12: the distance from the first point to the second in
meters
:param caps: the :ref:`capabilities <outmask>`
:return: a :class:`~geographiclib.geodesicline.GeodesicLine`

This function sets point 3 of the GeodesicLine to correspond to
point 2 of the direct geodesic problem.  The default value of *caps*
is STANDARD | DISTANCE_IN, allowing direct geodesic problem to be
solved.
"""
function DirectLine(geod::Geodesic, lat1, lon1, azi1, s12,
               caps = GeodesicCapability.STANDARD | GeodesicCapability.DISTANCE_IN)
    _GenDirectLine(geod, lat1, lon1, azi1, false, s12, caps)
end

"""Define a GeodesicLine object in terms of the direct geodesic
problem specified in terms of spherical arc length

:param lat1: latitude of the first point in degrees
:param lon1: longitude of the first point in degrees
:param azi1: azimuth at the first point in degrees
:param a12: spherical arc length from the first point to the second
in degrees
:param caps: the :ref:`capabilities <outmask>`
:return: a :class:`~geographiclib.geodesicline.GeodesicLine`

This function sets point 3 of the GeodesicLine to correspond to
point 2 of the direct geodesic problem.  The default value of *caps*
is STANDARD | DISTANCE_IN, allowing direct geodesic problem to be
solved.
"""
function ArcDirectLine(geod::Geodesic, lat1, lon1, azi1, a12,
               caps = GeodesicCapability.STANDARD | GeodesicCapability.DISTANCE_IN)
    _GenDirectLine(geod, lat1, lon1, azi1, true, a12, caps)
end

"""Define a GeodesicLine object in terms of the invese geodesic problem

:param lat1: latitude of the first point in degrees
:param lon1: longitude of the first point in degrees
:param lat2: latitude of the second point in degrees
:param lon2: longitude of the second point in degrees
:param caps: the :ref:`capabilities <outmask>`
:return: a :class:`~geographiclib.geodesicline.GeodesicLine`

This function sets point 3 of the GeodesicLine to correspond to
point 2 of the inverse geodesic problem.  The default value of *caps*
is STANDARD | DISTANCE_IN, allowing direct geodesic problem to be
solved.
"""
function InverseLine(geod::Geodesic, lat1, lon1, lat2, lon2,
        caps=GeodesicCapability.STANDARD | GeodesicCapability.DISTANCE_IN)
    a12, _, salp1, calp1, _, _, _, _, _, _ = Geodesics._GenInverse(geod, lat1, lon1, lat2, lon2, 0)
    azi1 = Math.atan2d(salp1, calp1)
    if (caps & (Geodesics.OUT_MASK & Geodesics.DISTANCE_IN)) > 0
        caps |= Geodesics.DISTANCE
    end
    line = GeodesicLine(geod, lat1, lon1, azi1, caps=caps, salp1=salp1, calp1=calp1)
    line = SetArc(line, a12)
    line
end

"""Private: general form of DirectLine"""
function _GenDirectLine(geod::Geodesic, lat1, lon1, azi1, arcmode::Bool, s12_a12,
                   caps=GeodesicCapability.STANDARD | GeodesicCapability.DISTANCE_IN)
    # Automatically supply DISTANCE_IN if necessary
    if !arcmode
        caps |= Geodesics.DISTANCE_IN
    end
    line = GeodesicLine(geod, lat1, lon1, azi1, caps=caps)
    line = if arcmode
        SetArc(line, s12_a12)
    else
        SetDistance(line, s12_a12)
    end
    line
end

