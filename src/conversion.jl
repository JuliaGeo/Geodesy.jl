
############################
### Identity conversions ###
############################

ECEF(ecef::ECEF, datum) = ecef
LLA(lla::LLA, datum) = lla
ENU(enu::ENU, datum) = enu


################################
### LLA <-> ECEF coordinates ###
################################

ECEF(lla::LLA, datum) = transform(ECEFfromLLA(datum), lla)

LLA(ecef::ECEF, datum) = transform(LLAfromECEF(datum), ecef)

################################
### ECEF <-> ENU coordinates ###
################################

ENU(ecef::ECEF, origin_lla::LLA, datum) = transform(ENUfromECEF(origin_lla, datum), ecef)
ENU(ecef::ECEF, origin_ecef::ECEF, datum) = transform(ENUfromECEF(origin_ecef, datum), ecef)

ENU(lla::LLA, origin_lla::LLA, datum) =  transform(ENUfromLLA(origin_lla, datum), lla)
ENU(lla::LLA, origin_ecef::ECEF, datum) =  transform(ENUfromLLA(origin_ecef, datum), lla)

#=
# Given a reference point for linarization
function ENU(ecef::ECEF, lla_ref::LLA, datum::Ellipsoid)
    ϕdeg, λdeg = lla_ref.lat, lla_ref.lon

    ecef_ref = ECEF(lla_ref, datum)
    ∂x = ecef.x - ecef_ref.x
    ∂y = ecef.y - ecef_ref.y
    ∂z = ecef.z - ecef_ref.z

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # R = [     -sinλ       cosλ  0.0
    #      -cosλ*sinϕ -sinλ*sinϕ cosϕ
    #       cosλ*cosϕ  sinλ*cosϕ sinϕ]
    #
    # east, north, up = R * [∂x, ∂y, ∂z]
    east  = ∂x * -sinλ      + ∂y * cosλ       + ∂z * 0.0
    north = ∂x * -cosλ*sinϕ + ∂y * -sinλ*sinϕ + ∂z * cosϕ
    up    = ∂x * cosλ*cosϕ  + ∂y * sinλ*cosϕ  + ∂z * sinϕ

    return ENU(east, north, up)
end
ENU(ecef::ECEF, lla_ref::LLA, datum) = ENU(ecef, lla_ref, ellipsoid(datum))

##############################
### LLA to ENU coordinates ###
##############################

# Given a reference point for linarization
function ENU(lla::LLA, lla_ref::LLA, datum::Ellipsoid)
    ecef = ECEF(lla, datum)
    return ENU(ecef, lla_ref, datum)
end
ENU(lla::LLA, lla_ref::LLA, datum) = ENU(lla, lla_ref, ellipsoid(datum))
=#
