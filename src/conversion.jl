
###############################
### LLA to ECEF coordinates ###
###############################

function ECEF(lla::LLA)
    lat, lon, alt = lla.lat, lla.lon, lla.alt
    d = WGS84
    e² = d.e * d.e

    sinlat, coslat = sind(lat), cosd(lat)


    N = d.a / sqrt(1 - e² * sinlat^2)  # Radius of curvature (meters)

    x = (N + alt) * coslat * cosd(lon)
    y = (N + alt) * coslat * sind(lon)
    z = (N * (1 - e²) + alt) * sinlat

    return ECEF(x, y, z)
end

###############################
### ECEF to LLA coordinates ###
###############################

function LLA(ecef::ECEF)
    x, y, z = ecef.x, ecef.y, ecef.z
    d = WGS84
    e² = d.e * d.e

    p = sqrt(x*x + y*y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e_prime^2 * d.b * sin(θ)^3, p - e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    return LLA(ϕ*180/π, λ*180/π, h)
end

###############################
### ECEF to ENU coordinates ###
###############################

# Given a reference point for linarization
function ENU(ecef::ECEF, lla_ref::LLA)
    ϕ = lla_ref.lat
    λ = lla_ref.lon

    ecef_ref = ECEF(lla_ref)
    ∂x = ecef.x - ecef_ref.x
    ∂y = ecef.y - ecef_ref.y
    ∂z = ecef.z - ecef_ref.z

    # Compute rotation matrix
    sinλ, cosλ = sind(λ), cosd(λ)
    sinϕ, cosϕ = sind(ϕ), cosd(ϕ)

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

# Given Bounds object for linearization
function ENU(ecef::ECEF, bounds::Bounds{LLA})
    lla_ref = center(bounds)
    return ENU(ecef, lla_ref)
end

##############################
### LLA to ENU coordinates ###
##############################

# Given a reference point for linarization
function ENU(lla::LLA, lla_ref::LLA)
    ecef = ECEF(lla)
    return ENU(ecef, lla_ref)
end

# Given Bounds object for linearization
function ENU(lla::LLA, bounds::Bounds{LLA})
    ecef = ECEF(lla)
    return ENU(ecef, bounds)
end

#################################
### LLA to ENU Bounds objects ###
#################################

function ENU(bounds::Bounds{LLA}, lla_ref::LLA = center(bounds))
    top_left_LLA = LLA(bounds.max_y, bounds.min_x)
    bottom_right_LLA = LLA(bounds.min_y, bounds.max_x)

    top_left_ENU = ENU(top_left_LLA, lla_ref)
    bottom_right_ENU = ENU(bottom_right_LLA, lla_ref)

    return Bounds{ENU}(bottom_right_ENU.north,
                       top_left_ENU.north,
                       top_left_ENU.east,
                       bottom_right_ENU.east)
end

########################
### Helper Functions ###
########################

### Get center point of Bounds region ###
function center{T}(bounds::Bounds{T})
    y_ref = (bounds.min_y + bounds.max_y) / 2
    x_ref = (bounds.min_x + bounds.max_x) / 2

    return T(XY(x_ref, y_ref))
end
