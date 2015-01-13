
###############################
### LLA to ECEF coordinates ###
###############################

function ECEF(lla::LLA)
    lat, lon, alt = lla.lat, lla.lon, lla.alt
    d = WGS84

    N = d.a / sqrt(1 - d.e*d.e * sind(lat)^2)  # Radius of curvature (meters)

    x = (N + alt) * cosd(lat) * cosd(lon)
    y = (N + alt) * cosd(lat) * sind(lon)
    z = (N * (1 - d.e*d.e) + alt) * sind(lat)

    return ECEF(x, y, z)
end

###############################
### ECEF to LLA coordinates ###
###############################

function LLA(ecef::ECEF)
    x, y, z = ecef.x, ecef.y, ecef.z
    d = WGS84

    p = sqrt(x*x + y*y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e_prime^2 * d.b * sin(θ)^3, p - d.e*d.e*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e*d.e * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    return LLA(ϕ*180/π, λ*180/π, h)
end

###############################
### ECEF to ENU coordinates ###
###############################

# Given a reference point for linarization
function ENU(ecef::ECEF, lla_ref::LLA)
    # Reference point to linearize about
    ϕ = lla_ref.lat
    λ = lla_ref.lon

    ecef_ref = ECEF(lla_ref)
    ecef_vec = [ecef.x - ecef_ref.x; ecef.y - ecef_ref.y; ecef.z - ecef_ref.z]

    # Compute rotation matrix
    R = [-sind(λ) cosd(λ) 0;
         -cosd(λ)*sind(ϕ) -sind(λ)*sind(ϕ) cosd(ϕ);
         cosd(λ)*cosd(ϕ) sind(λ)*cosd(ϕ) sind(ϕ)]
    ned = R * ecef_vec

    # Extract elements from vector
    east = ned[1]
    north = ned[2]
    up = ned[3]

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
