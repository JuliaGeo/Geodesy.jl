
# cartesian distance

distance(a::ENU, b::ENU) = _distance(a.east, a.north, a.up,
                                     b.east, b.north, b.up)

distance(a::ECEF, b::ECEF) = _distance(a.x, a.y, a.z,
                                       b.x, b.y, b.z)

function _distance(x1, y1, z1, x2, y2, z2)
    Δx = x2 - x1
    Δy = y2 - y1
    Δz = z2 - z1

    return hypot(hypot(Δx,  Δy), Δz)
end

# haversine (spherical)

function haversine_distance{T <: LL}(a::T, b::T)
    r = 6_371_009 # reduces average rather than maximum error

    # Using deg2rad is faster and a tiny bit less accurate than sind/cosd --
    #   already very rough, might as well keep it fast.
    ϕ₁, λ₁ = deg2rad(a.lat), deg2rad(a.lon)
    ϕ₂, λ₂ = deg2rad(b.lat), deg2rad(b.lon)

    hϕ = haversin(ϕ₂ - ϕ₁)
    hλ = haversin(λ₂ - λ₁)

    h = hϕ + cos(ϕ₁)*cos(ϕ₂)*hλ

    return 2r * asin(sqrt(h))
end

function haversin(θ)
    s = sin(θ/2)
    return s*s
end

# vicenty's (ellipsoidal)

distance{T <: LL}(a::T, b::T) = vicentys_inverse(a, b)[1]
