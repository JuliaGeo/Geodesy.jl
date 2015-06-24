# Vicenty's inverse (two points -> ellipsoidal distance and azimuths)
const ANTIPODAL_ERROR =
    ArgumentError("Points too near antipodal for distance algorithm to converge")

function vicentys_inverse{T <: LL}(a::T, b::T)
    d = ellipsoid(T)
    f = 1 - d.b/d.a

    # Reduced latitude (latitude on the auxiliary sphere)
    U₁ = atan(d.b / d.a * tand(a.lat))
    U₂ = atan(d.b / d.a * tand(b.lat))

    sinU₁, cosU₁ = sin(U₁), cos(U₁)
    sinU₂, cosU₂ = sin(U₂), cos(U₂)

    cosU₁cosU₂ = cosU₁ * cosU₂
    sinU₁sinU₂ = sinU₁ * sinU₂
    sinU₁cosU₂ = sinU₁ * cosU₂
    cosU₁sinU₂ = cosU₁ * sinU₂

    ΔL = deg2rad(b.lon - a.lon)

    # Needed in scope after the loop
    oldλ = sinλ = cosλ = sinσ = cosσ = σ = cos²α = cos2σm = cos²2σm = Inf
    itr = 0

    # Initialize λ to ΔL for first iteration
    λ = ΔL

    for itr=1:200
        oldλ = λ

        sinλ, cosλ = sin(λ), cos(λ)

        x = cosU₂ * sinλ
        y = cosU₁sinU₂ - sinU₁cosU₂ * cosλ
        sinσ = sqrt(x*x + y*y)
        cosσ = sinU₁sinU₂ + cosU₁cosU₂ * cosλ
        σ = atan2(sinσ, cosσ)

        sinα = cosU₁cosU₂ * sinλ / sinσ
        cos²α = 1 - sinα*sinα
        cos2σm = cosσ - 2sinU₁sinU₂ / cos²α
        cos²2σm = cos2σm * cos2σm
        C = f / 16 * cos²α * (4 + f*(4-3cos²α))

        λ = ΔL + (1-C) * f * sinα *
            (σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2cos²2σm)))

        abs(λ - oldλ) < 1e-12 && break
    end

    # Failure to converge happens when abs(λ) > π after the first iteration
    itr == 200 && throw(ANTIPODAL_ERROR)

    a², b² = d.a * d.a, d.b * d.b
    u² = cos²α * (a² - b²) / b²
    A = 1 + u² / 16384 * (4096 + u² * (-768 + u² * (320 - 175u²)))
    B = u² / 1024 * (256 + u² * (-128 + u² * (74 - 47u²)))
    Δσ = B * sinσ * (cos2σm + B / 4 *
        (cosσ * (-1 + 2cos²2σm) -
         B / 6 * cos2σm * (-3 + 4sinσ*sinσ) * (-3 + 4cos²2σm)))

    s = d.b * A * (σ - Δσ) # ellipsoidal distance
    α₁ = atan2(cosU₂ * sinλ,  cosU₁sinU₂ - sinU₁cosU₂ * cosλ) # start azimuth
    α₂ = atan2(cosU₁ * sinλ, -sinU₁cosU₂ + cosU₁sinU₂ * cosλ) # end azimuth

    return s, rad2deg(α₁), rad2deg(α₂)
end
