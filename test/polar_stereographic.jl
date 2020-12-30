@testset "Polar stereographic" begin

    tm = Geodesy.TransverseMercator(wgs84)

    # Exact north pole
    lat = 90.0
    lon = 180.0
    x = 2000000.0 - 2e6
    y = 2000000.0 - 2e6
    northpole = true

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    # Data generated independely from Charles Karney's http://geographiclib.sourceforge.net/cgi-bin/GeoConvert?input=89+1&zone=0&prec=9&option=Submit

    lat = 89.0
    lon = 1.0
    x = 2001937.679954440 - 2e6
    y = 1888990.389757469 - 2e6
    northpole = true

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = 88.0
    lon = 91.0
    x = 2222035.446703544 - 2e6
    y = 2003875.643138576 - 2e6
    northpole = true

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = 87.0
    lon = -179.0
    x = 1994185.827037653 - 2e6
    y = 2333093.745927457 - 2e6
    northpole = true

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = 86.0
    lon = -89.0
    x = 1555799.234876368 - 2e6
    y = 1992246.446803603 - 2e6
    northpole = true

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = -89.0
    lon = 1.0
    x = 2001937.679954440 - 2e6
    y = 2111009.610242531 - 2e6
    northpole = false

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = -88.0
    lon = 91.0
    x = 2222035.446703544 - 2e6
    y = 1996124.356861424 - 2e6
    northpole = false

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = -87.0
    lon = -179.0
    x = 1994185.827037653 - 2e6
    y = 1666906.254072543 - 2e6
    northpole = false

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    lat = -86.0
    lon = -89.0
    x = 1555799.234876368 - 2e6
    y = 2007753.553196397 - 2e6
    northpole = false

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)

    # Exact south pole
    lat = -90.0
    lon = 0.0
    x = 2000000.0 - 2e6
    y = 2000000.0 - 2e6
    northpole = false

    (x2, y2, γ, k) = Geodesy.polarst_fwd(northpole, 0.994, tm, lat, lon)
    @test (x2, y2) ≈ (x, y)

    (lat2, lon2, γ, k) = Geodesy.polarst_inv(northpole, 0.994, tm, x, y)
    @test (lat2, lon2) ≈ (lat, lon)
end
