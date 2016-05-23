@testset "UTM functions" begin
    # Test zone prediction
    @test Geodesy.utm_zone(1.0, 1.0) == (31, true)
    @test Geodesy.utm_zone(-1.0, -1.0) == (30, false)
    @test Geodesy.utm_zone(89.0, -1.0) == (0, true)
    @test Geodesy.utm_zone(-89.0, 1.0) == (0, false)

    # Svalbard
    @test Geodesy.utm_zone(75.0, 8.0) == (31, true)
    @test Geodesy.utm_zone(83.9, 10.0) == (33, true)
    @test Geodesy.utm_zone(72.1, 20.0) == (33, true)
    @test Geodesy.utm_zone(75.0, 22.0) == (35, true)
    @test Geodesy.utm_zone(75.0, 32.0) == (35, true)
    @test Geodesy.utm_zone(75.0, 34.0) == (37, true)

    # Noway
    @test Geodesy.utm_zone(60.0, 2.9) == (31, true)
    @test Geodesy.utm_zone(60.0, 3.1) == (32, true)

    # From types
    @test Geodesy.utm_zone(LatLon(1.0, 1.0)) == (31, true)
    @test Geodesy.utm_zone(LLA(1.0, 1.0)) == (31, true)
    @test Geodesy.utm_zone(ECEF(LLA(1.0, 1.0), wgs84), wgs84) == (31, true)

    # central meridians
    @test Geodesy.utm_meridian(31) == 3
    @test_throws Exception Geodesy.utm_meridian(0)
    @test_throws Exception Geodesy.utm_meridian(61)

end
