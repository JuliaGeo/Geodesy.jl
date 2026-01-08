@testset "Point constructors" begin
    # Constructors
    lla = LLA(1., 1., 0.)
    @test LLA(1., 1.) == lla
    @test LLA(lat = 1., lon = 1., alt = 0.) == lla
    @test LLA(1//1, 1, 0.) == lla
    @test LLA(1, 1) == LLA{Int}(1, 1, 0)

    latlon = LatLon(1., 1.)
    @test LatLon(lat = 1., lon = 1.) == latlon
    @test LatLon(ECEF(lla, wgs84), wgs84) ≈ latlon
    @test LatLon(lla) == latlon
    @test LatLon(1//1, 1.) == latlon
    @test LatLon(1, 1) == LatLon{Int}(1, 1)
    @test_throws MethodError LatLon(missing, missing)

    @test LLA(latlon) == lla

    utm = UTM(10., 10., 0.)
    @test UTM(10., 10.) == utm

    utmz = UTMZ(10., 10., 0., 1, true, 'D')
    @test UTMZ(10., 10., 1, true, 'D') == utmz
    @test UTMZ(utm, 1, true, 'D') == utmz

    @test UTM(utmz) == utm
end # @testset
