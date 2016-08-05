@testset "Datums and ellipsoids" begin
    # Check the ellipsoids
    @test ellipsoid(WGS84) == wgs84_ellipsoid
    @test ellipsoid(NAD27) == clarke1866
    @test ellipsoid(NAD83) == grs80
    @test ellipsoid(OSGB36) == airy1830
    @test ellipsoid(GDA94) == grs80

    # Check transverse-Mercator pre-calculations
    @test Geodesy.TransverseMercator(WGS84) == Geodesy.wgs84_tm
    @test Geodesy.TransverseMercator(NAD27) == Geodesy.clarke1866_tm
    @test Geodesy.TransverseMercator(OSGB36) == Geodesy.airy1830_tm
    @test Geodesy.TransverseMercator(GDA94) == Geodesy.grs80_tm

    # The LLAfromECEF aren't pre-cached, as yet. Probably should just create a
    # large ellispoid thingie for all of them to share...
    @test (LLAfromECEF(wgs84); true)
    @test (LLAfromECEF(nad27); true)
    @test (LLAfromECEF(nad83); true)
    @test (LLAfromECEF(osgb36); true)
    @test (LLAfromECEF(gda94); true)
end
