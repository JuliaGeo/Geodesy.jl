@testset "Datums and ellipsoids" begin
    # Check the ellipsoids
    @test ellipsoid(wgs84) == Geodesy.wgs84_el
    @test ellipsoid(grs80) == Geodesy.grs80_el
    @test ellipsoid(nad27) == Geodesy.nad27_el
    @test ellipsoid(osgb36) == Geodesy.osgb36_el

    # Check transverse-Mercator pre-calculations
    @test Geodesy.TransverseMercator(wgs84) == Geodesy.wgs84_tm
    @test Geodesy.TransverseMercator(grs80) == Geodesy.grs80_tm
    @test Geodesy.TransverseMercator(nad27) == Geodesy.nad27_tm
    @test Geodesy.TransverseMercator(osgb36) == Geodesy.osgb36_tm

    # The LLAfromECEF aren't pre-cached, as yet. Probably should just create a
    # large ellispoid thingie for all of them to share...
    @test (LLAfromECEF(wgs84); true)
    @test (LLAfromECEF(grs80); true)
    @test (LLAfromECEF(nad27); true)
    @test (LLAfromECEF(osgb36); true)
end
