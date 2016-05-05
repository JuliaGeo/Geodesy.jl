@testset "UTM functions" begin
    # Test zone prediction
    # TODO: test Norway, Svalbard
    @test Geodesy.utm_zone(1.0,1.0) == (31, true)
    @test Geodesy.utm_zone(-1.0,-1.0) == (30, false)

    @test Geodesy.utm_meridian(31) == 3
end
