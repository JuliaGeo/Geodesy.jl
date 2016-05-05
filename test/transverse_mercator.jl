@testset "Transverse Mercator projection and inverse" begin

    @test Geodesy.transverse_mercator(0.0, 3.0, 1.0, 1.0, Geodesy.ellipsoid(wgs84)) ≈ (-222650.79679758303,110642.22941557941,-0.03491928038359307,1.0006134677755945)

    @test Geodesy.transverse_mercator_inv(0.0, 3.0, -222650.79679758303, 110642.22941557941, Geodesy.ellipsoid(wgs84)) ≈ (1.0,1.0,-0.03491928038359307,1.0006134677755945)

end
