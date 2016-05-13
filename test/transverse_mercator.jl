@testset "Transverse Mercator projection and inverse" begin

    # File has 5000 entries of lon between ±50°, lat between ±90° (with UTM scaling 0.9996)
    # Adapted from Charles Karney's test data (http://zenodo.org/record/32470#.VzF-Krp97CI)
    f = open("TMcoords_lite.bin") # Binary file, Float64: lat, lon, x, y, γ, k
    tmdat = read(f, Float64, (5000, 6))
    close(f)

    tm = Geodesy.TransverseMercator(wgs84)

    for i = 1:size(tmdat,1)
        lat = tmdat[i,1]
        lon = tmdat[i,2]
        x = tmdat[i,3]
        y = tmdat[i,4]
        γ = tmdat[i,5]
        k = tmdat[i,6]

        (x2, y2, γ2, k2) = Geodesy.transverse_mercator_forward(0.0, lat, lon, 0.9996, tm)

        @test (x2, y2, γ2, k2) ≈ (x, y, γ, k)

        (lat2, lon2, γ2, k2) = Geodesy.transverse_mercator_reverse(0.0, x, y, 0.9996, tm)

        @test (lat2, lon2, γ2, k2) ≈ (lat, lon, γ, k)
    end
end
