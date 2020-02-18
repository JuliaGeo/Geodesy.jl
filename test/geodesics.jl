@testset "Geodesics" begin
    gc = GreatCircle(wgs84)
    datum = wgs84
    ellps = ellipsoid(datum)

    # Construction
    @test gc isa GreatCircle
    @test gc.datum == datum

    # Invocation
    @test_throws UndefKeywordError gc(LLA(0, 0, 0), dist=0)
    @test_throws ArgumentError gc(LLA(0, 0, 0), azi=0)
    @test_throws ArgumentError gc(LLA(0, 0, 0), azi=0, dist=0, angle=0)

    # Output
    @test propertynames(gc(LLA(0, 0, 0), LLA(1, 1, 1))) == (:azi, :baz, :dist, :angle)
    @test propertynames(gc(LLA(0, 0, 0), azi=0, dist=1)) == (:lon, :lat, :baz, :dist, :angle)

    # Conversion between LLA and other types
    point1 = LLA(lat = 180rand() - 90, lon = 360rand(), alt=100_000rand())
    point2 = LLA(lat = 180rand() - 90, lon = 360rand(), alt=100_000rand())
    origin1 = point1
    zone1 = UTMZ(point1, datum).zone
    isnorth1 = point1.lat >= 0
    origin2 = point2
    zone2 = UTMZ(point2, datum).zone
    isnorth2 = point2.lat >= 0
    point_types1 = (identity,
                   x -> ECEF(x, datum),
                   LatLon,
                   x -> ENU(x, origin1, datum),
                   x -> UTMZ(x, datum),
                   x -> UTM(x, zone1, isnorth1, datum))
    point_types2 = (identity,
                   x -> ECEF(x, datum),
                   LatLon,
                   x -> ENU(x, origin2, datum),
                   x -> UTMZ(x, datum),
                   x -> UTM(x, zone2, isnorth2, datum))
    azi = 720rand() - 360
    dist = 100_000rand()
    arc = 720rand()

    for type1 in point_types1
        tpoint1 = type1(point1)
        # Approximate equality because of numerical error in conversion
        # Forward
        if tpoint1 isa ENU
            @test gc(tpoint1, origin1, azi=azi, dist=dist) ≈ gc(point1, azi=azi, dist=dist)
            @test gc(tpoint1, origin1, azi=azi, angle=arc) ≈ gc(point1, azi=azi, angle=arc)
        elseif tpoint1 isa UTM
            @test gc(tpoint1, zone1, isnorth1, azi=azi, angle=dist) ≈ gc(point1, azi=azi, angle=dist)
            @test gc(tpoint1, zone1, isnorth1, azi=azi, angle=arc) ≈ gc(point1, azi=azi, angle=arc)
        else
            @test gc(tpoint1, azi=azi, dist=dist) ≈ gc(point1, azi=azi, dist=dist)
            @test gc(tpoint1, azi=azi, angle=arc) ≈ gc(point1, azi=azi, angle=arc)
        end
        # Inverse
        for type2 in point_types2
            result = gc(point1, point2)
            tpoint2 = type2(point2)
            if tpoint1 isa ENU
                if tpoint2 isa ENU
                    @test gc(tpoint1, origin1, tpoint2, origin2) ≈ result
                elseif tpoint2 isa UTM
                    @test gc(tpoint1, origin1, tpoint2, zone2, isnorth2) ≈ result
                else
                    @test gc(tpoint1, tpoint2, origin1) ≈ result
                end
            elseif tpoint1 isa UTM
                if tpoint2 isa ENU
                    @test gc(tpoint1, zone1, isnorth1, tpoint2, origin2) ≈ result
                elseif tpoint2 isa UTM
                    @test gc(tpoint1, zone1, isnorth1, tpoint2, zone2, isnorth2) ≈ result
                else
                    @test gc(tpoint1, tpoint2, zone1, isnorth1) ≈ result
                end
            else
                if tpoint2 isa ENU
                    @test gc(tpoint1, tpoint2, origin2) ≈ result
                elseif tpoint2 isa UTM
                    @test gc(tpoint1, tpoint2, zone2, isnorth2) ≈ result
                else
                    @test gc(tpoint1, tpoint2) ≈ result
                end
            end
        end
    end
end
