using Geodesy
using Base.Test
using Geodesy: XY, onBounds, boundaryPoint

# Construction

x, y = (rand(2) - .5) * 10_000

@test_throws TypeError Bounds{ECEF}(y, y + 1, x, x + 1)

Bounds{ENU}(y, y + 1, x, x + 1)

x, y = (rand() - .5) * 360, (rand() - .5) * 180

Bounds{LLA}(y, y + 1, x, x + 1)

Bounds(y, y, x, x)
Bounds(y, y, x, x - 1)
@test_throws ArgumentError Bounds(y, y - 1, x, x)

@test_throws ArgumentError Bounds(-90.1, 90.1, 0, 1)
@test_throws ArgumentError Bounds(0, 1, -180.1, 180.1)

# center

@test center(Bounds(0, 0, 179, -178)) == LLA(0, -179.5)

# inBounds

function test_bounds{T}(bounds::Bounds{T})
    min_x, min_y, max_x, max_y = bounds.min_x, bounds.min_y, bounds.max_x, bounds.max_y

    XY = Geodesy.XY

    @test inBounds(T(XY(min_x, min_y)), bounds)
    @test !inBounds(T(XY(min_x - eps(min_x), min_y)), bounds)
    @test !inBounds(T(XY(min_x, min_y - eps(min_y))), bounds)

    @test inBounds(T(XY(max_x, max_y)), bounds)
    @test !inBounds(T(XY(max_x + eps(max_x), max_y)), bounds)
    @test !inBounds(T(XY(max_x, max_y + eps(max_y))), bounds)
end

test_bounds(Bounds{ENU}(1.1, 2.2, 3.3, 4.4))
test_bounds(Bounds(1.1, 2.2, 3.3, 4.4))
test_bounds(Bounds(1.1, 2.2, 179.3, -179.4))

# boundaryPoint

for bounds in (Bounds(0, 1, 78, 79), Bounds(0, 1, 179.9, -179.1), Bounds{ENU}(-1, 1, -1, 1))
    T = isa(bounds, Bounds{LLA}) ? LLA : ENU

    c = center(bounds)
    cx, cy = getX(c), getY(c)

    for _ = 1:1_000
        in_both =    T(XY(cx + rand() - 0.5, cy + rand() - 0.5))
        in_x1 =      T(XY(cx + rand() - 0.5, cy + rand() + 1.0))
        in_x2 =      T(XY(cx + rand() - 0.5, cy - rand() - 1.0))
        in_y =       T(XY(cx + rand() + 1.0, cy + rand() - 0.5))
        in_neither = T(XY(cx + rand() + 1.0, cy + rand() + 1.0))

        @test onBounds(boundaryPoint(in_both, in_x1, bounds), bounds)
        @test onBounds(boundaryPoint(in_x2, in_both, bounds), bounds)
        @test onBounds(boundaryPoint(in_both, in_y, bounds), bounds)
        @test onBounds(boundaryPoint(in_neither, in_both, bounds), bounds)

    end
end
