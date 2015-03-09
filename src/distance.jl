
distance(a::ENU, b::ENU) = distance(a.east, a.north, a.up,
                                    b.east, b.north, b.up)

distance(a::ECEF, b::ECEF) = distance(a.x, a.y, a.z,
                                      b.x, b.y, b.z)

function distance(x1, y1, z1, x2, y2, z2)
    Δx = x2 - x1
    Δy = y2 - y1
    Δz = z2 - z1

    return hypot(hypot(Δx,  Δy), Δz)
end
