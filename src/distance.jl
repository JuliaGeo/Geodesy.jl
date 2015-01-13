distance(a::ENU, b::ENU) = distance(a.east, a.north, a.up,
                                    b.east, b.north, b.up)

distance(a::ECEF, b::ECEF) = distance(a.x, a.y, a.z,
                                      b.x, b.y, b.z)

# Cartesian coordinates
function distance(x0, y0, z0, x1, y1, z1)
    return sqrt((x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2)
end
