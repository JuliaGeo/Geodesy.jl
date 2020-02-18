# Displaying custom types

const SHOW_MAXLEN = maximum(x -> first(x) != '_' ? length(x) : 0,
    String(f) for T in (Geodesic, GeodesicLine, Result, Polygon) for f in fieldnames(T))

function Base.show(io::IO, x::Union{Geodesic, GeodesicLine, Result, Polygon})
    println(io, typeof(x), ":")
    for f in fieldnames(typeof(x))
    	first(String(f)) == '_' && continue
        v = getfield(x, f)
        println(io, lpad(String(f), SHOW_MAXLEN), ": ", v === nothing ? "<>" : v)
    end
end
