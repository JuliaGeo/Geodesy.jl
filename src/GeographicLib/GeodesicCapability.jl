module GeodesicCapability

const CAP_NONE = 0
const CAP_C1   = 1 << 0
const CAP_C1p  = 1 << 1
const CAP_C2   = 1 << 2
const CAP_C3   = 1 << 3
const CAP_C4   = 1 << 4
const CAP_ALL  = Int(0x1F)
const CAP_MASK = CAP_ALL
const OUT_ALL  = Int(0x7F80)
const OUT_MASK = Int(0xFF80)             # Includes LONG_UNROLL
const EMPTY         = 0
const LATITUDE      = 1 << 7  | CAP_NONE
const LONGITUDE     = 1 << 8  | CAP_C3
const AZIMUTH       = 1 << 9  | CAP_NONE
const DISTANCE      = 1 << 10 | CAP_C1
const STANDARD      = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE
const DISTANCE_IN   = 1 << 11 | CAP_C1 | CAP_C1p
const REDUCEDLENGTH = 1 << 12 | CAP_C1 | CAP_C2
const GEODESICSCALE = 1 << 13 | CAP_C1 | CAP_C2
const AREA          = 1 << 14 | CAP_C4
const LONG_UNROLL   = 1 << 15
const LONG_NOWRAP   = LONG_UNROLL       # For backwards compatibility only
const ALL           = OUT_ALL | CAP_ALL

end # module
