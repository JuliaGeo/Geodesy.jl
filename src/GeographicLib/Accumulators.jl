module Accumulators

import ..Math

"""Like math.fsum, but allows a running sum"""
mutable struct Accumulator{T<:AbstractFloat}
    _s::T
    _t::T
    function Accumulator(_s::T1, _t::T2) where {T1,T2}
    	T = promote_type(float(T1), float(T2))
    	Accumulator{T}(T(_s), T(_t))
    end
    Accumulator{T}(_s, _t) where T = new{T}(T(_s), T(_t))
end

"""Constructor"""
Accumulator(y=0.0) = Accumulator(y, 0.0)
Accumulator(acc::Accumulator) = Accumulator(acc._s, acc._t)

Base.:(==)(a::Accumulator, b::Accumulator) = a._s == b._s && a._t == b._t

"""Set value from argument"""
Set!(self::Accumulator, y::Accumulator) = ((self._s, self._t) = (y._s, y._t); self)
Set!(self::Accumulator, y) = ((self._s, self._t) = (y, 0); self)


"""Add a value"""
function Add!(self, y)
  # Here's Shewchuk's solution...
  # hold exact sum as [s, t, u]
  y, u = Math.sum(y, self._t)             # Accumulate starting at
  self._s, self._t = Math.sum(y, self._s) # least significant end
  # Start is _s, _t decreasing and non-adjacent.  Sum is now (s + t + u)
  # exactly with s, t, u non-adjacent and in decreasing order (except
  # for possible zeros).  The following code tries to normalize the
  # result.  Ideally, we want _s = round(s+t+u) and _u = round(s+t+u -
  # _s).  The follow does an approximate job (and maintains the
  # decreasing non-adjacent property).  Here are two "failures" using
  # 3-bit floats:
  #
  # Case 1: _s is not equal to round(s+t+u) -- off by 1 ulp
  # [12, -1] - 8 -> [4, 0, -1] -> [4, -1] = 3 should be [3, 0] = 3
  #
  # Case 2: _s+_t is not as close to s+t+u as it shold be
  # [64, 5] + 4 -> [64, 8, 1] -> [64,  8] = 72 (off by 1)
  #                    should be [80, -7] = 73 (exact)
  #
  # "Fixing" these problems is probably not worth the expense.  The
  # representation inevitably leads to small errors in the accumulated
  # values.  The additional errors illustrated here amount to 1 ulp of
  # the less significant word during each addition to the Accumulator
  # and an additional possible error of 1 ulp in the reported sum.
  #
  # Incidentally, the "ideal" representation described above is not
  # canonical, because _s = round(_s + _t) may not be true.  For
  # example, with 3-bit floats:
  #
  # [128, 16] + 1 -> [160, -16] -- 160 = round(145).
  # But [160, 0] - 16 -> [128, 16] -- 128 = round(144).
  #
  if self._s == 0            # This implies t == 0,
    self._s = u               # so result is u
  else
    self._t += u              # otherwise just accumulate u to t.
  end
  self
end

"""Return sum + y"""
function Sum(self, y = 0.0)
  if y == 0.0
    return self._s
  else
    b = Accumulator(self)
    Add!(b, y)
    return b._s
  end
end

"""Negate sum"""
function Negate!(self)
  self._s *= -1
  self._t *= -1
  self
end

end # module
