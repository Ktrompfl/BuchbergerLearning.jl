
@doc raw"""
    li(f::MPolyRingElem, o::MonomialOrdering)

Return the index of the leading term of $f$ w.r.t. to the monomial ordering $o$.
"""
function li(f::MPolyRingElem, o::MonomialOrdering)
  @assert !iszero(f)
  return Oscar.index_of_leading_term(f, o)
end

@doc raw"""
    lc(f::MPolyRingElem, o::MonomialOrdering)

Return the leading coefficient of $f$ w.r.t. to the monomial ordering $o$.
"""
lc(f::MPolyRingElem, o::MonomialOrdering) = coeff(f, li(f, o))

@doc raw"""
    lm(f::MPolyRingElem, o::MonomialOrdering)

Return the leading monomial of $f$ w.r.t. to the monomial ordering $o$.
"""
lm(f::MPolyRingElem, o::MonomialOrdering) = monomial(f, li(f, o))

@doc raw"""
    lt(f::MPolyRingElem, o::MonomialOrdering)

Return the leading term of $f$ w.r.t. to the monomial ordering $o$.
"""
lt(f::MPolyRingElem, o::MonomialOrdering) = term(f, li(f, o))

@doc raw"""
    lx(f::MPolyRingElem, o::MonomialOrdering)

Return the exponent vector of the leading term of $f$ w.r.t. to the monomial ordering $o$.
"""
lx(f::MPolyRingElem, o::MonomialOrdering) = exponent_vector(f, li(f, o))

@doc raw"""
    degree(f::MPolyRingElem)

Return the degree of $LM(f)$.
"""
function degree(f::MPolyRingElem, o::MonomialOrdering)
  @assert !iszero(f)
  k = li(f, o)
  n = ngens(parent(f))
  return sum(exponent(f, k, i) for i in 1:n)
end

@doc raw"""
    iscoprime(f1::T, f2::T, o::MonomialOrdering) where {T<:MPolyRingElem}

Return whether the lead monomials of $f1$ and $f2$ are coprime.
"""
function iscoprime(f1::T, f2::T, o::MonomialOrdering) where {T<:MPolyRingElem}
  @assert !iszero(f1)
  @assert !iszero(f2)

  k = li(f1, o)
  l = li(f2, o)
  n = ngens(parent(f1))

  return all(exponent(f1, k, i) * exponent(f2, l, i) == 0 for i in 1:n)
end

@doc raw"""
    ismonic(f::MPolyRingElem, o::MonomialOrdering)

Return whether $f$ is monic, i.e. $LC(f) = 1$.
"""
ismonic(f::MPolyRingElem, o::MonomialOrdering) = !iszero(f) && isone(lc(f, o))

@doc raw"""
    reduces(f::T, g::T, o::MonomialOrdering) where {T<:MPolyRingElem}

Return whether $f$ reduces $g$, i.e. if $LM(f) | LM(g)$.
"""
function reduces(f::T, g::T, o::MonomialOrdering) where {T<:MPolyRingElem}
  @assert !iszero(f)
  @assert !iszero(g)
  k = li(f, o)
  l = li(g, o)
  n = ngens(parent(f))

  return all(exponent(f, k, i) <= exponent(g, l, i) for i in 1:n)
end

@doc raw"""
    normalize(f::T, o::MonomialOrdering) where {T<:MPolyRingElem}

Normalize non-zero $f$. Return the resulting monic polynomial.
"""
function normalize(f::T, o::MonomialOrdering) where {T<:MPolyRingElem}
  @assert !iszero(f)
  return inv(lc(f, o)) * f
end
