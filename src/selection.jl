export RandomStrategy, FirstStrategy, LastStrategy, DegreeStrategy, CodegreeStrategy,
  NormalStrategy, StrangeStrategy, SugarStrategy, SpiceStrategy, TrueDegree

# TODO: use cmp rather than <

abstract type SelectionStrategy end

@doc raw"""
    select(S::SelectionStrategy, G::Vector{T}, Q::Vector{Tuple{Int,Int}}, o::MonomialOrdering) where {T<:MPolyRingElem}

Select a critical pair from $Q$ according to the provided strategy.
Return the index of the selected pair in $Q$.
"""
function select end

struct RandomStrategy <: SelectionStrategy
  random::Xoshiro
  RandomStrategy(seed::Int) = new(Xoshiro(seed))
end

function select(S::RandomStrategy, _, Q::Vector{Tuple{Int,Int}}, _)
  @assert !isempty(Q)
  inc!(STATS, :critical_pairs_selected)
  return rand(S.random, eachindex(Q))
end

abstract type TieStrategy <: SelectionStrategy end

function select(
  S::TieStrategy, G::Vector{T}, Q::Vector{Tuple{Int,Int}}, o::MonomialOrdering
) where {T<:MPolyRingElem}
  @assert !isempty(Q)
  inc!(STATS, :critical_pairs_selected)

  i = 1
  for j in 2:length(Q)
    if compare(S, G, Q[j]..., Q[i]..., o)
      i = j
    end
  end
  return i
end

struct FirstStrategy <: TieStrategy end

function compare(::FirstStrategy, _, i1::Int, j1::Int, i2::Int, j2::Int, _)
  return (j1, i1) < (j2, i2)
end

struct LastStrategy <: TieStrategy end

function compare(::LastStrategy, _, i1::Int, j1::Int, i2::Int, j2::Int, _)
  return (j1, i1) > (j2, i2)
end

struct DegreeStrategy <: TieStrategy end

function compare(
  ::DegreeStrategy, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # compute degree of lcm(LM(G[i]), LM(G[j]))
  d1 = sum(max.(lx(G[i1], o), lx(G[j1], o)))
  d2 = sum(max.(lx(G[i2], o), lx(G[j2], o)))
  return (d1, j1, i1) < (d2, j2, i2)
end

struct CodegreeStrategy <: TieStrategy end

function compare(
  ::CodegreeStrategy, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # compute degree of lcm(LM(G[i]), LM(G[j]))
  d1 = sum(max.(lx(G[i1], o), lx(G[j1], o)))
  d2 = sum(max.(lx(G[i2], o), lx(G[j2], o)))
  return (d1, j1, i1) > (d2, j2, i2)
end

struct NormalStrategy <: TieStrategy end

function compare(
  ::NormalStrategy, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # compute lcm(LM(G[i]), LM(G[j]))
  l1 = lcm(lm(G[i1], o), lm(G[j1], o))
  l2 = lcm(lm(G[i2], o), lm(G[j2], o))
  return (l1, j1, i1) < (l2, j2, i2)
end

struct StrangeStrategy <: TieStrategy end

function compare(
  ::StrangeStrategy, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # compute lcm(LM(G[i]), LM(G[j]))
  l1 = lcm(lm(G[i1], o), lm(G[j1], o))
  l2 = lcm(lm(G[i2], o), lm(G[j2], o))
  return (l1, j1, i1) > (l2, j2, i2)
end

# TODO: implement sugar, add a hint to glaze f
sugar(f::MPolyRingElem) = error("$f has no sugar attached.")

struct SugarStrategy <: TieStrategy end

function compare(
  ::SugarStrategy, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # compute degree of lcm(LM(G[i]), LM(G[j]))
  d1 = sum(max.(lx(G[i1], o), lx(G[j1], o)))
  d2 = sum(max.(lx(G[i2], o), lx(G[j2], o)))
  # compute sugar of spoly(G[i], G[j])
  s1 = max(d1 - sum(lx(G[i1], o)) + sugar(G[i1]), d1 - sum(lx(G[j1], o)) + sugar(G[j1]))
  s2 = max(d2 - sum(lx(G[i2], o)) + sugar(G[i2]), d2 - sum(lx(G[j2], o)) + sugar(G[j2]))
  # compute lcm(LM(G[i]), LM(G[j]))
  l1 = lcm(lm(G[i1], o), lm(G[j1], o))
  l2 = lcm(lm(G[i2], o), lm(G[j2], o))
  return (s1, l1, j1, i1) < (s2, l2, j2, i2)
end

struct SpiceStrategy <: TieStrategy end

function compare(
  ::SpiceStrategy, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # compute degree of lcm(LM(G[i]), LM(G[j]))
  d1 = sum(max.(lx(G[i1], o), lx(G[j1], o)))
  d2 = sum(max.(lx(G[i2], o), lx(G[j2], o)))
  # compute sugar of spoly(G[i], G[j])
  s1 = max(d1 - sum(lx(G[i1], o)) + sugar(G[i1]), d1 - sum(lx(G[j1], o)) + sugar(G[j1]))
  s2 = max(d2 - sum(lx(G[i2], o)) + sugar(G[i2]), d2 - sum(lx(G[j2], o)) + sugar(G[j2]))
  # compute lcm(LM(G[i]), LM(G[j]))
  l1 = lcm(lm(G[i1], o), lm(G[j1], o))
  l2 = lcm(lm(G[i2], o), lm(G[j2], o))
  return (s1, l1, j1, i1) > (s2, l2, j2, i2)
end

struct TrueDegree <: TieStrategy end

function compare(
  ::TrueDegree, G::Vector{T}, i1::Int, j1::Int, i2::Int, j2::Int, o::MonomialOrdering
) where {T<:MPolyRingElem}
  # prefer pairs where the degree of lcm(LM(G[i]), LM(G[j])) is small
  d1 = sum(max.(lx(G[i1], o), lx(G[j1], o)))
  d2 = sum(max.(lx(G[i2], o), lx(G[j2], o)))
  d1 < d2 && return true
  d1 > d2 && return false

  # prefer pairs where the number of monomials in spoly(G[i], G[j]) is small
  s1 = spoly(G[i1], G[j1])
  s2 = spoly(G[i2], G[j2])
  l1 = length(s1)
  l2 = length(s2)
  l1 < l2 && return true
  l1 > l2 && return false

  # prefer pairs where the degree of spoly(G[i], G[j]) is small
  d1 = degree(s1, o)
  d2 = degree(s2, o)
  d1 < d2 && return true
  d1 > d2 && return false

  return false
end
