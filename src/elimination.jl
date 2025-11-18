export NoStrategy, ProductCriterion, GebauerMoellerInstallation

abstract type EliminationStrategy end

@doc raw"""
    update!(S::EliminationStrategy, Q::Vector{Tuple{Int,Int}}, G::Vector{T}, f::T, o::MonomialOrdering) where {T<:MPolyRingElem}

Update the critical pairs $Q$ when adding $f$ to $G$.
Eliminate critical pairs according to the provided strategies.
"""
function update! end

struct NoStrategy <: EliminationStrategy end

function update!(
  ::NoStrategy, Q::Vector{Tuple{Int,Int}}, G::Vector{T}, f::T, o::MonomialOrdering
) where {T<:MPolyRingElem}
  @assert !iszero(f)
  k = length(G) + 1

  for i in eachindex(G)
    push!(Q, (i, k))
  end

  # stats
  inc!(STATS, :critical_pairs_created, length(G))
end

struct ProductCriterion <: EliminationStrategy end

function update!(
  ::ProductCriterion, Q::Vector{Tuple{Int,Int}}, G::Vector{T}, f::T, o::MonomialOrdering
) where {T<:MPolyRingElem}
  @assert !iszero(f)
  k = length(G) + 1

  # stats
  l0 = length(Q)

  # product criterion: ignore (i, k) if LM(G[i]) and LM(f) are coprime
  for (i, g) in enumerate(G)
    if !iscoprime(f, g, o)
      push!(Q, (i, k))
    end
  end

  # stats
  created = length(Q) - l0
  inc!(STATS, :critical_pairs_pruned, length(G) - created)
  inc!(STATS, :critical_pairs_created, created)
end

struct GebauerMoellerInstallation <: EliminationStrategy end

function update!(
  ::GebauerMoellerInstallation, Q::Vector{Tuple{Int,Int}}, G::Vector{T}, f::T,
  o::MonomialOrdering,
) where {T<:MPolyRingElem}
  @assert !iszero(f)
  k = lastindex(G) + 1
  lf = lm(f, o)

  # TODO: persist these data structures in the strategy struct
  H = Vector{T}(undef, length(G))  # H[i] = lcm(LM(f), LM(G[i]))
  lcms = Dict{T,Vector{Int}}()
  min_lcms = T[]
  for (i, g) in enumerate(G)
    lfg = lcm(lf, lm(g, o))
    H[i] = lfg
    # TODO: move the product criterion to here, as we don't even need to consider the lcms for the dropped pairs
    push!(get!(Vector{Int}, lcms, lfg), i)
  end

  # stats
  l0 = length(Q)

  # B-criterion: drop (i,j) if LM(f) | lcm(LM(G[i]),LM(G[j])) and lcm(LM(f),LM(G[i])) != lcm(LM(G[i]),LM(G[j])) != lcm(LM(f),LM(G[j]))
  filter!(Q) do pair
    (i, j) = pair
    # note: if this becomes performance critical, allocations here can be easily removed by basic performing operations on exponent vectors
    l = lcm(lm(G[i], o), lm(G[j], o))
    return !reduces(lf, l, o) || H[i] == l || H[j] == l
  end

  # stats
  inc!(STATS, :critical_pairs_discarded, l0 - length(Q))
  l0 = length(Q)

  # M-criterion: ignore (i, k) if H[i] is not minimal in H w.r.t. divisibility
  # F-criterion: ignore (i, k) if there exists j < i with H[j] = H[i]
  for l in sort(collect(keys(lcms)); by=f -> degree(f, o))
    if !any(reduces(m, l, o) for m in min_lcms)
      # l is minimal in H w.r.t. divisibility
      push!(min_lcms, l)

      # product criterion: ignore (i, k) if LM(G[i]) and LM(f) are coprime
      if !any(iscoprime(G[i], f, o) for i in lcms[l])
        # i is minimal with H[i] = l
        i = first(lcms[l])

        # TODO: modify here if any order is required on Q, e.g. if Q should be sorted
        push!(Q, (i, k))
      end
    end
  end

  # stats
  created = length(Q) - l0
  inc!(STATS, :critical_pairs_pruned, length(G) - created)
  inc!(STATS, :critical_pairs_created, created)
end
