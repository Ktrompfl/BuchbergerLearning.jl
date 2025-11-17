module BuchbergerLearning

export buchberger

using Oscar
using Random  # note: for random ideals/polynomials from various distributions look into Distribution.jl
using Base.ScopedValues

include("rpoly.jl")
include("polynomials.jl")
include("elimination.jl")
include("reduction.jl")
include("selection.jl")

const CriticalPair = Tuple{Int,Int}

Base.@kwdef mutable struct Stats
  # reductions
  reductions_total::Int = 0
  reductions_nonzero::Int = 0

  # critical pair selection
  critical_pairs_pruned::Int = 0
  critical_pairs_created::Int = 0
  critical_pairs_selected::Int = 0
  critical_pairs_discarded::Int = 0

  s_polynomials::Int = 0

  # generators
  generator_max_degree_pre_reduction::Int = 0
  generator_max_degree_post_reduction::Int = 0
  generator_count_pre_reduction::Int = 0
  generator_count_post_reduction::Int = 0

  # operations
  polynomial_additions::Int = 0
  monomial_additions::Int = 0
end

const STATS = ScopedValue(Stats())

function inc!(ctx, s::Symbol, val::Int=1)
  stats = ctx[]
  if hasfield(Stats, s)
    setfield!(stats, s, getfield(stats, s) + val)
  end
end

@doc raw"""
    spoly(f1::T, f2::T, o::MonomialOrdering) where {T<:MPolyRingElem}

Return the S-polynomial for $f1$ and $f2$.
"""
function spoly(f1::T, f2::T, o::MonomialOrdering) where {T<:MPolyRingElem}
  @assert !iszero(f1)
  @assert !iszero(f2)

  inc!(STATS, :s_polynomials)

  C = MPolyBuildCtx(parent(f1))
  e1 = lx(f1, o)
  e2 = lx(f2, o)
  # TODO: computing the inverses can be avoided by by instead computing LC(fj) ⋅ lcm(LM(fi), LM(fj)) / LM(fi)
  l1 = finish(push_term!(C, inv(lc(f1, o)), max.(e1, e2) .- e1))  # lcm(LM(f1), LM(f2)) / LT(f1)
  l2 = finish(push_term!(C, inv(lc(f2, o)), max.(e1, e2) .- e2))  # lcm(LM(f1), LM(f2)) / LT(f2)

  # s = l1 * f1 - l2 * f2
  l1 = mul!(l1, f1)
  s = submul!(l1, l2, f2)
  inc!(STATS, :additions)

  return s
end

@doc raw"""
  add_generator!(G::Vector{T}, R::TrieNode, f::T, o::MonomialOrdering) where {T<:MPolyRingElem}

Add $f$ to the generators $G$ and reducers $R$.
"""
function add_generator!(
  G::Vector{T}, R::TrieNode, f::T, o::MonomialOrdering
) where {T<:MPolyRingElem}
  @assert !iszero(f)

  # add f to generators G
  inc!(STATS, :generators)
  push!(G, f)

  # add f to reducers R
  k = length(G)
  add_reducer!(R, f, k, o)
end

@doc raw"""
    buchberger(I::MPolyIdeal{T}; complete_reduction::Bool=true, elimination_strategy::EliminationStrategy, selection_strategy::SelectionStrategy) where {T<:MPolyRingElem}

If the internal ordering of $I$ is global, return a Gröbner basis of $I$ with respect to this ordering.

!!! note
    The returned Gröbner basis is reduced if `complete_reduction = true`.
"""
function buchberger(
  I::MPolyIdeal{T};
  complete_reduction::Bool=true,
  o::MonomialOrdering=default_ordering(base_ring(I)),
  elimination_strategy::EliminationStrategy,
  selection_strategy::SelectionStrategy,
) where {T<:MPolyRingElem}
  is_global(o) || error("Ordering must be global")

  # box all generators (cache leading indices w.r.t. monomial ordering)
  S = reduction_polynomial_ring(base_ring(I), o)
  I = map(S, gens(I))

  # TODO: glaze all polynomials in I for sugar / spice selection strategies

  G = Vector{elem_type(S)}()  # generators
  sizehint!(G, length(I))

  R = init_reducers()  # reducers
  Q = CriticalPair[]  # critical pairs

  stats = Stats()

  # use scoped value so we don't have to pass down stats the entire call chain
  # warning: this is not thread safe, every thread needs its own scope
  @with STATS => stats begin
    for g in I
      iszero(g) && continue

      # update the critical pairs
      update!(elimination_strategy, Q, G, g, o)

      # add g to G and R
      add_generator!(G, R, g, o)
    end

    while !isempty(Q)
      (i, j) = popat!(Q, select(selection_strategy, G, Q, o))
      s = spoly(G[i], G[j], o)
      g = reduce!(s, G, R, o; full=true)

      iszero(g) && continue

      inc!(STATS, :reductions_nonzero)

      # update the critical pairs
      update!(elimination_strategy, Q, G, g, o)

      # add g to G and R
      add_generator!(G, R, g, o)
    end

    stats.generator_count_pre_reduction = length(G)
    stats.generator_max_degree_pre_reduction = maximum(degree(g, o) for g in G)

    if complete_reduction
      # transform G to the reduced Gröbner basis of ⟨G⟩
      error("interreduction not implemented for current reduction process")
      # interreduce!(G)
    end

    stats.generator_count_post_reduction = length(G)
    stats.generator_max_degree_post_reduction = maximum(degree(g, o) for g in G)
  end

  # unbox generators
  S = handle(S)
  G = map(handle, G)

  # build output
  J = Oscar.IdealGens(S, G, o)
  J.isGB = true
  J.isReduced = complete_reduction

  return J
end

end # module BuchbergerLearning
