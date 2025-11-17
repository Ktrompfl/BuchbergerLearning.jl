import Base: iterate

struct TrieNode
  polys::Vector{Int}
  children::Dict{Int,TrieNode}

  TrieNode() = new(Int[], Dict{Int,TrieNode}())
end

@doc raw"""
  add_reducer!(node::TrieNode, f::MPolyRingElem, idx::Int, o::MonomialOrdering)

Add the reducer $f$ at index $idx$ to the trie.
"""
function add_reducer!(node::TrieNode, f::MPolyRingElem, idx::Int, o::MonomialOrdering)
  @assert idx > 0
  n = ngens(parent(f))
  l = li(f, o)

  for x in 1:n
    node = get!(TrieNode, node.children, exponent(f, l, x))
  end

  push!(node.polys, idx)
end

@doc raw"""
  add_reducer!(node::TrieNode, f::MPolyRingElem, idx::Int, o::MonomialOrdering)

Remove the reducer $f$ at index $idx$ from the trie, if it exists.
"""
function remove_reducer!(node::TrieNode, f::MPolyRingElem, idx::Int, o::MonomialOrdering)
  @assert idx > 0
  n = ngens(parent(f))
  l = li(f, o)

  for x in 1:n
    node = get(node.children, exponent(f, l, x), nothing)
    isnothing(node) && return nothing
  end

  deleteat!(node.polys, findfirst(==(idx), node.polys))
end

@doc raw"""
  find_reducers(trie::TrieNode, f::MPolyRingElem, o::MonomialOrdering)

Return an iterator over the indices of all polynomials in the trie reducing $f$.
"""
function find_reducers(trie::TrieNode, f::MPolyRingElem, o::MonomialOrdering)
  l = li(f, o)
  return LexicographicalDivisorIterator(trie, f, l)
end

function init_reducers()
  return TrieNode()
end

struct LexicographicalDivisorIterator{T<:MPolyRingElem}
  root::TrieNode
  poly::T
  index::Int  # index of leading term
end

# state: vector of (node, x, child_iter_state)
function iterate(it::LexicographicalDivisorIterator)
  stack = Tuple{TrieNode,Int,Union{Nothing,Some{Any}}}[(it.root, 0, nothing)]
  return iterate(it, stack)
end

function iterate(it::LexicographicalDivisorIterator, stack)
  n = ngens(parent(it.poly))

  while !isempty(stack)
    node, depth, child_iter_state = pop!(stack)

    if depth == n
      @assert !isempty(node.polys)  # no empty leaves

      nxt = if child_iter_state === nothing
        iterate(node.polys)
      else
        iterate(node.polys, something(child_iter_state))
      end

      nxt === nothing && continue

      # no backtrack
      poly, state = nxt

      # update current frame iterator
      push!(stack, (node, depth, Some(state)))

      return poly, stack
    end

    nxt = if child_iter_state === nothing
      iterate(node.children)
    else
      iterate(node.children, something(child_iter_state))
    end

    while nxt !== nothing
      (exp, child), state = nxt
      if exp <= exponent(it.poly, it.index, depth + 1)

        # no backtrack
        push!(stack, (node, depth, Some(state)))

        # descend
        push!(stack, (child, depth + 1, nothing))
        break
      end

      nxt = iterate(node.children, state)
    end
  end

  return nothing
end

function reduce!(
  f::T, G::Vector{T}, R::TrieNode, o::MonomialOrdering; full::Bool=true
) where {T<:MPolyRingElem}
  inc!(STATS, :reductions)
  r = MPolyBuildCtx(parent(f))

  while !iszero(f)
    inc!(STATS, :reductions_steps)
    # reduce lead term of f with the first found reducer g ∈ R
    # TODO: add more sophisticated strategies
    reducer = iterate(find_reducers(R, f, o))  # right now just take first reducer, but this can be subject to strategies
    if isnothing(reducer)
      # lead term not reducible
      full || return f
      # inc!(STATS, :additions, 2)  - these are only additions because of poor implementation, so they shouldn't count
      push_term!(r, lc(f, o), lx(f, o))  # r += LT(f)
      f = sub!(f, lt(f, o)) # f = tail(f)
    else
      i, _ = reducer
      g = G[i]

      # reduce lead term
      inc!(STATS, :additions)
      f = submul!(f, divexact(lt(f, o), lt(g, o)), g)  # TODO: can provide tmp polynomial for efficiency
    end
  end

  return finish(r)
end
