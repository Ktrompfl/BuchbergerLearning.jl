
@doc raw"""
    interreduce!(G::Vector{T}) where {T<:MPolyRingElem}

Transform the Gröbner basis $G$ to the reduced Gröbner basis in-place. Return $G$.
"""
function interreduce!(G::Vector{T}, o::MonomialOrdering) where {T<:MPolyRingElem}
  i = 1
  while i <= length(G)
    m = lx(G[i], o)
    rest = @view G[[begin:(i - 1); (i + 1):end]]
    g = reduce!(G[i], rest, true)

    if iszero(g)
      G[i] = G[end]
      pop!(G)
    elseif lx(g, o) != m
      # lead monomial decreased, restart
      G[i] = g
      i = 1
    else
      i += 1
    end
  end

  # normalize generators
  for i in eachindex(G)
    G[i] = normalize(G[i], o)
  end

  # sort generators by leading monomial
  sort!(G; by=lm)

  return G
end
