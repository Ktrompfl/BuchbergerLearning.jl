
using Oscar

import Base: show, +, -, *, ^, ==, inv, isone, iszero, isequal, one, zero,
  deepcopy_internal, hash, length

# TODO: how to handle this correctly?
import Oscar.AbstractAlgebra: MPolyBuildCtx, push_term!, finish, parent_type, elem_type,
  base_ring_type,
  base_ring, parent, is_domain_type, is_exact_type, symbols, number_of_variables, gens, gen,
  internal_ordering, degrees, total_degree, is_gen, coefficients, exponent_vectors,
  characteristic, is_unit, canonical_unit, divides, divexact, divrem, div, gcd, lcm, remove,
  valuation, zero!, one!, mul!, add!, sub!, neg!, inv!, addmul!, submul!, coeff, monomial,
  term, exponent, exponent_vector

import Oscar: index_of_leading_term, default_ordering

export reduction_polynomial_ring

mutable struct ReductionMPolyRing{T<:RingElement,R<:MPolyRing{T}} <: MPolyRing{T}
  ring::R
  ordering::MonomialOrdering{R}

  function ReductionMPolyRing(
    ring::R, ordering::MonomialOrdering{R}, cached::Bool=false
  ) where {T<:RingElem,R<:MPolyRing{T}}
    # get_cached! attaches a finalizer, which requires the struct to be mutable
    return AbstractAlgebra.get_cached!(
      ReductionMPolyID, (ring, ordering), cached
    ) do
      new{T,R}(
        ring, ordering
      )
    end::ReductionMPolyRing{T,R}
  end
end

const ReductionMPolyID = AbstractAlgebra.CacheDictType{
  Tuple{MPolyRing,MonomialOrdering},ReductionMPolyRing
}()

mutable struct ReductionMPoly{T<:RingElement,P<:MPolyRingElem{T}} <: MPolyRingElem{T}
  parent::ReductionMPolyRing{T}
  poly::P
  lead_index::Int
end

function ReductionMPoly(
  ring::ReductionMPolyRing{T}, poly::P
) where {T<:RingElement,P<:MPolyRingElem{T}}
  idx = iszero(poly) ? 0 : index_of_leading_term(poly, ring.ordering)  # lead index of 0 is undefined
  return ReductionMPoly(ring, poly, idx)
end

# build from terms

struct ReductionMPolyBuildCtx
  ring::ReductionMPolyRing
  builder::MPolyBuildCtx
end

MPolyBuildCtx(R::ReductionMPolyRing) = ReductionMPolyBuildCtx(R, MPolyBuildCtx(handle(R)))

function push_term!(m::ReductionMPolyBuildCtx, c::RingElem, v::Vector{Int})
  push_term!(m.builder, c, v)
  return m
end

function finish(m::ReductionMPolyBuildCtx)
  f = finish(m.builder)
  return ReductionMPoly(m.ring, f)
end

# Data type and parent object methods

ordering(R::ReductionMPolyRing) = R.ordering
handle(R::ReductionMPolyRing) = R.ring
handle(f::ReductionMPoly) = f.poly

parent_type(::Type{ReductionMPoly{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  ReductionMPolyRing{T,parent_type(P)}

elem_type(::Type{ReductionMPolyRing{T,P}}) where {T<:RingElement,P<:MPolyRing{T}} =
  ReductionMPoly{T,elem_type(P)}

base_ring_type(::Type{ReductionMPolyRing{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  base_ring_type(P)

base_ring(R::ReductionMPolyRing) = base_ring(R.ring)

parent(f::ReductionMPoly) = f.parent

is_domain_type(::Type{ReductionMPoly{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  is_domain_type(P)

is_exact_type(::Type{ReductionMPoly{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  is_exact_type(P)

hash(f::ReductionMPoly, h::UInt) = hash(handle(f), h)

deepcopy_internal(f::ReductionMPoly, dict::IdDict) =
  ReductionMPoly(parent(f), deepcopy_internal(handle(f), dict), f.lead_index)

symbols(R::ReductionMPolyRing) = symbols(R.ring)

number_of_variables(R::ReductionMPolyRing) = number_of_variables(handle(R))

gens(R::ReductionMPolyRing) = map(R, gens(handle(R)))

gen(R::ReductionMPolyRing, i::Int) = R(gen(handle(R), i))

internal_ordering(R::ReductionMPolyRing) = internal_ordering(handle(R))

length(f::ReductionMPoly) = length(handle(f))

degrees(f::ReductionMPoly) = degrees(handle(f))

total_degree(f::ReductionMPoly) = total_degree(handle(f))

is_gen(f::ReductionMPoly) = is_gen(handle(f))

coefficients(f::ReductionMPoly) = coefficients(handle(f))

exponent_vectors(f::ReductionMPoly) = exponent_vectors(handle(f))

# TODO: monomials and terms would need a custom iterator type, but a generic version is provided

# Basic manipulation

characteristic(R::ReductionMPolyRing) = characteristic(handle(R))

zero(R::ReductionMPolyRing) = ReductionMPoly(R, zero(handle(R)), 0)

one(R::ReductionMPolyRing) = ReductionMPoly(R, one(handle(R)), 1)

iszero(f::ReductionMPoly) = iszero(handle(f))

isone(f::ReductionMPoly) = isone(handle(f))

is_unit(f::ReductionMPoly) = is_unit(handle(f))

canonical_unit(f::ReductionMPoly) = ReductionMPoly(parent(f), canonical_unit(handle(f)), 1)

# String I/O

function show(io::IO, f::ReductionMPoly)
  show(io, handle(f))
end

# safe operations

function -(f::ReductionMPoly)
  return ReductionMPoly(parent(f), -handle(f))
end

function +(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(parent(f), handle(f) + handle(g))
end

function -(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(parent(f), handle(f) - handle(g))
end

function *(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(parent(f), handle(f) * handle(g))
end

function *(f::ReductionMPoly{T}, c::T) where {T<:RingElem}
  return ReductionMPoly(parent(f), handle(f) * c)
end

function *(c::T, f::ReductionMPoly{T}) where {T<:RingElem}
  return ReductionMPoly(parent(f), c * handle(f))
end

function ^(f::ReductionMPoly, e::Int)
  return ReductionMPoly(parent(f), handle(f)^e)
end

# comparison

function ==(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return handle(f) == handle(g)
end

function isequal(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return isequal(handle(f), handle(g))
end

# division

function divides(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return divides(handle(f), handle(g))
end

function divexact(f::ReductionMPoly, g::ReductionMPoly; check::Bool=true)
  check_parent(f, g)
  return ReductionMPoly(parent(f), divexact(handle(f), handle(g); check=check))
end

function divrem(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  (q, r) = divrem(handle(f), handle(g))
  return (ReductionMPoly(parent(f), q), ReductionMPoly(parent(f), r))
end

function div(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(parent(f), div(handle(f), handle(g)))
end

function gcd(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(parent(f), gcd(handle(f), handle(g)))
end

function lcm(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(parent(f), lcm(handle(f), handle(g)))
end

function remove(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  (v, q) = remove(handle(f), handle(g))
  return (v, ReductionMPoly(parent(f), q))
end

function valuation(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return valuation(handle(f), handle(g))
end

function divexact(f::ReductionMPoly, c::Integer)
  return ReductionMPoly(parent(f), divexact(handle(f), c))
end

function divexact(f::ReductionMPoly, c::Rational)
  return ReductionMPoly(parent(f), divexact(handle(f), c))
end

function divexact(f::ReductionMPoly{T}, c::T) where {T<:RingElem}
  return ReductionMPoly(parent(f), divexact(handle(f), c))
end

# inverse

function inv(f::ReductionMPoly)
  return ReductionMPoly(parent(f), inv(handle(f)))
end

# unsafe operations

function zero!(f::ReductionMPoly)
  f.poly = zero!(handle(f))
  # lead index of 0 is undefined
  return f
end

function one!(f::ReductionMPoly)
  f.poly = one!(handle(f))
  f.lead_index = 1
  return f
end

function _recompute_lead_index!(f::ReductionMPoly)
  iszero(f) && return nothing  # lead index of 0 is undefined
  f.lead_index = index_of_leading_term(handle(f), parent(f).ordering)
end

function mul!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly)
  f.poly = mul!(handle(f), handle(g), handle(h))
  _recompute_lead_index!(f)
  return f
end

function add!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly)
  f.poly = add!(handle(f), handle(g), handle(h))
  _recompute_lead_index!(f)
  return f
end

function sub!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly)
  f.poly = sub!(handle(f), handle(g), handle(h))
  _recompute_lead_index!(f)
  return f
end

function neg!(f::ReductionMPoly, g::ReductionMPoly)
  f.poly = neg!(handle(f), handle(g))
  # lead index unchanged
  return f
end

function inv!(f::ReductionMPoly, g::ReductionMPoly)
  f.poly = inv!(handle(f), handle(g))
  _recompute_lead_index!(f)
  return f
end

function addmul!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly, t::ReductionMPoly)
  f.poly = addmul!(handle(f), handle(g), handle(h), handle(t))
  _recompute_lead_index!(f)
  return f
end

function submul!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly, t::ReductionMPoly)
  f.poly = submul!(handle(f), handle(g), handle(h), handle(t))
  _recompute_lead_index!(f)
  return f
end

# accessors

function index_of_leading_term(f::ReductionMPoly, o::MonomialOrdering)
  iszero(f) && error("index of leading term undefined for zero")
  return o == ordering(parent(f)) ? f.lead_index : index_of_leading_term(handle(f), o)
end

function coeff(f::ReductionMPoly, n::Int)
  return coeff(handle(f), n)
end

function coeff(f::ReductionMPoly, exps::Vector{Int})
  return coeff(handle(f), exps)
end

function monomial(f::ReductionMPoly, n::Int)
  return ReductionMPoly(parent(f), monomial(handle(f), n), 1)
end

function term(f::ReductionMPoly, n::Int)
  return ReductionMPoly(parent(f), term(handle(f), n), 1)
end

function exponent(f::ReductionMPoly, i::Int, j::Int)
  return exponent(handle(f), i, j)
end

function exponent_vector(f::ReductionMPoly, i::Int)
  return exponent_vector(handle(f), i)
end

# function setcoeff!(a::MyMPoly, exps::Vector{Int}, c::S) where {S<:RingElement}
# end

# function fit!(f::MyMPoly{T}, n::Int) where {T<:RingElem}
# end

# function setcoeff!(a::MyMPoly{T}, i::Int, c::T) where {T<:RingElement}
# end

# function setcoeff!(a::MyMPoly{T}, i::Int, c::U) where {T<:RingElement,U<:Integer}
# end

# function combine_like_terms!(a::MyMPoly{T}) where {T<:RingElement}
# end

# function set_exponent_vector!(
#   a::MyMPoly{T}, i::Int, exps::Vector{Int}
# ) where {T<:RingElement}
# end

# function sort_terms!(a::MyMPoly{T}) where {T<:RingElement}
# end

# Random generation

# RandomExtensions.maketype(R::ConstPolyRing, _) = elem_type(R)

# rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{ConstPoly,ConstPolyRing}}) =
#   sp[][1](rand(rng, sp[][2]))

# rand(rng::AbstractRNG, R::ConstPolyRing, n::AbstractUnitRange{Int}) = R(rand(rng, n))

# rand(R::ConstPolyRing, n::AbstractUnitRange{Int}) = rand(Random.default_rng(), R, n)

# promotion rules

Base.promote_rule(
  ::Type{ReductionMPoly{T,P}}, ::Type{ReductionMPoly{T,P}}
) where {T<:RingElement,P<:MPolyRingElem{T}} = ReductionMPoly{T,P}

# polynomial constructors

function (R::ReductionMPolyRing)()
  return zero(R)
end

function (R::ReductionMPolyRing)(c::Integer)
  return ReductionMPoly(R, handle(R)(c), 1)
end

function (R::ReductionMPolyRing{T})(c::T) where {T<:RingElement}
  return ReductionMPoly(R, handle(R)(c), 1)
end

function (R::ReductionMPolyRing)(f::MPolyRingElem)
  handle(R) != parent(f) && error("Unable to coerce element")
  return ReductionMPoly(R, f)
end

function (R::ReductionMPolyRing{T})(f::ReductionMPoly{T}) where {T<:RingElement}
  R != parent(f) && error("Unable to coerce element")
  return f
end

# ring constructor

function reduction_polynomial_ring(
  R::T, ordering::MonomialOrdering{T}=default_ordering(R), cached::Bool=true
) where {T<:MPolyRing}
  return ReductionMPolyRing(R, ordering, cached)
end
