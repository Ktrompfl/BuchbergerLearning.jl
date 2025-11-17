
using AbstractAlgebra

using Random: Random, SamplerTrivial
using AbstractAlgebra.RandomExtensions: RandomExtensions, Make2, AbstractRNG

import AbstractAlgebra: parent_type, elem_type, base_ring, base_ring_type, parent,
  is_domain_type,
  is_exact_type, canonical_unit, isequal, divexact, zero!, mul!, add!,
  get_cached!, is_unit, characteristic, Ring, RingElem

import Oscar: MonomialOrdering

import Base:
  show, +, -, *, ^, ==, inv, isone, iszero, one, zero, rand,
  deepcopy_internal, hash

const ReductionMPolyID = AbstractAlgebra.CacheDictType{MPolyRing,ReductionMPolyRing}()

struct ReductionMPolyRing{T<:RingElement,R<:MPolyRing{T}} <: MPolyRing{T}
  ring::R
  ordering::MonomialOrdering

  function ReductionMPolyRing(R::MPolyRing, cached::Bool=false)
    return get_cached!(
      ReductionMPolyID, R, cached
    ) do
      new{elem_type(R),typeof(R)}(R)
    end::ReductionMPolyRing{elem_type(R),typeof(R)}
  end
end

mutable struct ReductionMPoly{T<:RingElement,P::MPolyRingElem{T}} <:
               MPolyRingElem{T}
  const parent::ReductionMPolyRing{T}
  poly::P
  lead_index::Int
end

function ReductionMPoly(
  ring::ReductionMPolyRing{T}, poly::P
) where {T<:RingElement,P::MPolyRingElem{T}}
  idx = index_of_leading_term(poly, ring.ordering)
  return ReductionMPoly(ring, poly, idx)
end

# build from terms

mutable struct ReductionMPolyBuildCtx
  ring::ReductionMPolyRing
  builder::MPolyBuildCtx

  ReductionMPolyBuildCtx(R::ReductionMPolyRing) = new(R, MPolyBuildCtx(R.ring))
end

function push_term!(m::ReductionMPolyBuildCtx, c::RingElem, v::Vector{Int})
  push_term!(m.builder, c, v)
end

function finish(m::ReductionMPolyBuildCtx)
  f = finish(m.builder)
  return ReductionMPoly(m.ring, f)
end

# Data type and parent object methods

parent_type(::Type{ReductionMPoly{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  ReductionMPolyRing{T}  # FIXME

elem_type(::Type{ReductionMPolyRing{T}}) where {T<:RingElement} = ReductionMPoly{T}  # FIXME

base_ring_type(::Type{ConstPolyRing{T}}) where {T<:RingElement} = parent_type(T)  # FIXME

base_ring(R::ReductionMPolyRing) = base_ring(R.ring)

parent(f::ReductionMPoly) = f.parent

is_domain_type(::Type{ReductionMPoly{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  is_domain_type(P)

is_exact_type(::Type{ReductionMPoly{T,P}}) where {T<:RingElement,P<:MPolyRingElem{T}} =
  is_exact_type(P)

hash(f::ReductionMPoly, h::UInt) = hash(f.poly, h)

deepcopy_internal(f::ReductionMPoly, dict::IdDict) =
  ReductionMPoly(f.ring, deepcopy_internal(f.poly, dict), f.lead_index)

symbols(R::ReductionMPolyRing) = symbols(R.ring)

number_of_variables(R::ReductionMPolyRing) = number_of_variables(R.ring)

gens(R::ReductionMPolyRing) = map(R, gens(R.ring))

gen(R::ReductionMPolyRing, i::Int) = R(gen(R.ring, i))

internal_ordering(R::ReductionMPolyRing) = internal_ordering(R.ring)

length(f::ReductionMPoly) = length(f.poly)

degrees(f::ReductionMPoly) = degrees(f.poly)

total_degree(f::ReductionMPoly) = total_degree(f.poly)

is_gen(f::ReductionMPoly) = is_gen(f.poly)

coefficients(f::ReductionMPoly) = coefficients(f.poly)

exponent_vectors(f::ReductionMPoly) = exponent_vectors(f.poly)

# TODO: monomials and terms would need a custom iterator type, but a generic version is provided

# Basic manipulation

characteristic(R::ReductionMPolyRing) = characteristic(R.ring)

zero(R::ReductionMPolyRing) = ReductionMPoly(R, zero(R.ring), 0)

one(R::ReductionMPolyRing) = ReductionMPoly(R, one(R.ring), 1)

iszero(f::ReductionMPoly) = iszero(f.poly)

isone(f::ReductionMPoly) = isone(f.poly)

is_unit(f::ReductionMPoly) = is_unit(f.poly)

canonical_unit(f::ReductionMPoly) = ReductionMPoly(f.parent, canonical_unit(f.poly), 1)

# String I/O

function show(io::IO, f::ReductionMPoly)
  show(io, f.poly)
end

# safe operations

function -(f::ReductionMPoly)
  return ReductionMPoly(f.parent, -f.poly)
end

function +(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(f.parent, f.poly + g.poly)
end

function -(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(f.parent, f.poly - g.poly)
end

function *(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return ReductionMPoly(f.parent, f.poly * g.poly)
end

function ^(f::ReductionMPoly, e::Int)
  return ReductionMPoly(f.parent, f.poly^e)
end

# comparison

function ==(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return f.poly == g.poly
end

function isequal(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return isequal(f.poly, g.poly)
end

# division

function divides(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return divides(f.poly, g.poly)
end

function divexact(f::ReductionMPoly, g::ReductionMPoly; check::Bool=true)
  check_parent(f, g)
  return ReductionMPoly(f.parent, divexact(f.poly, g.poly; check=check))
end

function remove(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  (v, q) = remove(f.poly, q.poly)
  return (v, ReductionMPoly(f.parent, q))
end

function valuation(f::ReductionMPoly, g::ReductionMPoly)
  check_parent(f, g)
  return valuation(f.poly, q.poly)
end

function divexact(f::ReductionMPoly, c::Integer)
  return ReductionMPoly(f.parent, divexact(f, c))
end

function divexact(f::ReductionMPoly, c::Rational)
  return ReductionMPoly(f.parent, divexact(f, c))
end

function divexact(f::ReductionMPoly{T}, c::T)
  return ReductionMPoly(f.parent, divexact(f, c))
end

# inverse

function inv(f::ReductionMPoly)
  return ReductionMPoly(f.parent, inv(f.poly))
end

# unsafe operations

function zero!(f::ReductionMPoly)
  f.poly = zero!(f.poly)
  # lead index of 0 is undefined
  return f
end

function one!(f::ReductionMPoly)
  f.poly = one!(f.poly)
  f.lead_index = 1
  return f
end

function mul!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly)
  f.poly = mul!(f.poly, g.poly, h.poly)
  f.lead_index = index_of_leading_term(poly, f.parent.ordering)
  return f
end

function add!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly)
  f.poly = add!(f.poly, g.poly, h.poly)
  f.lead_index = index_of_leading_term(poly, f.parent.ordering)
  return f
end

function sub!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly)
  f.poly = sub!(f.poly, g.poly, h.poly)
  f.lead_index = index_of_leading_term(poly, f.parent.ordering)
  return f
end

function neg!(f::ReductionMPoly, g::ReductionMPoly)
  f.poly = neg!(f.poly, g.poly)
  # lead index unchanged
  return f
end

function inv!(f::ReductionMPoly, g::ReductionMPoly)
  f.poly = inv!(f.poly, g.poly)
  f.lead_index = index_of_leading_term(poly, f.parent.ordering)
  return f
end

function addmul!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly, t::ReductionMPoly)
  f.poly = addmul!(f.poly, g.poly, h.poly, t.poly)
  f.lead_index = index_of_leading_term(poly, f.parent.ordering)
  return f
end

function submul!(f::ReductionMPoly, g::ReductionMPoly, h::ReductionMPoly, t::ReductionMPoly)
  f.poly = submul!(f.poly, g.poly, h.poly, t.poly)
  f.lead_index = index_of_leading_term(poly, f.parent.ordering)
  return f
end

function coeff(f::ReductionMPoly, n::Int)
  return coeff(f.poly, n)
end

function coeff(f::ReductionMPoly, exps::Vector{Int})
  return coeff(f.poly, exps)
end

function monomial(f::ReductionMPoly, n::Int)
  return ReductionMPoly(f.parent, monomial(f.poly, n), 1)
end

function term(f::MyMPoly{T}, n::Int) where {T<:RingElem}
  return ReductionMPoly(f.parent, term(f.poly, n), 1)
end

function exponent(f::ReductionMPoly, i::Int, j::Int)
  return exponent(f.poly, i, j)
end

function exponent_vector(f::ReductionMPoly, i::Int)
  return exponent_vector(f.poly, i)
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

# Promotion rules

# promote_rule(::Type{ConstPoly{T}}, ::Type{ConstPoly{T}}) where {T<:RingElement} =
#   ConstPoly{T}

# function promote_rule(::Type{ConstPoly{T}}, ::Type{U}) where {T<:RingElement,U<:RingElement}
#   promote_rule(T, U) == T ? ConstPoly{T} : Union{}
# end

# # Constructors

# function (R::ConstPolyRing{T})() where {T<:RingElement}
#   p = R.ring()
#   return SugarMPoly(R, p, total_degree(p))
# end

# function (R::ConstPolyRing{T})(c::Integer) where {T<:RingElement}
#   p = R.ring(c)
#   return SugarMPoly(R, p, total_degree(p))
# end

# # Needed to prevent ambiguity
# function (R::SugarMPolyRing{T})(c::T) where {T<:Integer}
#   p = R.ring(c)
#   return ReductionMPoly(R, p)
# end

# function (R::ReductionMPolyRing{T})(c::T) where {T<:RingElement}
#   p = R.ring(c)
#   return ReductionMPoly(R, p)
# end

# function (R::ReductionMPolyRing{T})(f::ReductionMPoly{T}) where {T<:RingElement}
#   R != parent(f) && error("Unable to coerce element")
#   return f
# end

# # Parent constructor

# function reduction_polynomial_ring(R::MPolyRing, cached::Bool=true)
#   return ReductionMPolyRing(R, cached)
# end
