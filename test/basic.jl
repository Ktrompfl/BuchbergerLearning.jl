using BuchbergerLearning, Oscar, Test

@testset "S-polynomials" begin
  R, (x, y) = polynomial_ring(GF(32003), [:x, :y]; internal_ordering=:degrevlex)
  o = degrevlex(R)

  @test BuchbergerLearning.spoly(x^2, y^2, o) == R(0)
  @test BuchbergerLearning.spoly(x * y + x, y^2 + y, o) == R(0)
  @test BuchbergerLearning.spoly(x^2 + 1, y^2 + 1, o) == y^2 - x^2
  @test BuchbergerLearning.spoly(y^2 + 1, x^2 + 1, o) == x^2 - y^2
  @test BuchbergerLearning.spoly(x^2 + y, x * y + x, o) == y^2 - x^2
  @test BuchbergerLearning.spoly(x^2 + y, y^2 + y, o) == y^3 - x^2 * y
  @test BuchbergerLearning.spoly(2x^2 - 3y, 5x * y - 4, o) == 12802x + 16000y^2
end

@testset "reduction" begin
  R, (x, y) = polynomial_ring(GF(32003), [:x, :y]; internal_ordering=:degrevlex)
  o = degrevlex(R)

  function trie(G)
    root = BuchbergerLearning.TrieNode()
    for (i, g) in enumerate(G)
      BuchbergerLearning.add_reducer!(root, g, i, o)
    end
    return root
  end

  @test BuchbergerLearning.reduce!(x^2, [x], trie([x]), o; full=true) == R(0)
  @test BuchbergerLearning.reduce!(x^2, [y], trie([y]), o; full=true) == x^2
  @test BuchbergerLearning.reduce!(x^2 + y, [y], trie([y]), o; full=true) == x^2
  @test BuchbergerLearning.reduce!(x^2 + y, [y], trie([y]), o; full=false) == x^2 + y
  @test BuchbergerLearning.reduce!(x^2 + y, [x, y], trie([x, y]), o; full=true) == R(0)
  @test BuchbergerLearning.reduce!(x^2 + y, [x, y], trie([x, y]), o; full=false) == R(0)
  @test BuchbergerLearning.reduce!(
    x * y, [x - 1, y - 1], trie([x - 1, y - 1]), o; full=true
  ) ==
    R(1)
  @test BuchbergerLearning.reduce!(
    x * y, [y - 1, x - 1], trie([y - 1, x - 1]), o; full=true
  ) ==
    R(1)
end

@testset "normalization" begin
  R, (x, y) = polynomial_ring(GF(32003), [:x, :y]; internal_ordering=:degrevlex)
  o = degrevlex(R)

  @test BuchbergerLearning.normalize(x^2 + 2y, o) == x^2 + 2y
  @test BuchbergerLearning.normalize(2x^2 + y, o) == x^2 + 16002y
end

@testset "interreduction" begin
  # TODO
end

@testset "basic" begin
  gb(I) = Oscar.groebner_basis(I; algorithm=:buchberger, complete_reduction=true)

  R, (x, y, z) = polynomial_ring(GF(32003), [:x, :y, :z]; internal_ordering=:degrevlex)

  I = ideal(R, [y - x^2, z - x^3])

  G = buchberger(I; elimination_strategy=NoStrategy(), selection_strategy=DegreeStrategy())
  @test elements(G) == elements(gb(I)) broken = true  # the output is not reduced right now
end

@testset "elimination" begin end
