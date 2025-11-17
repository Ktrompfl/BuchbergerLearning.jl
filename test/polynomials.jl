using BuchbergerLearning, Oscar, Test

@testset "polynomials" begin
  K = GF(32003)
  R, (x, y) = polynomial_ring(K, [:x, :y]; internal_ordering=:lex)
  S, (a, b) = polynomial_ring(K, [:a, :b]; internal_ordering=:degrevlex)

  o1 = degrevlex(R)
  o2 = lex(S)

  @testset "lc" begin
    @test BuchbergerLearning.lc(x^2 + 2y + 3, o1) == K(1)
    @test BuchbergerLearning.lc(2y + x^2 + 3, o1) == K(1)
    @test BuchbergerLearning.lc(2y^3 + 1, o1) == K(2)
    @test BuchbergerLearning.lc(-y^3, o1) == K(-1)
  end

  @testset "lm" begin
    @test BuchbergerLearning.lm(x^2 + 2y + 3, o1) == x^2
    @test BuchbergerLearning.lm(2y + x^2 + 3, o1) == x^2
    @test BuchbergerLearning.lm(2y^3 + 1, o1) == y^3
    @test BuchbergerLearning.lm(x^3 + 5y^3, o1) == x^3
  end

  @testset "lt" begin
    @test BuchbergerLearning.lt(x^2 + 2y + 3, o1) == x^2
    @test BuchbergerLearning.lt(2y + x^2 + 3, o1) == x^2
    @test BuchbergerLearning.lt(2y^3 + 1, o1) == 2y^3
    @test BuchbergerLearning.lt(x^3 + 5y^3, o1) == x^3
  end

  @testset "lx" begin
    @test BuchbergerLearning.lx(x^2 + 2y + 3, o1) == [2, 0]
    @test BuchbergerLearning.lx(2y + x^2 + 3, o1) == [2, 0]
    @test BuchbergerLearning.lx(2x * y^3 + 1, o1) == [1, 3]
    @test BuchbergerLearning.lx(x^3 + 5y^3, o1) == [3, 0]
  end

  @testset "degree" begin
    # degrevlex
    @test BuchbergerLearning.degree(x^5, o1) == 5
    @test BuchbergerLearning.degree(x^5 + y^4, o1) == 5
    @test BuchbergerLearning.degree(y^5 + x^4, o1) == 5

    # lex
    @test BuchbergerLearning.degree(a^5, o2) == 5
    @test BuchbergerLearning.degree(a^5 + b^4, o2) == 5
    @test BuchbergerLearning.degree(b^5 + a^4, o2) == 4
  end

  @testset "iscoprime" begin
    @test BuchbergerLearning.iscoprime(x^2, y^2, o1) == true
    @test BuchbergerLearning.iscoprime(x^2, x * y, o1) == false
    @test BuchbergerLearning.iscoprime(x^2 + y, y^2, o1) == true
    @test BuchbergerLearning.iscoprime(2x, 3x, o1) == false
  end

  @testset "ismonic" begin
    @test BuchbergerLearning.ismonic(x^4, o1) == true
    @test BuchbergerLearning.ismonic(2x^4, o1) == false
    @test BuchbergerLearning.ismonic((-x) * (-y), o1) == true
    @test BuchbergerLearning.ismonic(x^3 + 3y + 2, o1) == true
    @test BuchbergerLearning.ismonic(4x^2 + y + 1, o1) == false
  end

  @testset "reduces" begin
    @test BuchbergerLearning.reduces(x^2, x^4, o1) == true
    @test BuchbergerLearning.reduces(x^4, x^2, o1) == false
    @test BuchbergerLearning.reduces(x^3, x^3, o1) == true
    @test BuchbergerLearning.reduces(x * y, x^4, o1) == false
    @test BuchbergerLearning.reduces(x^2 + y, x^3, o1) == true
    @test BuchbergerLearning.reduces(y^2, x^3 + y^2, o1) == false
  end
end
