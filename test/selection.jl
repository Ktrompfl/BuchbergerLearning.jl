using BuchbergerLearning, Oscar, Test

@testset "Selection strategies" begin
  @testset "first" begin
    R, (x, y) = polynomial_ring(GF(2), [:x, :y]; internal_ordering=:degrevlex)
    o = degrevlex(R)
    S = FirstStrategy()
    G = [x^2, y, x * y, x]

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2

    Q = [(1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1

    Q = [(2, 3), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2

    Q = [(2, 3), (1, 4)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1

    Q = [(1, 4), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
  end

  @testset "last" begin
    R, (x, y) = polynomial_ring(GF(2), [:x, :y]; internal_ordering=:degrevlex)
    o = degrevlex(R)
    S = LastStrategy()
    G = [x^2, y, x * y, x]

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1

    Q = [(1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2

    Q = [(2, 3), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1

    Q = [(2, 3), (1, 4)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2

    Q = [(1, 4), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
  end

  @testset "degree" begin
    R, (x, y) = polynomial_ring(GF(2), [:x, :y]; internal_ordering=:degrevlex)
    o = degrevlex(R)
    S = DegreeStrategy()
    G = [x^3, x * y, y]
    H = [y, x * y, x^3]

    Q = [(1, 2), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(2, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(2, 3), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 2
  end

  @testset "codegree" begin
    R, (x, y) = polynomial_ring(GF(2), [:x, :y]; internal_ordering=:degrevlex)
    o = degrevlex(R)
    S = CodegreeStrategy()
    G = [x^3, x * y, y]
    H = [y, x * y, x^3]

    Q = [(1, 2), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(2, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(2, 3), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 1
  end

  @testset "normal" begin
    # degrevlex
    R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]; internal_ordering=:degrevlex)
    o = degrevlex(R)
    S = NormalStrategy()
    G = [x, y, z, x * y]
    H = [x * y, z, y, x]

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 2), (1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 3
    @test BuchbergerLearning.select(S, H, Q, o) == 3

    Q = [(1, 2), (1, 4), (2, 4)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 3

    # lex
    R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]; internal_ordering=:lex)
    o = lex(R)
    S = NormalStrategy()
    G = [x, y, z, x * y]
    H = [x * y, z, y, x]

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 2), (1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 3
    @test BuchbergerLearning.select(S, H, Q, o) == 3

    Q = [(1, 2), (1, 4), (2, 4)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 3
  end

  @testset "strange" begin
    # degrevlex
    R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]; internal_ordering=:degrevlex)
    o = degrevlex(R)
    S = StrangeStrategy()
    G = [x, y, z, x * y]
    H = [x * y, z, y, x]

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 2), (1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 2), (1, 4), (2, 4)]
    @test BuchbergerLearning.select(S, G, Q, o) == 3
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    # lex
    R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]; internal_ordering=:lex)
    o = lex(R)
    S = StrangeStrategy()
    G = [x, y, z, x * y]
    H = [x * y, z, y, x]

    Q = [(1, 2), (1, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 3), (1, 2)]
    @test BuchbergerLearning.select(S, G, Q, o) == 2
    @test BuchbergerLearning.select(S, H, Q, o) == 2

    Q = [(1, 2), (1, 3), (2, 3)]
    @test BuchbergerLearning.select(S, G, Q, o) == 1
    @test BuchbergerLearning.select(S, H, Q, o) == 1

    Q = [(1, 2), (1, 4), (2, 4)]
    @test BuchbergerLearning.select(S, G, Q, o) == 3
    @test BuchbergerLearning.select(S, H, Q, o) == 1
  end
end
