using BuchbergerLearning
using Test
using Aqua
using JET

@testset "BuchbergerLearning.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(BuchbergerLearning)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(BuchbergerLearning; target_defined_modules = true)
    end
    # Write your tests here.
end
