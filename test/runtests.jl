using Test
using ControlSystemsBase

using ControlTimingSafety

@testset "ControlTimingSafety.jl" begin
    @testset "Automaton constructors" begin
        sysd = ss([0.5 0.1; 0.02 0.9], [0.5; 0.1], [1 -1], [0], 1)
        K = [0.3 0.2]

        @testset "Hold&Skip-Next" begin
            sysd_hs = hold_skip_next(sysd, K)
            @test sysd_hs.Φ[1] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 0 0 0 0 0; 0 0 0 0 0; -0.3 -0.2 0 0 0]
            @test sysd_hs.Φ[2] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 0 0 0 0 0; 0 0 0 0 0; 0 0 -0.3 -0.2 0]
            @test sysd_hs.Φ[3] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 1]
            @test sysd_hs.Φ[4] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]
        end
        @testset "Zero&Skip-Next" begin
            sysd_zs = zero_skip_next(sysd, K)
            @test sysd_zs.Φ[1] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 0 0 0 0 0; 0 0 0 0 0; -0.3 -0.2 0 0 0]
            @test sysd_zs.Φ[2] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 0 0 0 0 0; 0 0 0 0 0; 0 0 -0.3 -0.2 0]
            @test sysd_zs.Φ[3] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0]
            @test sysd_zs.Φ[4] == [0.5 0.1 0 0 0.5; 0.02 0.9 0 0 0.1; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0]
        end
        @testset "Hold&Kill" begin
            sysd_hk = hold_kill(sysd, K)
            @test sysd_hk.Φ[1] == [0.5 0.1 0.5; 0.02 0.9 0.1; -0.3 -0.2 0]
            @test sysd_hk.Φ[2] == [0.5 0.1 0.5; 0.02 0.9 0.1; 0 0 1]
        end
        @testset "Zero&Kill" begin
            sysd_zk = zero_kill(sysd, K)
            @test sysd_zk.Φ[1] == [0.5 0.1 0.5; 0.02 0.9 0.1; -0.3 -0.2 0]
            @test sysd_zk.Φ[2] == [0.5 0.1 0.5; 0.02 0.9 0.1; 0 0 0]
        end
    end

    include("schedule_synthesis_tests.jl")
end
