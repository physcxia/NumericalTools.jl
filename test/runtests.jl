using NumericalTools
using Test

@testset "geomspace" begin
    @test_throws ArgumentError geomspace(1, 2, 1)
    @test_throws ArgumentError geomspace(1.0, 2.0, -2)
    @test length(geomspace(1, 2)) == 50
    v = geomspace(1e-20, 1e20, 41);
    for i = 1:length(v)
        @test v[i] ≈ 1e-20 * 10.0^(i - 1)
    end
    v = geomspace(1e-20, 1e20, 40; endpoint=false);
    for i = 1:length(v)
        @test v[i] ≈ 1e-20 * 10.0^(i - 1)
    end
end

@testset "linspace" begin
    @test_throws ArgumentError linspace(1, 2, 1)
    @test_throws ArgumentError linspace(1.0, 2.0, -2)
    @test length(linspace(1, 2)) == 50
    v = linspace(-1, 1.0, 11; endpoint=true)
    for i = 1:length(v)
        @test isapprox(v[i], -1 + (i - 1) * 0.2, atol=eps(eltype(v)))
    end
    v = linspace(-1.0, 1, 10; endpoint=false)
    for i = 1:length(v)
        @test isapprox(v[i], -1 + (i - 1) * 0.2, atol=eps(eltype(v)))
    end
end

@testset "sqrtm1" begin
    function bigsqrtm1(x)
        x::BigFloat = x
        return Float64(sqrt(1 + x) - 1)
    end
    v = geomspace(1e-30, 10)
    @test bigsqrtm1.(v) ≈ sqrtm1.(v)
end
