using NumericalTools
using Test

@testset "geomspace" begin
    @test_throws ArgumentError geomspace(1, 2, 1)
    @test_throws ArgumentError geomspace(1.0, 2.0, -2)
    @test length(geomspace(1, 2)) == 50

    for v in [
        geomspace(1e-20, 1e20, 41),
        geomspace(1e-20, 1e20, 40; endpoint=false),
        geomspace(-1e-10, -1e10, 21),
    ]
        @test v ≈ [v[1] * 10.0^(i - 1) for i = 1:length(v)]
    end
    for v in [
        geomspace(1e20, 1e-20, 41),
        geomspace(1e20, 1e-20, 40; endpoint=false),
        geomspace(-1e10, -1e-10, 21),
    ]
        @test v ≈ [v[1] * 10.0^(1 - i) for i = 1:length(v)]
    end
end

@testset "linspace" begin
    @test_throws ArgumentError linspace(1, 2, 1)
    @test_throws ArgumentError linspace(1.0, 2.0, -2)
    @test length(linspace(1, 2)) == 50

    for v in [
        linspace(-1, 1.0, 11; endpoint=true),
        linspace(-1.0, 1, 10; endpoint=false),
    ]
        @test isapprox(v, [-1 + (i - 1) * 0.2 for i = 1:length(v)], atol=2eps(eltype(v)))
    end
    for v in [
        linspace(1.0, -1, 11; endpoint=true),
        linspace(1, -1.0, 10; endpoint=false),
    ]
        @test isapprox(v, [1 - (i - 1) * 0.2 for i = 1:length(v)], atol=2eps(eltype(v)))
    end
end

@testset "sqrtm1" begin
    @test_throws DomainError sqrtm1(-2.0)
    function bigsqrtm1(x)
        x::BigFloat = x
        return Float64(sqrt(1 + x) - 1)
    end

    vplus = geomspace(1e-30, 10)
    @test bigsqrtm1.(vplus) ≈ sqrtm1.(vplus)
    @test isapprox(bigsqrtm1.(vplus), sqrtm1.(vplus), rtol=eps())

    vminus = geomspace(-1, -1e-30)
    @test bigsqrtm1.(vminus) ≈ sqrtm1.(vminus)
    @test isapprox(bigsqrtm1.(vminus), sqrtm1.(vminus), rtol=eps())

    @test sqrtm1(-1) == -1
    @test sqrtm1(0) == 0
end
