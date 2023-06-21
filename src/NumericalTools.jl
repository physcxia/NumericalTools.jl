module NumericalTools

export geomspace, linspace, sqrtm1, loginterpolator, Throw

using Interpolations: linear_interpolation, Throw

@doc raw"""
    geomspace(start, stop, num=50; endpoint=true)

Generate a Vector of geometric sequence,

```math
v_i = a \left(\frac{b}{a}\right)^{(i - 1)/n}
```

for ``i = 1, 2, ⋯, N``, where ``N = `` `num`, ``a = `` `start`, ``b = `` `stop`, and
``n = `` `num` ``- 1`` if `endpoint` = `true`, otherwise ``n = `` `num`.

# Arguments
- `start::Number`: The starting value of the sequence.
- `stop::Number`: The final value of the sequence.
- `num::Integer=50`: Number of samples to generate.

# Keywords
- `endpoint::Bool=true`: If true, `stop` is the last sample. Otherwise, it is not included.

# Example

```jldoctest
julia> geomspace(1, 1e4, 5)
5-element Vector{Float64}:
     1.0
    10.0
   100.0
  1000.0
 10000.0
```

"""
function geomspace(start::Number, stop::Number, num::Integer=50; endpoint=true)
    if (num <= 1) throw(ArgumentError("num <= 1")) end
    q = endpoint ? (stop/start)^(1/(num-1)) : (stop/start)^(1/num);
    res_type = typeof(q * oneunit(promote_type(typeof(start), typeof(stop))))
    res = Vector{res_type}(undef, num)
    res[1] = start;
    for i = 2:num
        @inbounds res[i] = res[i-1] * q;
    end
    return res
end


@doc raw"""
    linspace(start, stop, num=50; endpoint=true)

Generate a Vector of arithmetic sequence,

```math
v_i = a + (i - 1) \frac{b - a}{n}
```

for ``i = 1, 2, ⋯, N``, where ``N = `` `num`, ``a = `` `start`, ``b = `` `stop`, and
``n = `` `num` ``- 1`` if `endpoint` = `true`, otherwise ``n = `` `num`.

# Arguments
- `start::Number`: The starting value of the sequence.
- `stop::Number`: The final value of the sequence.
- `num::Integer=50`: Number of samples to generate.

# Keywords
- `endpoint::Bool=true`: If true, `stop` is the last sample. Otherwise, it is not included.

# Example
```jldoctest
julia> linspace(1.0, 5.0, 5)
5-element Vector{Float64}:
 1.0
 2.0
 3.0
 4.0
 5.0
```
"""
function linspace(start::Number, stop::Number, num::Integer=50; endpoint=true)
    if (num <= 1) throw(ArgumentError("num <= 1")) end
    d = endpoint ? (stop - start) / (num - 1) : (stop - start) / num;
    res = Vector{typeof(d)}(undef, num)
    res[1] = start;
    for i = 2:num
        @inbounds res[i] = res[i-1] + d;
    end
    return res
end


@doc raw"""
    sqrtm1(x)
    sqrtm1(x, a)

These functions calculate

```math
    \text{sqrtm1}(x) = \sqrt{1 + x} - 1,
```

for ``x > -1``, where ``x`` is dimensionless, and

```math
    \text{sqrtm1}(x, a) = \sqrt{a^2 + x} - a,
```

where ``x`` and ``a^2`` have the same dimension. They are designed to handle situations
where the absolute value of `x` is much smaller than 1 and `a`, where `a` is positive.

# Example
```jldoctest
julia> sqrt(1 + 1e-16) - 1
0.0

julia> sqrtm1(1e-16)
5.0e-17

julia> sqrt(1e16^2 + 1) - 1e16
0.0

julia> sqrtm1(1, 1e16)
5.0e-17

julia> m = 1; v = 1e-10; p = m*v; T = 1/2*m*v^2
5.0000000000000005e-21

julia> sqrt(p^2 + m^2) - m
0.0

julia> sqrtm1(p^2, m)
5.0000000000000005e-21
```

"""
sqrtm1(x::Integer) = sqrt(1 + x) - 1
function sqrtm1(x)
    if x < -1
        throw(DomainError("$x < -1"))
    end
    if abs(x) < 2eps(typeof(x))
        return x / 2
    end
    return sqrt(1 + x) - 1
end
function sqrtm1(x, a)
    if a > zero(a)
        return a * sqrtm1(x / a^2)
    end
    return sqrt(a^2 + x) - a
end


"""
    loginterpolator(x, y; method="loglog", extrapolation_bc=nothing)

Log-log and semi-log interpolator. Return an interpolator function.

Note that `y` is allowed to contain non-positive values in "loglog" and "ylog" modes, which
are replaced with `zero(eltype(y))`. Since log(0) = -Inf and the calculation involving `Inf`
may produce `NaN`, the interpolated function checks for and replaces `NaN` with zero.
See [this issue] (https://github.com/JuliaMath/Interpolations.jl/issues/459) for
interpolation with `Inf`. In general, non-positive values is invalid and should be avoid for
log interpolation, but this treatment makes life easier.

# Arguments
- `x::AbstractVector`: x array.
- `y::AbstractVector`: y array.

# Keywords
- `method="loglog"`: {"loglog", "xlog", "ylog"}, optional.
    Choose which axis to set in log scale.
- `extrapolation_bc=nothing`: Pass to [`linear_interpolation`]
    (http://juliamath.github.io/Interpolations.jl/latest/api/#Interpolations.linear_interpolation).

# Example
```jldoctest
julia> x = geomspace(1e-1, 10); y = x.^2; itp = loginterpolator(x, y);

julia> itp(1)
1.0

julia> itp(0.5)
0.25

julia> itp(5)
24.999999999999996
```

"""
function loginterpolator(
    x::AbstractVector,
    y::AbstractVector;
    method::String="loglog", extrapolation_bc=nothing
)::Function
    xunit = oneunit(eltype(x))
    yunit = oneunit(eltype(y))
    if extrapolation_bc isa Number
        extrapolation_bc = extrapolation_bc / yunit
    end

    if method == "loglog"
        y = @. ifelse(y < zero(y), zero(y), y)
        extrapolation = isnothing(extrapolation_bc) ? typemin(yunit)/yunit : extrapolation_bc
        loglog = linear_interpolation(
            log.(x / xunit), log.(y / yunit), extrapolation_bc=extrapolation)

        function wrapper_loglog(x)
            if x <= zero(x)
                @warn "x <= 0, zero returned" x
                return zero(yunit)
            end
            res = exp(loglog(log(x / xunit)))
            return isnan(res) ? zero(yunit) : res * yunit
        end

        return wrapper_loglog
    end

    if method == "ylog"
        y = @. ifelse(y < zero(y), zero(y), y)
        extrapolation = isnothing(extrapolation_bc) ? typemin(yunit)/yunit : extrapolation_bc
        ylog = linear_interpolation(x, log.(y ./ yunit), extrapolation_bc=extrapolation)

        function wrapper_ylog(x)
            res = exp(ylog(x))
            return (isnan(res) ? zero(yunit) : res * yunit)
        end

        return wrapper_ylog
    end

    if method == "xlog"
        extrapolation = isnothing(extrapolation_bc) ? zero(yunit)/yunit : extrapolation_bc
        xlog = linear_interpolation(
            log.(x ./ xunit), y ./ yunit, extrapolation_bc=extrapolation)

        function wrapper_xlog(x)
            if x <= zero(x)
                @warn "x <= 0, zero returned" x
                return zero(yunit)
            end
            return xlog(log(x / xunit)) * yunit
        end

        return wrapper_xlog
    end

    throw(ArgumentError("Unkown method: $method"))
end


end  # module NumericalTools
