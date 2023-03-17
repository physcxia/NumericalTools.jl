module NumericalTools

export geomspace, linspace, sqrtm1

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
    res = Vector{typeof(q)}(undef, num)
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

This function calculates

```math
    \text{sqrtm1}(x) = \sqrt{1 + x} - 1
```

for ``x > -1``, and takes care of small `x`.

# Example
```jldoctest
julia> sqrt(1 + 1e-16) - 1
0.0

julia> sqrtm1(1e-16)
5.0e-17
```
"""
function sqrtm1(x::AbstractFloat)
    if x < -1 throw(DomainError("$x < -1")) end
    if abs(x) < 2eps(typeof(x)) return x / 2 end
    return sqrt(1 + x) - 1
end
sqrtm1(x::Number) = sqrt(1 + x) - 1


end  # module NumericalTools
