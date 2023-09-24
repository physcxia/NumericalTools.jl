# NumericalTools

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://physcxia.github.io/NumericalTools.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://physcxia.github.io/NumericalTools.jl/dev/)
[![Build Status](https://github.com/physcxia/NumericalTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/physcxia/NumericalTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/physcxia/NumericalTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/physcxia/NumericalTools.jl)

A simple Julia package offering some useful numerical tools for daily scientific research.

Currently implemented tools are:

- `geomspace`: Generate a geometric sequence like [numpy.geomspace](https://numpy.org/doc/stable/reference/generated/numpy.geomspace.html).
- `linspace`: Generate an arithmetic sequence like [numpy.linspace](https://numpy.org/doc/stable/reference/generated/numpy.linspace.html#numpy.linspace).
- `logspace`: Generate a log-linear sequence like [numpy.logspace](https://numpy.org/doc/stable/reference/generated/numpy.logspace.html#numpy.logspace).
- `sqrtm1`: Accurately calculate $f(x) = \sqrt{1 + x} - 1$ for small $x$.
- `loginterpolator`: A wrapper of [`linear_interpolation`](https://juliamath.github.io/Interpolations.jl/latest/convenience-construction/) for 1-D linear interpolation in log scale.

See [documentation](https://physcxia.github.io/NumericalTools.jl/dev/) for details.

## Installation

```julia
julia> ]
pkg> add NumericalTools
```
