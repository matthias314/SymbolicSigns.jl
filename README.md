# SymbolicSigns.jl

This package provides types to work with symbolic signs that often come up
in differential homological algebra and other graded contexts. Specifically,
given any commutative ring `R` and any type `T`, the package allows to define
a new commutative ring `SymbolicSignRing{T,R}` that extends `R` by expressions
of the form `(-1)^|x|` and `(-1)^(|x|*|y|)` where `|x|` and `|y|` represent
the symbolic degrees of `x::T` and `y::T`, respectively.

The symbolic signs `(-1)^(|x|*|y|)` and `(-1)^(|y|*|x|)` are treated as equal,
and `(-1)^(|x|*|x|)` is simplified to `(-1)^|x|`.

## Examples

Signs usually come up whenever two graded objects are swapped. Here is an example
with tensors, based on the package [LinearCombinations.jl](https://github.com/matthias314/LinearCombinations.jl).
We define the degree of a small letter to be its position in the alphabet.
Then we swap the factors of a
[`Tensor`](https://matthias314.github.io/LinearCombinations.jl/stable/tensor/#LinearCombinations.Tensor)
with the function
[`swap`](https://matthias314.github.io/LinearCombinations.jl/stable/tensor/#LinearCombinations.swap).
This incurs a sign which is the product of the degrees of the factors.
```julia
julia> using LinearCombinations

julia> LinearCombinations.deg(x::Char) = x-'a'+1

julia> deg('a'), deg('c')
(1, 3)

julia> t = Tensor('a', 'c')
a⊗c

julia> swap(t)
-c⊗a
```
We now do the same with symbolic signs. The computation of the sign out of the degrees
is done under the hood here, but it illustrates the functionality of the package.
Below we will construct signs explicitly.
```julia
julia> using LinearCombinations, SymbolicSigns

julia> LinearCombinations.deg(x::Char) = DegTerm(x)

julia> deg('a'), deg('c')
(|a|, |c|)

julia> t = Tensor('a', 'c')
a⊗c

julia> b = swap(t)
(-1)^(|a||c|)*c⊗a

julia> coefftype(b) == SymbolicSignRing{Char,Int}
true

julia> swap(b)
a⊗c
```
Similarly, signs come up when moving functions past tensor factors.
We define a linear function `f` of degree 1 and let `g` be the
tensor product of `f` with itself.
```julia
julia> @linear f; f(x::Char) = uppercase(x)
f (generic function with 2 methods)

julia> LinearCombinations.deg(::typeof(f)) = 1

julia> g = Tensor(f, f)
f⊗f

julia> g(t)
(-1)^(|a|)*A⊗C

julia> g(swap(t))
(-1)^(|a||c|+|c|)*C⊗A
```
We can also give `f` a symbolic degree.
```julia
julia> LinearCombinations.deg(::typeof(f)) = DegTerm('f')

julia> g(t)
(-1)^(|a||f|)*A⊗C
```
Here is an example showing how sign expressions are actually built.
```julia
julia> da = DegTerm('a')
|a|

julia> dc = DegTerm('c')
|c|

julia> e = dc + da*dc
|c|+|a||c|

julia> s1 = Sign(da)
(-1)^(|a|)

julia> s2 = 2*Sign(e)
2*(-1)^(|c|+|a||c|)

julia> s1-s2
-2*(-1)^(|c|+|a||c|)+(-1)^(|a|)

julia> s1*(s1-s2)
1-2*(-1)^(|a|+|c|+|a||c|)
```
