"""
    $(@__MODULE__)

This package provides types to work with symbolic signs that often come up
in differential homological algebra and other graded contexts. Specifically,
given any commutative ring `R` and any type `T`, the package allows to define
a new commutative ring `WithSigns{T,R}` that extends `R` by expressions
of the form `(-1)^|x|` and `(-1)^(|x|*|y|)` where `|x|` and `|y|` represent
the symbolic degrees of `x::T` and `y::T`, respectively.

The symbolic signs `(-1)^(|x|*|y|)` and `(-1)^(|y|*|x|)` are treated as equal,
and `(-1)^(|x|*|x|)` is simplified to `(-1)^|x|`.

See also [`WithSigns`](@ref), [`Deg`](@ref).
"""
module SymbolicSigns

export Deg, DegSum, Sign, signexp, WithSigns

using StructEqualHash, LinearCombinations, Modulo2

import Base: show, ==, hash, one, isone, iszero, iseven, isodd, +, -, *, ^,
    convert, length, iterate, eltype, copy, promote_rule

using LinearCombinations: linear_convert
import LinearCombinations: Zero, termcoeff, withsign, sign_type, show_summand

#
# Deg
#

"""
    Deg{T}

    Deg{T}() where T
    Deg(x::T) where T
    Deg(x::T, y::T) where T

This type represents a *degree term* representing the degree of an element from `T`
or the product of the degrees of two elements from `T`. There is also a version
with no elements from `T` to represent the constant `1` as a degree term.

The degree terms `Deg(x, y)` and `Deg(y, x)` are treated as equal,
and `Deg(x, x)` is simplified to `Deg(x)`.

See also [`DegSum`](@ref).

# Examples
```jldoctest
julia> dx = Deg(:x)
|x|

julia> typeof(dx)
Deg{Symbol}

julia> dy = Deg(:y)
|y|

julia> dx*dy
|x||y|

julia> dx*dy == Deg(:x, :y)
true

julia> dx*dy == Deg(:y, :x)   # order doesn't matter
true

julia> dx*dx   # note the simplification
|x|

julia> dz = Deg(:z); dx*dy*dz
ERROR: cannot multiply two Deg's with more than 2 factors in total
[...]

julia> d1 = Deg{Symbol}()
1

julia> d1 == one(dx)
true
```
"""
struct Deg{T}
    n::Int8
    x::T
    y::T
    Deg{T}(x, y) where T = x == y ? new{T}(1, x) : new{T}(2, x, y)
    # do we need this test for equality in applications?
    Deg{T}(x) where T = new{T}(1, x)
    Deg{T}() where T = new{T}(0)
end

Deg(x::T...) where T = Deg{T}(x...)

function show(io::IO, t::Deg{T}) where T
    print(io, "Deg")
    if get(io, :typeinfo, Any) != Deg{T} &&
            !((t.n == 1 && typeof(t.x) == T) || (t.n == 2 && typeof(t.x) == T && typeof(t.y) == T))
        print(io, '{', T, '}')
    end
    # we cannot use join because it would print 'x' and :x as "x"
    # TODO: shall we order x and y?
    print(io, '(')
    t.n >= 1 && show(io, t.x)
    t.n == 2 && (print(io, ", "); show(io, t.y))
    print(io, ')')
end

function show(io::IO, ::MIME"text/plain", t::Deg)
    if t.n == 2
        xr = string(t.x)
        yr = string(t.y)
        if xr <= yr
            print(io, '|', xr, "||", yr, '|')
        else
            print(io, '|', yr, "||", xr, '|')
        end
    elseif t.n == 1
        print(io, '|', t.x, '|')
    else
        print(io, '1')
    end
end

function convert(::Type{D}, t::Deg) where D <: Deg
    if t isa D
        t
    elseif t.n == 0
        D()
    elseif t.n == 1
        D(t.x)
    else
        D(t.x, t.y)
    end
end

length(t::Deg)::Int = t.n

eltype(t::Deg{T}) where T = T

function iterate(t::Deg, i = 1)
    if i > length(t)
        nothing
    elseif i == 1
        t.x, 2
    elseif i == 2
        t.y, 3
    end
end

function ==(t::Deg{T}, u::Deg{T}) where T
    if t.n != u.n
        false
    elseif t.n == 2
        (t.x == u.x && t.y == u.y) || (t.x == u.y && t.y == u.x)
    elseif t.n == 1
        t.x == u.x
    else
        true
    end
end

function hash(t::Deg, h::UInt)
    h = StructEqualHash.typehash(Deg, h)
    if t.n == 2
        h ⊻ hash(t.x) ⊻ hash(t.y)
    elseif t.n == 1
        h ⊻ hash(t.x)
    else
        h
    end
end

isone(t::Deg) = iszero(length(t))

one(::Type{T}) where T <: Deg = T()
one(::T) where T <: Deg = one(T)

function *(t::Deg{T}, u::Deg{T}) where T
    if t.n == 1 && u.n == 1
        Deg(t.x, u.x)
    elseif t.n == 0
        u
    elseif u.n == 0
        t
    else
        error("cannot multiply two Deg's with more than 2 factors in total")
    end
end

@linear_broadcastable Deg

#
# DegSum
#

"""
    DegSum{T} == Linear{Deg{T}, ZZ2}

A *sign exponent* is a linear combination of degree terms with coefficients in `ZZ2`.

See also [`LinearCombinations.Linear`](https://matthias314.github.io/LinearCombinations.jl/stable/linear/#LinearCombinations.Linear),
`Modulo2.ZZ2`, [`Deg`](@ref), [`Sign`](@ref).

# Examples
```jldoctest
julia> se = dx + dx*dy + 1
|x||y|+|x|+1

julia> typeof(se) == DegSum{Symbol}
true

julia> se == DegSum(dx => 1, dx*dy => -1, Deg{Symbol}() => -1)
true

julia> DegSum(dx) == DegSum(dx => 1)   # note the short form
true
```
"""
const DegSum{T} = Linear{Deg{T}, ZZ2}

DegSum(t::Pair{Deg{T}}...) where T = DegSum{T}(t)
DegSum{T}(n::Number) where T = Linear(one(Deg{T}) => ZZ2(n))

+(ts::Deg{T}...) where T = DegSum{T}(map(t -> t => one(ZZ2), ts))
+(t::Deg, n::Number) = isone(t) ? DegSum(t => ZZ2(n+1)) : DegSum(t => one(ZZ2), one(t) => ZZ2(n))

+(e::DegSum{T}, n::Number) where T = addmul(e, one(Deg{T}), ZZ2(n))

+(n::Number, t::Deg) = t+n
+(n::Number, e::DegSum) = e+n

-(t::Deg) = +t
-(t::Deg{T}, u::Deg{T}) where T = t+u
-(t::Deg, n::Number) = t+n
-(e::DegSum, n::Number) = e+n
-(n::Number, t::Deg) = t+n
-(n::Number, e::DegSum) = e+n

*(::Deg, ::Zero) = Zero()
*(::Zero, ::Deg) = Zero()
*(::DegSum, ::Zero) = Zero()
*(::Zero, ::DegSum) = Zero()
*(c::Number, t::Deg) = DegSum(t => c)
*(t::Deg, c::Number) = c*t

withsign(e::Union{DegSum,Deg}, c) = has_char2(c) ? c : Sign(e)*c

convert(::Type{DegSum}, t::Deg) = +t
convert(::Type{DegSum{T}}, t::Deg) where T = +convert(Deg{T}, t)
convert(::Type{DegSum{T}}, n::Number) where T = DegSum{T}(n)

#
# Sign
#

"""
    Sign{T}

This is a wrapper around `DegSum{T}` that allows to use a multiplicative notation.

See also [`DegSum`](@ref), [`WithSigns`](@ref).

# Examples
```jldoctest
julia> s1 = Sign(dx + dx*dy)
(-1)^(|x||y|+|x|)

julia> s2 = Sign(dy + 1)
(-1)^(|y|+1)

julia> s1*s2
(-1)^(|y|+|x||y|+|x|+1)

julia> isone(s1^2)
true

julia> convert(ZZ2, s1)
1
```
"""
struct Sign{T}
    e::DegSum{T}
end

Sign(t::Deg) = Sign(convert(DegSum, t))

function show(io::IO, s::Sign{T}) where T
    print(io, "Sign(")
    show(io, s.e)
    print(io, ')')
end

function show(io::IO, ::MIME"text/plain", s::Sign)
    if isone(s)
        print(io, '1')
    else
        e = signexp(s)
        print(io, "(-1)^")
        length(e) > 1 && print(io, '(')
        show(IOContext(io, :compact => true), MIME"text/plain"(), e)
        length(e) > 1 && print(io, ')')
    end
    nothing
end

@struct_equal_hash Sign{T} where T

signexp(s::Sign) = s.e

convert(::Type{S}, s::Sign) where S <: Sign = S(signexp(s))

copy(s::Sign) = Sign(copy(s.e))

length(s::Sign) = length(signexp(s))

eltype(s::Sign) = eltype(signexp(s))

iterate(s::Sign, state...) = iterate(signexp(s), state...)

iszero(::Sign) = false

-(s::Sign{T}) where T = WithSigns{T}(s => -1)

+(ss::Sign{T}...) where T = WithSigns{T}(map(s -> s => 1, ss))
-(s::Sign{T}, t::Sign{T}) where T = WithSigns{T}(s => 1, t => -1)

+(c, s::Sign{T}) where T = Linear(s => 1, one(Sign{T}) => c)
+(s::Sign, c) = c+s

-(c, s::Sign{T}) where T = Linear(s => -1, one(Sign{T}) => c)
-(s::Sign{T}, c) where T = Linear(s => 1, one(Sign{T}) => -c)

isone(s::Sign) = iszero(signexp(s))

one(::Type{Sign{T}}) where T = Sign(zero(DegSum{T}))
one(::T) where T <: Sign = one(T)

*(s1::Sign{T}, s2::Sign{T}) where T = Sign(signexp(s1) + signexp(s2))

function mul!(s1::Sign{T}, s2::Sign{T}) where T
    add!(signexp(s1), signexp(s2))
    s1
end

^(s::Sign, n::Integer) = isodd(n) ? s : one(s)

Modulo2.ZZ2(::Sign) = one(ZZ2)

for op in (:(+), :(-), :(*))
    @eval $op(c::ZZ2, s::Sign) = $op(c, ZZ2(s))
    @eval $op(s::Sign, c::ZZ2) = $op(ZZ2(s), c)
end

iseven(::Sign) = true
isodd(::Sign) = false

@linear_broadcastable Sign

function show_summand(io::IO, s::Sign, cs)
    show(io, MIME"text/plain"(), cs)
    isone(s) && return
    print(io, '*')
    show(io, MIME"text/plain"(), s)
end

#
# WithSigns
#

"""
    LinearCombinations.termcoeff((s, c)::Pair{Sign{T}}) where T -> Pair{Sign{T}}

Convert the term-coefficient pair `(s, c)` to one where the exponent of the sign
has no constant term. If the exponent of `s` has a constant term, then it is
removed and instead the sign of `c` changed.

This is important to make linear combinations of signs unique.

See also [`LinearCombinations.termcoeff`](https://matthias314.github.io/LinearCombinations.jl/stable/extensions/#LinearCombinations.termcoeff).

# Examples
```jldoctest
julia> using LinearCombinations: termcoeff

julia> dx = Deg(:x)
|x|

julia> termcoeff(Sign(dx) => 2)
(-1)^(|x|) => 2

julia> termcoeff(Sign(dx+1) => 2)
(-1)^(|x|) => -2
```
"""
function termcoeff(sc::Pair{Sign{T}}) where T
    s, c = sc
    x = Deg{T}()
    if isone(s.e[x])
        Sign(addmul(s.e, x, one(ZZ2))) => -c
    else
        sc
    end
end

"""
    WithSigns{T,R} = Linear{Sign{T},R}

This type extends the commutative ring `R` by elements of type `Sign{T}`.
Elements are formal linear combinations with terms of type `Sign{T}` and
coefficients of type `R`.

See also [`LinearCombinations.Linear`](https://matthias314.github.io/LinearCombinations.jl/stable/linear/#LinearCombinations.Linear),
[`Sign`](@ref).

# Examples
```jldoctest
julia> dx, dy = Deg(:x), Deg(:y)
(|x|, |y|)

julia> s1, s2 = Sign(dx + dx*dy), Sign(dy + 1)
((-1)^(|x||y|+|x|), (-1)^(|y|+1))

julia> a = Linear(s1 => 2, s2 => -1)
(-1)^(|y|)+2*(-1)^(|x||y|+|x|)

julia> convert(ZZ2, a)
1

julia> 3.0*a
6.0*(-1)^(|x|+|x||y|)+3.0*(-1)^(|y|)

julia> ZZ2(1)*a
1
```
Note that the extended ring has zero divisors:
```jldoctest
julia> a = Linear(Sign(dx) => 1, Sign(dy) => 1)
(-1)^(|y|)+(-1)^(|x|)

julia> b = Linear(Sign(dx) => 1, Sign(dy) => -1)
-(-1)^(|y|)+(-1)^(|x|)

julia> a*b
0
```
"""
const WithSigns{T,R} = Linear{Sign{T},R}

# needed because "WithSigns" means "WithSigns{T,R} where {T,R}"
WithSigns(x::Pair{S}...) where S <: Sign = Linear{S}(x...)

Modulo2.ZZ2(a::WithSigns) = sum(ZZ2, coeffs(a); init = zero(ZZ2))

iseven(a::WithSigns) = isone(ZZ2(a))
isodd(a::WithSigns) = iszero(ZZ2(a))

*(s::Sign{T}, c::R) where {T,R<:Number} = WithSigns{T,R}(s => c)
*(c::R, s::Sign{T}) where {T,R<:Number} = s*c

*(a::WithSigns{T}, s::Sign{T}) where T = LinearCombinations.mul(a, s)

for op in (:(+), :(-), :(*))
    @eval $op(a::WithSigns, b::WithSigns) = invoke($op, Tuple{AbstractLinear, AbstractLinear}, a, b)
    if op in (:(+), :(-))
        @eval $op(a::WithSigns{T,R}, s::Sign{T}) where {T,R} = addmul(a, s, $op(one(R)))
        @eval $op(s::Sign{T}, a::WithSigns{T}) where T = add!($op(a), s)
    end
    @eval function $op(a::WithSigns{T}, c) where T
        if has_char2(c)
            $op(ZZ2(a), c)
        elseif $op in (+, -) && c isa Number
            addmul(a, one(Sign{T}), $op(c))
        else
            R = coefftype(a)
            S = typeof(c)
            invoke($op, Tuple{AbstractLinear{Sign{T},R}, S}, a, c)
        end
    end
end

*(a::WithSigns, s::LinearCombinations.Sign) = isone(s) ? a : -a

function -(c, a::WithSigns{T,R}) where {T,R}
    if has_char2(c)
        ZZ2(a) + c
    elseif c isa Number
        addmul!((-one(R)*one(c))*a, one(Sign{T}), c)
    else
        S = typeof(c)
        invoke(-, Tuple{S, AbstractLinear{Sign{T},R}}, c, a)
    end
end

convert(::Type{SR}, a::WithSigns{T}) where {T, SR<:WithSigns{T}} = linear_convert(SR, a)
convert(::Type{SR}, s::Sign{T}) where {T, R, SR<:WithSigns{T,R}} = SR(s => one(R))
convert(::Type{SR}, c::Number) where {T, SR<:WithSigns{T}} = SR(one(Sign{T}) => c)
convert(::Type{SR}, s::LinearCombinations.Sign) where {T, R, SR<:WithSigns{T,R}} = SR(one(Sign{T}) => convert(R, s))

function promote_rule(::Type{WithSigns{T,R}}, ::Type{S}) where {T,R,S}
    U = promote_type(R, S)
    has_char2(S) ? U : WithSigns{T,U}
end

promote_rule(::Type{WithSigns{T,R}}, ::Type{Sign{T}}) where {T,R} = WithSigns{T,R}
promote_rule(::Type{WithSigns{T,R}}, ::Type{WithSigns{T,S}}) where {T,R,S} = WithSigns{T,promote_type(R,S)}

# it's important that Sign{T} is again the second argument
promote_rule(::Type{R}, ::Type{Sign{T}}) where {T,R} = has_char2(R) ? R : WithSigns{T,R}

sign_type(::Type{DegSum{T}}) where T = WithSigns{T,LinearCombinations.Sign}
sign_type(::Type{Deg{T}}) where T = WithSigns{T,LinearCombinations.Sign}

end # module
